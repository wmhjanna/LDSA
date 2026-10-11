#!/usr/bin/env python3
"""Joint moving-block bootstrap for exact-overlap LDSA-pollutant trends."""
import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import pandas as pd
from build_engine import library_path

engine_path = library_path()
if not engine_path.exists():
    raise SystemExit("Build the C++ engine first: python scripts/build_engine.py")
engine = ctypes.CDLL(str(engine_path)).weighted_sen
I = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
D = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
engine.argtypes = [ctypes.c_int]*3 + [I, I, D, D, D, I, I, I, D]
engine.restype = None

def sen_draws(t, y, counts):
    n = len(t)
    ii, jj = np.triu_indices(n, 1)
    slopes = (y[jj]-y[ii])/(t[jj]-t[ii])
    order = np.argsort(slopes, kind='stable')
    out = np.empty((len(counts),3), dtype=np.float64)
    engine(n,len(ii),len(counts),ii[order].astype('int32'),jj[order].astype('int32'),
           np.ascontiguousarray(slopes[order]),np.ascontiguousarray(t,dtype='float64'),np.ascontiguousarray(y,dtype='float64'),
           np.arange(n,dtype='int32'),np.argsort(y,kind='stable').astype('int32'),np.ascontiguousarray(counts,dtype='int32'),out)
    return out

def mbb_counts(n, length, B, rng):
    starts = rng.integers(0,n-length+1,(B,int(np.ceil(n/length))))
    indices = (starts[:,:,None] + np.arange(length)).reshape(B,-1)[:,:n]
    out = np.zeros((B,n),dtype=np.int32)
    np.add.at(out,(np.repeat(np.arange(B),n),indices.ravel()),1)
    return out

def calendar_counts(months,length,B,rng):
    # Sample original time-stamped bivariate records in fixed calendar windows.
    # Empty windows contain no observations; no missing values are imputed.
    windows = [np.flatnonzero((months>=s)&(months<s+length)) for s in range(int(months[-1])-length+2)]
    windows = [w for w in windows if len(w)]
    n=len(months); out=np.zeros((B,n),dtype=np.int32)
    for b in range(B):
        blocks=[];total=0
        while total<n:
            z=windows[rng.integers(len(windows))];blocks.append(z);total+=len(z)
        out[b]=np.bincount(np.concatenate(blocks)[:n],minlength=n)
    return out

def bh(p, by=False):
    p=np.asarray(p,float);order=np.argsort(p);n=len(p)
    a=p[order]*n/np.arange(1,n+1)*(sum(1/np.arange(1,n+1)) if by else 1)
    q=np.empty(n);q[order]=np.minimum(1,np.minimum.accumulate(a[::-1])[::-1]);return q

def stats(draws, point, starts):
    valid=np.isfinite(draws).all(1)
    delta=draws[valid,0]-draws[valid,1];B=len(delta)
    lo,hi=np.quantile(delta,[.025,.975])
    # Equal-tail percentile-bootstrap probability; approximate, not an exact permutation p value.
    p=min(1.,2*min((1+np.sum(delta<=0))/(B+1),(1+np.sum(delta>=0))/(B+1)))
    baseline_low=np.quantile(starts,[.025],axis=0)[0]
    reliable=bool(np.all(baseline_low>0) and valid.all())
    cov=np.cov(draws[valid].T,ddof=1)
    np.testing.assert_allclose(np.var(delta,ddof=1),cov[0,0]+cov[1,1]-2*cov[0,1],rtol=1e-10,atol=1e-10)
    return dict(Delta=point,CI_lower=lo,CI_upper=hi,Delta_SE=np.std(delta,ddof=1),
                Bootstrap_bias=np.mean(delta)-point,Approx_p=p,Inference_available=reliable,
                CI_excludes_zero=bool((lo>0 or hi<0) and reliable),
                LDSA_baseline_p025=baseline_low[0],Pollutant_baseline_p025=baseline_low[1],
                Nonpositive_baseline_draws=int((starts<=0).any(1).sum()),Finite_draws=B,
                Bootstrap_slope_correlation=float(np.corrcoef(draws[valid].T)[0,1]),
                Unpaired_SE_diagnostic=float(np.sqrt(cov[0,0]+cov[1,1])))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Prepared wide monthly CSV; see data/README.md")
    parser.add_argument("output", type=Path)
    parser.add_argument("--replicates", type=int, default=9999)
    parser.add_argument("--all-settings", action="store_true", help="Also run blocks 3/6/12 and calendar windows of 6 months")
    args = parser.parse_args()
    if args.replicates < 100:
        parser.error("At least 100 replicates are required.")
    source = pd.read_csv(args.input, float_precision="round_trip")
    required = {"Site", "Station_class", "Pollutant", "Date", "LDSA", "Pollutant_value"}
    if not required.issubset(source):
        parser.error("Missing columns: " + ", ".join(sorted(required - set(source))))
    source["Date"] = pd.to_datetime(source["Date"], errors="raise")
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "bootstrap_draws").mkdir(exist_ok=True)
    results, eligibility, audit = [], [], []
    started = time.time()
    for (metric, site), g in source.groupby(["Pollutant", "Site"], sort=True):
        g = g.sort_values("Date")
        dates = pd.DatetimeIndex(g.Date)
        if dates.hasnans or dates.to_period("M").duplicated().any():
            raise ValueError(f"{site}/{metric}: one valid record per common month is required")
        if g.Station_class.nunique() != 1:
            raise ValueError(f"{site}/{metric}: inconsistent station class")
        Y = g[["LDSA", "Pollutant_value"]].to_numpy(float)
        if not np.isfinite(Y).all():
            raise ValueError(f"{site}/{metric}: remove missing records jointly before fitting")
        n = len(g)
        months = np.asarray((dates.year - dates[0].year)*12 + dates.month - dates[0].month)
        coverage = n / (int(months[-1]) + 1)
        included = bool(n >= 60 and coverage >= .60)
        eligibility.append(dict(Site=site, Pollutant=metric, N_months=n, Coverage=coverage, Included=included))
        if not included:
            continue
        t = (dates - dates[0]).days.to_numpy(float)
        fit = np.stack([sen_draws(t, Y[:, j], np.ones((1, n), dtype=np.int32))[0] for j in range(2)])
        point = fit[0, 2] - fit[1, 2]
        default = int(np.rint(n**(1/3)))
        base = dict(Site=site, Station_class=g.Station_class.iloc[0], Pollutant=metric, N_months=n,
                    First_month=str(dates[0].date()), Last_month=str(dates[-1].date()), Coverage=coverage,
                    Missing_calendar_months=int(months[-1]+1-n), Gap_transitions=int(np.sum(np.diff(months)>1)),
                    LDSA_slope=fit[0, 2], Pollutant_slope=fit[1, 2],
                    LDSA_start_concentration=fit[0, 1], Pollutant_start_concentration=fit[1, 1])
        settings = [("Auto", default)]
        if args.all_settings:
            settings += [("Block3", 3), ("Block6", 6), ("Block12", 12), ("Calendar6", 6)]
        draws_save, index_cache = {}, {}
        for method, length in settings:
            seed = 20261008 + int(hashlib.sha256(f"{site}|{metric}|{method}".encode()).hexdigest()[:8], 16)
            rng = np.random.default_rng(seed)
            if method == "Calendar6":
                counts = calendar_counts(months, length, args.replicates, rng)
            elif length in index_cache:
                counts = index_cache[length]
            else:
                counts = mbb_counts(n, length, args.replicates, rng)
                index_cache[length] = counts
            assert (counts.sum(1) == n).all()
            # The same row multiplicities jointly resample the two indicators.
            z = np.stack([sen_draws(t, Y[:, j], counts) for j in range(2)], axis=1)
            draws, starts = z[:, :, 2], z[:, :, 1]
            if np.isfinite(draws).all(axis=1).sum() < 2:
                raise ValueError(f"{site}/{metric}: fewer than two finite joint bootstrap draws")
            result = dict(**base, Method=method, Block_length=length,
                          Block_unit="calendar months" if method == "Calendar6" else "retained observations",
                          Replicates=args.replicates, **stats(draws, point, starts))
            if np.any(fit[:, 1] <= 0):
                result["Inference_available"] = False
                result["CI_excludes_zero"] = False
            result["Status"] = "Available" if result["Inference_available"] else "Percentage normalization unstable"
            results.append(result)
            draws_save[method + "_slopes"] = draws
            draws_save[method + "_starts"] = starts
            if method == "Auto":
                for b in range(3):
                    idx = np.repeat(np.arange(n), counts[b])
                    for j, name in enumerate(["LDSA", metric]):
                        audit.append(dict(Site=site, Pollutant=metric, Metric=name, Replicate=b,
                            Slope_per_day=z[b,j,0], Start=z[b,j,1], Percent=z[b,j,2],
                            t=t[idx].tolist(), y=Y[idx,j].tolist()))
        safe_id = lambda s: "".join(c if c.isalnum() or c in "_.-" else "_" for c in str(s))
        np.savez_compressed(args.output / "bootstrap_draws" / f"{safe_id(site)}__{safe_id(metric)}.npz", **draws_save)
        print(site, metric, n, "months; elapsed", round(time.time()-started, 1), "s", flush=True)
    pd.DataFrame(eligibility).to_csv(args.output / "eligibility.csv", index=False)
    if not results:
        raise SystemExit("No pair meets >=60 common months and >=60% calendar coverage.")
    result = pd.DataFrame(results)
    for method, ids in result.groupby("Method").groups.items():
        p = result.loc[ids, "Approx_p"].where(result.loc[ids, "Inference_available"], 1.)
        result.loc[ids, "BH_q_input_family"] = bh(p)
        result.loc[ids, "BY_q_input_family"] = bh(p, by=True)
        result.loc[ids, "Family_size"] = len(ids)
    result.to_csv(args.output / "paired_difference_results.csv", index=False)
    (args.output / "independent_reference_samples.json").write_text(json.dumps(audit))
    (args.output / "run_metadata.json").write_text(json.dumps(dict(
        replicates=args.replicates, minimum_common_months=60, minimum_coverage=.60,
        numpy=np.__version__, pandas=pd.__version__, elapsed_seconds=time.time()-started,
        multiplicity="All eligible input pairs, separately for each block setting"), indent=2))

if __name__ == "__main__":
    main()
