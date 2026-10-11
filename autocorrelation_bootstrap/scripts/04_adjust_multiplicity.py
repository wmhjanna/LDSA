#!/usr/bin/env python3
"""Recompute BH and BY across all supplied eligible pairs, per block setting."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def adjust(p, by=False):
    p = np.asarray(p, float)
    order = np.argsort(p)
    n = len(p)
    scaled = p[order] * n / np.arange(1, n+1)
    if by:
        scaled *= np.sum(1 / np.arange(1, n+1))
    q = np.empty(n)
    q[order] = np.minimum(1, np.minimum.accumulate(scaled[::-1])[::-1])
    return q

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--expected-family-size", type=int)
    args = parser.parse_args()
    d = pd.read_csv(args.input, float_precision="round_trip")
    required = {"Method", "Site", "Pollutant", "N_months", "Coverage", "Approx_p", "Inference_available"}
    if not required.issubset(d):
        parser.error("Missing required paired-result columns")
    d = d[d.N_months.ge(60) & d.Coverage.ge(.60)].copy()
    if d.empty or d.duplicated(["Method", "Site", "Pollutant"]).any():
        parser.error("Empty or duplicated pair set")
    for _, ids in d.groupby("Method").groups.items():
        if args.expected_family_size is not None and len(ids) != args.expected_family_size:
            parser.error(f"Expected {args.expected_family_size} pairs, found {len(ids)}")
        available = d.loc[ids, "Inference_available"].astype(str).str.lower().eq("true")
        # Keep unstable pairs in the family with conservative p = 1.
        p = d.loc[ids, "Approx_p"].where(available, 1.).to_numpy(float)
        if not (np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()):
            parser.error("Invalid bootstrap probabilities")
        d.loc[ids, "BH_q_recomputed"] = adjust(p)
        d.loc[ids, "BY_q_recomputed"] = adjust(p, by=True)
        d.loc[ids, "Family_size"] = len(ids)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(args.output, index=False)
    print("Wrote", len(d), "pair-setting rows to", args.output)

if __name__ == "__main__":
    main()
