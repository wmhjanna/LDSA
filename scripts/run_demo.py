#!/usr/bin/env python3
"""Recalculate the observed BCN_UB record and verify archived study results."""
from pathlib import Path
import argparse, json, platform, sys, time
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from integration import integrate_ldsa
from coverage import compare_coverage
from temporal import paired_monthly, observation_breakpoint

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, default=ROOT/'demo_output')
    args = ap.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    folder = ROOT/'examples/data/BCN_UB'
    pnsd = pd.read_csv(folder/'pnsd_bin_concentrations.csv.gz')
    df = pd.read_csv(folder/'deposition_fractions.csv')
    diameters = df.diameter_nm.to_numpy()
    columns = [f'D_{d:g}nm_cm-3' for d in diameters]
    values = pnsd[columns].to_numpy(float)
    result = pd.DataFrame({'Time': pd.to_datetime(pnsd.Time, format='mixed')})
    for region in ['HA','TB','ALV','Total']:
        result[region] = integrate_ldsa(values, diameters, df[region].to_numpy())
    expected = pd.read_csv(ROOT/'examples/expected/BCN_UB_ldsa.csv.gz')
    pd.testing.assert_series_equal(result.Time, pd.to_datetime(expected.Time, format='mixed'), check_names=False)
    np.testing.assert_allclose(result.iloc[:,1:], expected.iloc[:,1:], rtol=1e-12, atol=1e-9, equal_nan=True)
    result.to_csv(args.output/'BCN_UB_ldsa.csv.gz', index=False, compression={'method':'gzip','mtime':0})
    monthly = result.set_index('Time').resample('MS').mean()
    monthly.to_csv(args.output/'BCN_UB_monthly.csv')
    paired = paired_monthly(result.Time,result.Total)
    breakpoint = observation_breakpoint(result.Time,result.Total,len(paired))
    assert len(paired)==75 and breakpoint==pd.Timestamp('2019-10-05 03:30:00')
    (args.output/'BCN_UB_breakpoint.json').write_text(json.dumps({'valid_paired_months':len(paired),'breakpoint_timestamp':str(breakpoint),'display_month':str(breakpoint.to_period('M'))},indent=2))
    v = result.loc[result.Time.dt.year.between(2015,2019),'Total'].dropna()
    spatial = pd.DataFrame([dict(Site='BCN_UB',Type='UB',Period='2015-2019',n=len(v),Mean=v.mean(),Median=v.median(),P25=v.quantile(.25),P75=v.quantile(.75),P95=v.quantile(.95))])
    pd.testing.assert_frame_equal(spatial,pd.read_csv(ROOT/'examples/expected/BCN_UB_spatial_statistics.csv'),check_dtype=False,rtol=1e-12,atol=1e-9)
    spatial.to_csv(args.output/'BCN_UB_spatial_statistics.csv',index=False)
    # Supplementary network-level example uses site summaries, not other stations' PNSD.
    results = compare_coverage(pd.read_csv(ROOT/'examples/data/site_statistics.csv'))
    for table, name in zip(results,['class_summary.csv','network_summary.csv','pairwise_comparisons.csv']):
        pd.testing.assert_frame_equal(table,pd.read_csv(ROOT/'examples/expected'/name),check_dtype=False,rtol=1e-10,atol=1e-10)
        table.to_csv(args.output/name,index=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for font_dir in [Path('/System/Library/Fonts/Supplemental'),Path('C:/Windows/Fonts')]:
        if font_dir.exists():
            for font in font_dir.glob('*.ttf'):
                if font.name.lower().startswith(('times new roman','times')):
                    font_manager.fontManager.addfont(str(font))
    family = 'Times New Roman' if any(f.name == 'Times New Roman' for f in font_manager.fontManager.ttflist) else 'DejaVu Serif'
    plt.rcParams.update({'font.family':family,'font.weight':'bold','axes.labelweight':'bold','axes.titleweight':'bold','pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(7.2,3.2))
    ax.plot(monthly.index,monthly.Total,color='#1487c9',lw=1.4)
    ax.set(xlabel='Year',ylabel=r'LDSA ($\mathrm{\mu m^{2}\ cm^{-3}}$)',title='BCN_UB: monthly mean total LDSA')
    ax.margins(x=.01);ax.grid(axis='y',alpha=.22);fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(args.output/f'BCN_UB_demo.{ext}',dpi=300)
    plt.close(fig)
    report={'status':'PASS','station':'BCN_UB','rows':len(result),'channels':len(diameters),'valid_LDSA_rows':int(result.Total.notna().sum()),'maximum_absolute_difference_um2_cm-3':float(np.nanmax(np.abs(result.iloc[:,1:].to_numpy()-expected.iloc[:,1:].to_numpy()))),'checked':['All BCN_UB timestamps and four LDSA output columns','BCN_UB 2015-2019 spatial statistics','BCN_UB change point and 75 paired-comparison months','65/53/49 coverage summaries and pairwise comparisons'],'seconds':round(time.perf_counter()-start,3),'python':platform.python_version(),'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'plot_font':family}
    (args.output/'validation.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__ == '__main__':main()
