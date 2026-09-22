from pathlib import Path
import json,itertools,hashlib,sys
sys.path.insert(0,'@PROJECT_ROOT@/work/ldsa_support/deps')
import pandas as pd,numpy as np
from scipy.stats import kruskal,mannwhitneyu
ROOT=Path('@PROJECT_ROOT@');W=ROOT/'work/record_exclusions/spatial'
O=Path('@DELIVERY_ROOT@/record_exclusion_analysis/diameter_coverage')
src=ROOT/'outputs/record_exclusions/analysis/site_spatial_statistics.csv';d=pd.read_csv(src)
coverage=pd.read_csv(O/'coverage_selection.csv');configs=pd.read_csv(O/'configuration_coverage.csv')
upper_excluded=set(configs.loc[configs.native_high_nm<500-1e-9,'Site']);keep=set(coverage.loc[coverage.keep,'Site'])
sets={'Main_65':d,'Upper_complete_53':d[~d.Site.isin(upper_excluded)],'Strict_49':d[d.Site.isin(keep)]}
assert [len(z) for z in sets.values()]==[65,53,49]
order=['TR','UB','SUB','RB'];summ=[];overall=[];pairs=[]
for scenario,z in sets.items():
 groups=[z.loc[z.Type.eq(t),'Mean'].to_numpy() for t in order];h,p=kruskal(*groups)
 overall.append(dict(scenario=scenario,n=len(z),mean=float(z.Mean.mean()),median=float(z.Mean.median()),H=float(h),df=3,p=float(p),epsilon_squared=float((h-3)/(len(z)-4))))
 for typ,v in zip(order,groups):summ.append(dict(scenario=scenario,Type=typ,n=len(v),mean=float(v.mean()),median=float(np.median(v)),P25=float(np.quantile(v,.25)),P75=float(np.quantile(v,.75))))
 local=[]
 for a,b in itertools.combinations(order,2):
  x=z.loc[z.Type.eq(a),'Mean'].to_numpy();y=z.loc[z.Type.eq(b),'Mean'].to_numpy();test=mannwhitneyu(x,y,alternative='two-sided',method='auto')
  local.append(dict(scenario=scenario,A=a,B=b,n_A=len(x),n_B=len(y),U=float(test.statistic),p_raw=float(test.pvalue),rank_biserial=float(2*test.statistic/(len(x)*len(y))-1)))
 idx=np.argsort([r['p_raw'] for r in local]);running=0
 for j,i in enumerate(idx):running=max(running,(len(local)-j)*local[i]['p_raw']);local[i]['p_Holm']=min(1,running)
 pairs.extend(local)
summary=pd.DataFrame(summ);overall=pd.DataFrame(overall);pairs=pd.DataFrame(pairs)
baseline=summary[summary.scenario.eq('Main_65')].set_index('Type')
summary['mean_change_percent']=[100*(r['mean']/baseline.loc[r['Type'],'mean']-1) for r in summ]
summary['median_change_percent']=[100*(r['median']/baseline.loc[r['Type'],'median']-1) for r in summ]
overall['mean_change_percent']=100*(overall['mean']/overall['mean'].iloc[0]-1);overall['median_change_percent']=100*(overall['median']/overall['median'].iloc[0]-1)
coverage['upper_excluded']=coverage.Site.isin(upper_excluded);coverage['Type']=coverage.Site.map(d.set_index('Site').Type)
coverage.to_csv(O/'site_eligibility.csv',index=False,encoding='utf-8-sig');summary.to_csv(O/'class_summary.csv',index=False,encoding='utf-8-sig');overall.to_csv(O/'network_summary.csv',index=False,encoding='utf-8-sig');pairs.to_csv(O/'pairwise_comparisons.csv',index=False,encoding='utf-8-sig')
labels={'Main_65':'Full network: 65 stations','Upper_complete_53':'Upper-bound coverage: 53 stations','Strict_49':'Full size coverage: 49 stations'}
report='# Particle-size-coverage sensitivity\n\n'
for scenario in sets:
 g=summary[summary.scenario.eq(scenario)]
 rank=' > '.join(g.sort_values('mean',ascending=False).Type)
 mrank=' > '.join(g.sort_values('median',ascending=False).Type)
 r=overall[overall.scenario.eq(scenario)].iloc[0]
 report+=f"{labels[scenario]}: mean ranking {rank}; median ranking {mrank}. Network mean {r['mean']:.2f} μm² cm⁻³; change {r.mean_change_percent:+.2f}%; Kruskal–Wallis H={r.H:.3f}, P={r.p:.6g}.\n\n"
report+='\n| Dataset | n | Mean | Median | Mean change (%) | H | P |\n|---|---:|---:|---:|---:|---:|---:|\n'
for r in overall.to_dict('records'):report+=f"| {labels[r['scenario']]} | {r['n']} | {r['mean']:.3f} | {r['median']:.3f} | {r['mean_change_percent']:+.2f}% | {r['H']:.3f} | {r['p']:.6g} |\n"
report+='\n| Station class | n: 65 → 49 | Mean: 65 → 49 | Median: 65 → 49 | Mean change (%) |\n|---|---:|---:|---:|---:|\n'
for t in order:
 a=summary[(summary.Type==t)&(summary.scenario=='Main_65')].iloc[0];z=summary[(summary.Type==t)&(summary.scenario=='Strict_49')].iloc[0]
 report+=f"| {t} | {a['n']}→{z['n']} | {a['mean']:.3f}→{z['mean']:.3f} | {a['median']:.3f}→{z['median']:.3f} | {z.mean_change_percent:+.2f}% |\n"
report+='\n| Comparison | Full network | Upper-bound coverage | Full size coverage |\n|---|---:|---:|---:|\n'
for a,b in itertools.combinations(order,2):
 v=[pairs[(pairs.scenario==s)&(pairs.A==a)&(pairs.B==b)].iloc[0].p_Holm for s in sets];report+=f'| {a} vs {b} | {v[0]:.6g} | {v[1]:.6g} | {v[2]:.6g} |\n'
report+='\nPairwise P values use the Holm correction.\n'
(O/'sensitivity_report.md').write_text(report)
(O/'Figure_size_coverage_sensitivity_caption.txt').write_text('Sensitivity of site-level LDSA comparisons to particle-size coverage.\n\na, Arithmetic means (symbols) and interquartile ranges across site-mean total LDSA values (horizontal lines) for all 65 sites (circles), 53 sites after excluding any contributing native configuration with an upper diameter below 500 nm (squares), and 49 sites meeting the strict-coverage criterion (triangles). Interquartile ranges describe between-site variation and are not confidence intervals. b, Percentage changes in the arithmetic mean of site-mean LDSA relative to the 65-site analysis, overall and by station class. Sites are weighted equally. The dashed line denotes no change. The samples overlap and are not treated as independent datasets.\n\nStatistics prioritize 2015–2019, using the nearest available year where that period has no data. Official invalid PNSD quality flags have been excluded while valid warning flags are retained. The strict subset requires all contributing native configurations to span 17–500 nm and excludes LND_TR, whose established LDSA integration starts at 31.65 nm. Station-class Kruskal–Wallis P values are '+', '.join(f'{r.p:.6g}' for r in overall.itertuples())+' for the 65-, 53- and 49-site sets, respectively.\n')
(O/'coverage_input_provenance.json').write_text(json.dumps(dict(source=str(src),sha256=hashlib.sha256(src.read_bytes()).hexdigest(),groups={k:len(v) for k,v in sets.items()}),ensure_ascii=False,indent=2))
print(overall.to_string(index=False));print(summary.to_string(index=False));print('SENSITIVITY_READY')
