"""Spatial comparisons for one 65-station network and coverage subsets."""
import itertools
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

def compare_coverage(d):
    """Return class, overall and Holm-adjusted pairwise summaries."""
    sets={'Main_65':d,'Upper_complete_53':d[d.upper_coverage],'Strict_49':d[d.full_coverage]}
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
    return summary, overall, pairs
