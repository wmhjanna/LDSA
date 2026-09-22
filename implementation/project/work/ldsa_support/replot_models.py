from pathlib import Path
import sys,textwrap
ROOT=Path('@PROJECT_ROOT@');sys.path.insert(0,str(ROOT/'work'));sys.path.insert(0,str(ROOT/'work/ldsa_support/deps'))
import model_cleaned_ldsa as m
import pandas as pd
models={(typ,phase):pd.read_csv(m.ANA/f'heldout_SHAP_{typ}_{phase}.csv') for typ in ['UB','TR'] for phase in ['Before','After']}
imp=pd.read_csv(m.ANA/'SHAP_importance.csv')
text=(ROOT/'work/model_cleaned_ldsa.py').read_text();body=text.split('    fig,axs=plt.subplots(2,2,figsize=(12,10),layout=\'constrained\')',1)[1].split('\nif __name__',1)[0]
code="fig,axs=plt.subplots(2,2,figsize=(12,10),layout='constrained')\n"+textwrap.dedent(body)
scope=vars(m).copy();scope.update(models=models,imp=imp);exec(compile(code,__file__,'exec'),scope)
