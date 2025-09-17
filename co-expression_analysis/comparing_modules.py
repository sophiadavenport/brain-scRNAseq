import pandas as pd
import numpy as np
import os
import itertools
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--liana_csv", required=True)
parser.add_argument("--dataset", required=True)
parser.add_argument("--celltype_list", nargs="+", required=True)
parser.add_argument("--correlation_df_path", required=True)
parser.add_argument("--all_celltype_modules_path", required=True)
args = parser.parse_args()

dataset=args.dataset
celltypes=args.celltype_list

all_celltype_modules=[]
indv_exp=[]
for celltype in celltypes:
    indv_mod_exp=pd.read_csv('results/'+dataset+'__'+celltype+'/'+dataset+'__'+celltype+"_individualmatrix.csv")
    indv_mod_exp['module']=celltype+'_'+indv_mod_exp['Unnamed: 0'].astype(str)
    indv_mod_exp=indv_mod_exp.set_index('module').drop('Unnamed: 0', axis=1)    
    indv_exp.append(indv_mod_exp)
    try:
        ks_stats=pd.read_csv('results/'+dataset+'__'+celltype+'/module_geneexp_dist/'+dataset+'_KSstats.csv')
        ks_stats['module']=celltype+'_'+ks_stats['module'].astype(str)
        ks_stats=ks_stats.rename(columns={'confidence': 'ks_confidence'})
        ks_stats=ks_stats[['module', 'median_ks_pval', 'ks_confidence']]
    except:
        print('no KS stats:', dataset, celltype)
    fisher=pd.read_csv('results/'+dataset+'__'+celltype+'/DEG_mod_association/'+dataset+'_fisherstats.csv')
    fisher['module']=celltype+'_'+fisher['module'].astype(str)
    fisher=fisher.rename(columns={'p_adj': 'fisher_p_adj', 'odds_ratio': 'fisher_odds_ratio'})
    fisher=fisher[['module', 'fisher_odds_ratio', 'fisher_p_adj']]
    all_celltype_modules.append(pd.merge(fisher, ks_stats, on='module').sort_values(by=['fisher_p_adj', 'median_ks_pval'], ascending=True))
    
dataset_indv_exp=pd.concat(indv_exp).T #need to transpose else will have individuals as correlations
dataset_modules=pd.concat(all_celltype_modules).sort_values(by=['median_ks_pval', 'fisher_p_adj'], ascending=True)

dataset_modules.to_csv(args.all_celltype_modules_path)
#spearman correlation between modules from different cell types
results = []
for m1, m2 in itertools.combinations(dataset_indv_exp.columns, 2):
    if m1.rsplit('_', 1)[0] == m2.rsplit('_', 1)[0]:
        continue #skip same celltype comparisons
    x=dataset_indv_exp[m1]
    y=dataset_indv_exp[m2]
    mask = x.notna() & y.notna()
    x=x[mask];y=y[mask]
    if x.nunique() < 2 or y.nunique() < 2: #skip constant vectors (spearmanr returns nan if constant)
        continue
    corr, pval = spearmanr(x, y, nan_policy='omit')
    results.append({"module1": m1, "module2": m2,  "correlation": float(corr), "pvalue": float(pval), "n": int(len(x))})   #n is number of samples used

cor_df=pd.DataFrame(results)
cor_df['padj']=multipletests(cor_df['pvalue'], method='fdr_bh')[1] #FDR correction across all tested pairs
cor_df=cor_df.sort_values(by='padj', ascending=True)

cor_df.to_csv(args.correlation_df_path)