import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--celltype", required=True)
parser.add_argument("--dataset", required=True)
parser.add_argument("--modules_csv_path", required=True)
parser.add_argument("--de_folder", required=True)
parser.add_argument("--fisher_csv_path", required=True)
parser.add_argument("--deg_cutoff", type=float, required=True)
parser.add_argument("--background_genes", required=True)
parser.add_argument("--celltype_folder_name", required=True)
args = parser.parse_args()

def fishers_enrichment_per_module(modules_df, background_genes_set, deg_df, celltype, deg_cutoff):
    results=[]
    deg_df=deg_df[deg_df['FDR']<deg_cutoff]
    if len(deg_df)==0:
        print(f"{args.celltype}: No differentially expressed genes at this cutoff value")
        return pd.DataFrame(columns=['celltype', 'module', 'overlap', 'odds_ratio', 'p_value', 'p_adj']) #returning empty pandas df if no test can be run
    deg_set=set(deg_df['gene'])
    for module in modules_df['module'].unique():
        module_genes=set(modules_df[modules_df['module']==module]['gene'])
        total_background=len(background_genes_set)
        a=len(deg_set&module_genes)
        b=len(deg_set-module_genes)
        c=len(module_genes-deg_set)
        d=total_background-(a+b+c)
        
        table=[[a, b],[c, d]]
        oddsratio, p_value=fisher_exact(table, alternative='greater') 
        results.append({'celltype':celltype,'module': module,'overlap': a,'odds_ratio': oddsratio,'p_value': p_value})
    results_df=pd.DataFrame(results)
    try:
        if results_df.empty or results_df['p_value'].isnull().all():
            raise ValueError("No valid p-values to adjust.")
        results_df['p_adj']=multipletests(results_df['p_value'], method='fdr_bh')[1]
    except Exception as e:
        print(f"{celltype}: Skipping p-value adjustment — {e}")
    return(results_df)

de_csv_path=f"{args.de_folder}/{args.dataset}_DEG/results/{args.dataset}_{args.celltype_folder_name}__{args.celltype}_edger_results.csv"
deg_df=pd.read_csv(de_csv_path)
background_genes_set=set(pd.read_csv(args.background_genes)['gene'])
modules_df=pd.read_csv(args.modules_csv_path)

mod_results=fishers_enrichment_per_module(modules_df=modules_df,background_genes_set=background_genes_set, deg_df=deg_df, celltype=args.celltype, deg_cutoff=args.deg_cutoff)
mod_results.to_csv(args.fisher_csv_path, index=False)