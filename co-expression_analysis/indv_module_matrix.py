import scanpy as sc
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--adata", required=True)
parser.add_argument("--modules", required=True)
parser.add_argument("--matrix", required=True)
parser.add_argument("--celltype_col", required=True)
parser.add_argument("--current_celltype", required=True)
args = parser.parse_args()

def get_individuallevel_module_exp(adata, module_df):
    modules = module_df['module'].unique()
    if 'Individual' in list(adata.obs.columns):
        individual_col='Individual'
        individuals = adata.obs['Individual'].unique()
    elif 'individual_id_anon' in list(adata.obs.columns):
        individual_col='individual_id_anon'
        individuals = adata.obs['individual_id_anon'].unique()
    else:
        individual_col='individualID'
        individuals = adata.obs['individualID'].unique()
    matrix = pd.DataFrame(index=modules, columns=individuals, dtype=float)
    for module_num in list(module_df['module'].unique()):
        gene_list=list(module_df[module_df['module']==module_num]['gene'])   
        mod_adata=adata[:, adata.var_names.isin(gene_list)]

        individual_means = {}

        for ind in list(mod_adata.obs[individual_col].unique()):
            ind_adata=mod_adata[mod_adata.obs[individual_col]==ind]
            gene_means = np.mean(ind_adata.layers['normalized_counts'], axis=0)
            module_mean = np.mean(gene_means)
            individual_means[ind] = module_mean
        
        matrix.loc[module_num] = pd.Series(individual_means)
    
    return(matrix)

adata = sc.read_h5ad(args.adata)

adata.obs["Celltype_clean"] = adata.obs[args.celltype_col].astype(str).str.replace(r" \+ ", "_", regex=True).str.replace(r" ", "_", regex=True).str.replace(r"\(", "", regex=True).str.replace(r"\)", "", regex=True).str.replace("/", "-")

cell_adata=adata[adata.obs["Celltype_clean"]==args.current_celltype, :]

modules_df=pd.read_csv(args.modules)

matrix=get_individuallevel_module_exp(cell_adata, modules_df)

matrix.to_csv(args.matrix)