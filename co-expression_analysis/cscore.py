import scanpy as sc
import numpy as np
import pandas as pd
from CSCORE import CSCORE
import anndata as ad
import argparse
import sys
import os
import gc

parser=argparse.ArgumentParser()
parser.add_argument("--adata", required=True)
parser.add_argument("--celltype_col", required=True)
parser.add_argument("--current_celltype", required=True)
parser.add_argument("--output_prefix", required=True)
args=parser.parse_args()

os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

adata=sc.read_h5ad(args.adata)

adata.obs["Celltype_clean"]=adata.obs[args.celltype_col].astype(str).str.replace(r" \+ ", "_", regex=True).str.replace(r" ", "_", regex=True).str.replace(r"\(", "", regex=True).str.replace(r"\)", "", regex=True).str.replace("/", "-")

cell_adata=adata[adata.obs["Celltype_clean"]==args.current_celltype, :]

print(cell_adata.shape)
cell_adata=cell_adata.copy() #clean copy

del adata #clearing up memory
gc.collect()

if cell_adata.X.max() < 100: #already normalized... resetting to raw counts
    cell_adata.X=cell_adata.layers["counts"]
if cell_adata.raw is None:
    cell_adata.raw=ad.AnnData(X=cell_adata.layers['counts'].copy(),var=cell_adata.var.copy(),obs=cell_adata.obs.copy())

#scale counts by seq depths
sc.pp.normalize_total(cell_adata, target_sum=1)
nonzero_expr=(cell_adata.X > 0)
frac_expressed=nonzero_expr.sum(axis=0).A1 / cell_adata.X.shape[0]
expr_frac_df=pd.DataFrame({'gene': cell_adata.var_names,'frac_expressed': frac_expressed})
expressed_genes_df=expr_frac_df[expr_frac_df['frac_expressed'] >= 0.10]
expressed_genes_indices=expressed_genes_df.index.astype(int).to_numpy()
print("number of genes in over 10pct of cells:", len(expressed_genes_df['gene']))
res=CSCORE(cell_adata, expressed_genes_indices)

#saving outputs for R...
np.savetxt(f"{args.output_prefix}_coexpr.csv", res[0], delimiter=",")
np.savetxt(f"{args.output_prefix}_pvals.csv", res[1], delimiter=",")
np.savetxt(f"{args.output_prefix}_teststats.csv", res[2], delimiter=",")

#Save selected gene names (same order as matrices)
genes_selected=expressed_genes_df['gene'].to_numpy()
np.savetxt(f"{args.output_prefix}_genes.csv", genes_selected, fmt='%s')