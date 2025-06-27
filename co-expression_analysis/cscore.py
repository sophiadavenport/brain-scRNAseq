import scanpy as sc
import numpy as np
import pandas as pd
from CSCORE import CSCORE
import anndata as ad
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--adata", required=True)
parser.add_argument("--celltype_col", required=True)
parser.add_argument("--current_celltype", required=True)
parser.add_argument("--output_prefix", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.adata)

cell_adata = adata[adata.obs[args.celltype_col]==args.current_celltype, :]

cell_adata.write("/lab/solexa_sun/shared_Data/scRNAseq/brain/combined_data/AD_microglia_6_25.h5ad")

if cell_adata.X.max() < 100: #already normalized... resetting to raw counts
    cell_adata.X=cell_adata.layers["counts"]
if cell_adata.raw is None:
    cell_adata.raw = ad.AnnData(X=cell_adata.layers['counts'].copy(),var=cell_adata.var.copy(),obs=cell_adata.obs.copy())

#scale counts by seq depths
sc.pp.normalize_total(cell_adata, target_sum=1)
mean_exp = (cell_adata.X).sum(axis=0).A1
mean_exp_df = pd.DataFrame({'gene': cell_adata.var_names, 'mean_expression': mean_exp})
top_genes_df = mean_exp_df.sort_values(by='mean_expression', ascending=False).head(5000) #selecting top 5000 genes by expression
top_genes_indices = top_genes_df.index.astype(int).to_numpy()

if "CON" in cell_adata.obs.Condition.unique():
    cell_adata_control = cell_adata[cell_adata.obs["Condition"] == 'CON', :]
elif "Control" in cell_adata.obs.Condition.unique():
    cell_adata_control = cell_adata[cell_adata.obs["Condition"] == 'Control', :]
else:
    print("Unsure control coding in dataset... cannot proceed")
    sys.exit(1)

cell_adata_disease = cell_adata[cell_adata.obs["Condition"] != 'Control', :]

res_control = CSCORE(cell_adata_control, top_genes_indices)
res_disease = CSCORE(cell_adata_disease, top_genes_indices)

#saving outputs for R...
np.savetxt(f"{args.output_prefix}_ctrl_coexpr.csv", res_control[0], delimiter=",")
np.savetxt(f"{args.output_prefix}_ctrl_pvals.csv", res_control[1], delimiter=",")
np.savetxt(f"{args.output_prefix}_ctrl_teststats.csv", res_control[2], delimiter=",")

np.savetxt(f"{args.output_prefix}_disease_coexpr.csv", res_disease[0], delimiter=",")
np.savetxt(f"{args.output_prefix}_disease_pvals.csv", res_disease[1], delimiter=",")
np.savetxt(f"{args.output_prefix}_disease_teststats.csv", res_disease[2], delimiter=",")

#Save selected gene names (same order as matrices)
genes_selected = top_genes_df['gene'].to_numpy()
np.savetxt(f"{args.output_prefix}_ctrl_genes.csv", genes_selected, fmt='%s')
np.savetxt(f"{args.output_prefix}_disease_genes.csv", genes_selected, fmt='%s')