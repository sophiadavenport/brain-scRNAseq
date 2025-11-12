import scanpy as sc
import numpy as np
import pandas as pd
from scipy import io as scipy_io
from scipy import sparse as sp_sparse
from scipy.sparse import csr_matrix
import anndata as ad
import argparse
import sys
import os

parser=argparse.ArgumentParser()
parser.add_argument("--adata", required=True)
parser.add_argument("--celltype_col", required=True)
parser.add_argument("--output_prefix", required=True)
args=parser.parse_args()

os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

def write_empty_outputs(prefix, current_celltype):
    os.makedirs(prefix, exist_ok=True)
    empty_matrix = csr_matrix((1,1), dtype=np.float32)
    scipy_io.mmwrite(f"{prefix}/{current_celltype}_counts.mtx", empty_matrix)
    np.savetxt(f"{prefix}/{current_celltype}_genes.tsv", np.array([]), fmt="%s")
    np.savetxt(f"{prefix}/{current_celltype}_barcodes.tsv", np.array([]), fmt="%s")
    # optional empty metadata file
    open(f"{prefix}/{current_celltype}_metadata.csv", "w").close()

adata=sc.read_h5ad(args.adata, backed='r')

adata.obs["Celltype_clean"]=adata.obs[args.celltype_col].astype(str).str.replace(r" \+ ", "_", regex=True).str.replace(r" ", "_", regex=True).str.replace(r"\(", "", regex=True).str.replace(r"\)", "", regex=True).str.replace("/", "-")

celltypes=adata.obs["Celltype_clean"].unique()

for current_celltype in celltypes:
    cell_adata=adata[adata.obs["Celltype_clean"]==current_celltype, :].to_memory()
    print('current_celltype:', cell_adata.shape)

    if cell_adata.n_obs < 100: #setting minimum number of cells to be 100
        print(f"Skipping {current_celltype}: only {cell_adata.n_obs} cells")
        write_empty_outputs(args.output_prefix, current_celltype)
    else:
        cell_adata_counts=cell_adata.layers['counts']
        
        scipy_io.mmwrite(f"{args.output_prefix}/{current_celltype}_counts.mtx", cell_adata_counts)

        cell_adata.var_names.to_series().to_csv(f"{args.output_prefix}/{current_celltype}_genes.tsv", sep="\t", header=False)
        cell_adata.obs_names.to_series().to_csv(f"{args.output_prefix}/{current_celltype}_barcodes.tsv", sep="\t", header=False)

        #Metadata fields if available
        meta_cols=[x for x in ["Age","Sex", "pct_mito", "pct_counts_mt", "Celltype_Class", "unified_subtype"] if x in cell_adata.obs.columns]
        cell_adata.obs[meta_cols].to_csv(f"{args.output_prefix}/{current_celltype}_metadata.csv")

