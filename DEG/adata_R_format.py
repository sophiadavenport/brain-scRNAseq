import scanpy as sc
import pandas as pd
import argparse
from scipy import io, sparse
from scipy.io import mmwrite

parser = argparse.ArgumentParser()
parser.add_argument("--adata", required=True)
parser.add_argument("--dataset", required=True)
parser.add_argument("--celltype", required=True)
parser.add_argument("--split_col", required=True)
parser.add_argument("--counts_path", required=True)
parser.add_argument("--metadata_path", required=True)
parser.add_argument("--barcodes_path", required=True)
parser.add_argument("--features_path", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.adata)

#Clean names
adata.obs[args.split_col] = adata.obs[args.split_col].astype(str).str.replace(" ", "_").str.replace("+", "").str.replace("/", "_").str.replace("__","_")

#Filter to the specified cell type
subset = adata[adata.obs[args.split_col] == args.celltype].copy()

#Use raw counts
if "counts" in subset.layers:
    subset.X = subset.layers["counts"]
else:
    print(args.dataset, args.celltype, 'Counts not available in layers. Current maximum counts is:', subset.X.max(), '\n')

#Metadata safety checks
if 'pct_counts_mt' not in subset.obs.columns:
    if 'percent.mt' in subset.obs.columns:
        subset.obs.rename(columns={'percent.mt': 'pct_counts_mt'}, inplace=True)
    elif 'percent_mt' in subset.obs.columns:
        subset.obs.rename(columns={'percent_mt': 'pct_counts_mt'}, inplace=True)
    else:
        print(f'[WARN] Could not find mitochondrial column in {args.dataset} for {args.celltype}')

if 'n_genes_by_counts' not in subset.obs.columns:
    if 'nFeature_RNA' in subset.obs.columns:
        subset.obs.rename(columns={'nFeature_RNA': 'n_genes_by_counts'}, inplace=True)
    else:
        print(f'[WARN] Could not find feature count column in {args.dataset} for {args.celltype}')

#Save count matrix (.mtx format)
matrix = subset.X.T #transposing (edgeR wants genes x cells)
mmwrite(args.counts_path, matrix)

#Saving metadata, features, and barcode
subset.obs_names.to_series().to_csv(args.barcodes_path, index=False, header=False)
subset.var_names.to_series().to_csv(args.features_path, index=False, header=False)
subset.obs.to_csv(args.metadata_path)