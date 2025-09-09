#script will run liana but will not filter any of the results, filtering to be done downstream
import scanpy as sc
import pandas as pd
from liana.mt import rank_aggregate
import argparse
import os
from liana.method import cellphonedb, connectome, cellchat

parser=argparse.ArgumentParser()
parser.add_argument('--adata', required=True, help='Path to h5ad single-cell dataset')
parser.add_argument('--celltype_col', default='Celltype_Class', help='Column in adata.obs for cell type labels')
parser.add_argument('--output_csv', required=True, help='Path to output CSV of LIANA results')
args=parser.parse_args()

adata=sc.read_h5ad(args.adata)

#cellphonedb(adata, groupby='Celltype_Class', resource_name='consensus', expr_prop=0.1, verbose=True, key_added='cpdb_res', use_raw=False, n_perms=10)
#connectome(adata, groupby='Celltype_Class', resource_name='consensus',expr_prop=0.1, verbose=True, key_added='conn_res', use_raw=False)

rank_aggregate(adata,groupby=args.celltype_col,use_raw=False, resource_name='consensus', key_added='liana_res')

liana_results=adata.uns['liana_res']

os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)
liana_results.to_csv(args.output_csv, index=False)
print(f"LIANA results saved to {args.output_csv}")
