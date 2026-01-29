#run in scanpy env
import scanpy as sc
import argparse
import scanpy_helper_functions as sh # type: ignore

parser = argparse.ArgumentParser()
parser.add_argument("--joined_adata", required=True, help="Joined adata from join_adataspy (does not have umap rerun yet)")
parser.add_argument("--joined_outfile", required=True)
args = parser.parse_args()

adata_combined=sc.read_h5ad(args.joined_adata)

print('input adata:', adata_combined)
cur_max=adata_combined.X.max()
print('current max count:', cur_max)
print('adata unique_subtype:', adata_combined.obs.unified_subtype.unique())
print('adata Celltype_Class:', adata_combined.obs.Celltype_Class.unique())

if cur_max > 100:
    print('adata input with raw counts, setting to logged before continuing...')
    sc.pp.normalize_total(adata_combined)
    sc.pp.log1p(adata_combined)

sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000) #rerun highly variable genes with combined data
sc.pp.pca(adata_combined) #rerun PCA with combined data before harmony integration
sc.pl.pca_variance_ratio(adata_combined, n_pcs=50, log=True) #pca plot
sc.external.pp.harmony_integrate(adata_combined, key='datasource')
sc.pp.neighbors(adata_combined, use_rep='X_pca_harmony')
sc.tl.umap(adata_combined)

adata_combined=sh.set_celltype_colors(adata=adata_combined, celltype_class_col='Celltype_Class', celltype_subtype_col='unified_subtype')

adata_combined.write(args.joined_outfile)