import scanpy as sc
import scanpy_helper_functions as sh # type: ignore
import numpy as np
import anndata as an
import argparse
import gc

parser = argparse.ArgumentParser()
parser.add_argument("--AD", required=True)
parser.add_argument("--MS", required=True)
parser.add_argument("--BD_SZ", required=True)
parser.add_argument("--joined_outfile", required=True)
parser.add_argument("--umap_path", required=True)
args = parser.parse_args()

AD=sc.read_h5ad(args.AD)
MS=sc.read_h5ad(args.MS)
adatas_list=[AD, MS]
data_sources=['Mathys', 'Macnair']
for i, adata in enumerate(adatas_list):
    data_source=data_sources[i]
    assert 'Celltype' in adata.obs, f"Missing 'Celltype' in {adata}"
    assert 'Subtype' in adata.obs, f"Missing 'Subtype' in {adata}"
    assert 'datasource' in adata.obs, f"Missing 'datasource' in {adata}"
    assert 'Condition' in adata.obs, f"Missing 'Condition' in {adata}"
    assert 'unified_celltype' in adata.obs, f"Missing 'Condition' in {adata}"
    assert 'unified_celltype_broad' in adata.obs, f"Missing 'Condition' in {adata}"
    assert 'Sex' in adata.obs, f"Missing 'Condition' in {adata}"
    #ensure no duplicated genes:
    duplicated_genes = adata.var.index[adata.var.index.duplicated()]
    if len(duplicated_genes) > 0:
        print('adjusting an adata for duplicated genes...')
        adata.var_names_make_unique()
     #ensure no overlapping indexes
    adata.obs['Barcode'] = adata.obs_names
    adata.obs_names = [f"{data_source}-" + barcode for barcode in adata.obs_names]
print('First joining adatas...\n')
adata_combined = an.concat(adatas_list, join='inner', keys=data_sources) #only same column names in obs will be kept!
print('First combined adata shape: ',adata_combined.shape, '\n')

del MS
del AD
gc.collect()

BD_SZ=sc.read_h5ad(args.BD_SZ)

if 'data_sources' not in BD_SZ.obs:
    BD_SZ.obs['data_sources']=BD_SZ.obs['data_source']
if 'datasource' not in BD_SZ.obs:
    BD_SZ.obs['datasource']='Han_Ruzicka'

assert 'Celltype' in BD_SZ.obs, f"Missing 'Celltype' in {BD_SZ}"
assert 'Subtype' in BD_SZ.obs, f"Missing 'Subtype' in {BD_SZ}"
assert 'datasource' in BD_SZ.obs, f"Missing 'datasource' in {BD_SZ}"
assert 'Condition' in BD_SZ.obs, f"Missing 'Condition' in {BD_SZ}"
assert 'unified_celltype' in BD_SZ.obs, f"Missing 'Condition' in {BD_SZ}"
assert 'unified_celltype_broad' in BD_SZ.obs, f"Missing 'Condition' in {BD_SZ}"
assert 'Sex' in BD_SZ.obs, f"Missing 'Condition' in {BD_SZ}"

adatas_list=[adata_combined, BD_SZ]
print('Trying join 2...\n')
adata_combined = an.concat(adatas_list, join='inner', keys=None) #trying with no key for second join ()
print('Join 2 successful:', adata_combined.shape, '\n')
print('Adata combined:', adata_combined, '\n')

del BD_SZ
gc.collect()

sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000) #rerun highly variable genes with combined data
sc.pp.pca(adata_combined) #rerun PCA with combined data before harmony integration
sc.external.pp.harmony_integrate(adata_combined, key='datasource')
sc.pp.neighbors(adata_combined, use_rep='X_pca_harmony')
sc.tl.umap(adata_combined)

adata_combined.write(args.joined_outfile)

###UMAPs
sc.pl.umap(adata_combined, color='unified_celltype', size=1, save=args.umap_path, show=False)
sc.pl.umap(adata_combined, color='unified_celltype_broad', size=1, save="results/ADBDMSSZ_full_unified_celltype__broad_umap.png",show=False)