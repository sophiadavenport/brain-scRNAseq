import scanpy as sc
import scanpy_helper_functions as sh # type: ignore
import numpy as np
import anndata as an
import argparse
import gc

parser = argparse.ArgumentParser()
parser.add_argument("--AD", required=True)
parser.add_argument("--MS", required=True)
parser.add_argument("--BD", required=True)
parser.add_argument("--SZ", required=True)
parser.add_argument("--joined_outfile", required=True)
args = parser.parse_args()

def assert_correct_cols_checkobsnames(adata, name, author):
    assert 'Celltype' in adata.obs.columns, f"Missing 'Celltype' in {name}"
    assert 'datasource' in adata.obs.columns, f"Missing 'datasource' in {name}"
    assert 'Condition' in adata.obs.columns, f"Missing 'Condition' in {name}"
    assert 'unified_subtype' in adata.obs.columns, f"Missing 'unified_subtype' in {name}"
    assert 'Celltype_Class' in adata.obs.columns, f"Missing 'Celltype_Class' in {name}"
    assert 'Sex' in adata.obs.columns, f"Missing 'Sex' in {name}"
    #ensure no overlapping indexes
    adata.obs['Barcode'] = adata.obs_names
    adata.obs_names = [f"{author}-" + barcode for barcode in adata.obs_names]
    return(adata)

AD=sc.read_h5ad(args.AD)
MS=sc.read_h5ad(args.MS)
adatas_list=[AD, MS]
data_sources=['Mathys', 'Macnair']
assert_correct_cols_checkobsnames(adata=AD, name='AD', author='Mathys')
assert_correct_cols_checkobsnames(adata=MS, name='MS', author='Macnair')
print('First joining adatas...\n')
adata_combined = an.concat(adatas_list, join='inner', keys=data_sources) #only same column names in obs will be kept!
print('First combined adata shape: ',adata_combined.shape, '\n')

del MS
del AD
gc.collect()

BD=sc.read_h5ad(args.BD)

if 'datasource' not in BD.obs:
    BD.obs['datasource']='Han'

assert_correct_cols_checkobsnames(adata=BD, name='BD', author='Han')

adatas_list=[adata_combined, BD]
print('Trying join 2...\n')
adata_combined = an.concat(adatas_list, join='inner', keys=None) #trying with no key for second join ()
print('Join 2 successful:', adata_combined.shape, '\n')
print('Adata combined:', adata_combined, '\n')

del BD
gc.collect()

SZ=sc.read_h5ad(args.SZ)

if 'datasource' not in SZ.obs:
    SZ.obs['datasource']='Ruzicka'

assert_correct_cols_checkobsnames(adata=SZ, name='SZ', author='Ruzicka')

adatas_list=[adata_combined, SZ]
print('Trying join 3...\n')
adata_combined = an.concat(adatas_list, join='inner', keys=None) #trying with no key for second join ()
print('Join 3 successful:', adata_combined.shape, '\n')
print('Adata combined:', adata_combined, '\n')

del SZ
gc.collect()

adata_combined.write(args.joined_outfile)