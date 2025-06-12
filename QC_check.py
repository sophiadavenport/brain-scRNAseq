import scanpy_helper_functions as sh
import scanpy as sc
import argparse
import gc
import pandas as pd
from pybiomart import Server, Dataset
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--MS_counts", required=True)
parser.add_argument("--MS_metadata", required=True)
parser.add_argument("--MS_genenames", required=True)
parser.add_argument("--MS_outfile", required=True)
parser.add_argument("--AD_h5ad", required=True)
parser.add_argument("--AD_metadata", required=True)
parser.add_argument("--AD_outfile", required=True)
parser.add_argument("--BD_SZ", required=True)
parser.add_argument("--BD_SZ_outfile", required=True)
parser.add_argument("--keep_genes_file", required=True)
args = parser.parse_args()

##Adding metadata to AD:
AD=sc.read_h5ad(args.AD_h5ad)
AD_metadata=pd.read_csv(args.AD_metadata)
print('Starting AD: ', AD.shape)
print('\nadding additional metadata info for AD...\n')
AD.obs=pd.merge(AD.obs, AD_metadata, how='left')
AD.obs['Condition'] = np.where(AD.obs['Pathologic_diagnosis_of_AD'] == 'yes', 'AD','Control') #adding condition column

try:
    len(AD.obs.cell_type_high_resolution)
    print('celltype in adata.obs for AD\n')
except:
    print('AD is still missing celltypes\n')

##Getting MS into adata format 
MS=sh.mtx_to_adata(counts_path=args.MS_counts, 
                   gene_names_path=args.MS_genenames, 
                   metadata_path=args.MS_metadata,
                   genename_col='symbol', barcodes_col='cell_id', 
                   data_name=None, directory_path=None, raw_counts_path=None, spatial_path=None)

if MS.var_names[0].upper().startswith('ENSG'): #if ensemblIDs are set as var names change to gene names
    print('MS has ensembl ids as var_names...switching...\n')
    MS=sh.add_gene_names_to_adata(MS)
    print('MS var names reset\n')

AD_MS=sh.check_adata_format(adatas=[MS, AD], batches=['sample_source', 'batch'],
                                      data_sources=['Macnair', 'Mathys'],
                                      celltype_colnames=['type_broad', 'major.celltype'], 
                                      subtype_colnames=['type_fine', 'cell_type_high_resolution'], 
                                      condition_colnames=['diagnosis', 'Condition'],
                                      age_colnames=['age_at_death', 'age_death_x'],
                                      sex_colnames=['sex', 'msex_x'])

MS_form=AD_MS[0]
#MS_form.write(args.MS_output)
AD_form=AD_MS[1]

#reducing memory by removing unneccesary objects downstream
del MS
del AD
del AD_MS
gc.collect()

BD_SZ=sc.read_h5ad(args.BD_SZ)

#removing genes not included in all adatas
common_genes = set(BD_SZ.var_names) & set(AD_form.var_names) & set(MS_form.var_names)
common_genes_df = pd.DataFrame(sorted(common_genes), columns=["gene"])
common_genes_df.to_csv(args.keep_genes_file, index=False)
print('Common Genes:',len(common_genes),'\n')
AD_form = AD_form[:, AD_form.var_names.isin(common_genes)].copy()
AD_form.write(args.AD_outfile)
print('Formatted AD Shape:', AD_form.shape)
del AD_form
gc.collect()
MS_form = MS_form[:, MS_form.var_names.isin(common_genes)].copy()
MS_form.write(args.MS_outfile)
print('Formatted MS Shape:', MS_form.shape)
del MS_form
gc.collect()
BD_SZ = BD_SZ[:, BD_SZ.var_names.isin(common_genes)].copy()
BD_SZ.write(args.BD_SZ_outfile)
print('Formatted BD_SZ Shape:', BD_SZ.shape)