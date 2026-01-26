import sys
import scanpy as sc
import scanpy_helper_functions as sh
import numpy as np
import anndata as an
import sklearn
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--BD_adata", required=True)
parser.add_argument("--SZ_counts", required=True)
parser.add_argument("--SZ_metadata", required=True)
parser.add_argument("--SZ_genenames", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

BD=sc.read_h5ad(args.BD_adata)

SZ=sh.mtx_to_adata(counts_path=args.SZ_counts, 
                   gene_names_path=args.SZ_genenames, 
                   metadata_path=args.SZ_metadata,
                   genename_col='Gene', barcodes_col='Barcodes', 
                   data_name=None, directory_path=None, raw_counts_path=None, spatial_path=None)

for adata in [BD, SZ]:
    if adata.X.max()<20:
        adata.layers['counts']=adata.X #saving copy of raw counts to adata object under counts
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    if 'Celltype' not in adata.obs.columns: #standardizing celltype column name
        if 'celltype' in adata.obs.columns:
            adata.obs.rename(columns={'celltype': 'Celltype'}, inplace=True)
        elif 'cell_group' in adata.obs.columns:
            adata.obs.rename(columns={'cell_group': 'Celltype'}, inplace=True)
        else:
            print('error celltype column not found to be standardized')
    if 'Subtype' not in adata.obs.columns: #standardizing subtype column name
        if 'subtype' in adata.obs.columns:
            adata.obs.rename(columns={'subtype': 'Subtype'}, inplace=True)
        elif 'Celltype_inferred' in adata.obs.columns:
            adata.obs.rename(columns={'Celltype_inferred': 'Subtype'}, inplace=True)
        else:
            adata.obs['Subtype']=np.full(adata.obs.shape[0], "NA")
    if 'Condition' not in adata.obs.columns: #standardizing condition column name
        if 'Phenotype' in adata.obs.columns:
            adata.obs.rename(columns={'Phenotype': 'Condition'}, inplace=True)
        elif 'phenotype' in adata.obs.columns:
            adata.obs.rename(columns={'phenotype': 'Condition'}, inplace=True)
        elif 'Diagnosis' in adata.obs.columns:
            adata.obs.rename(columns={'Diagnosis': 'Condition'}, inplace=True)
        else:
            print('error condition/phenotype column not found to be standardized')
    if 'Age' not in adata.obs.columns:
        if 'age' in adata.obs.columns:
            adata.obs.rename(columns={'age': 'Age'}, inplace=True)
        else:
            adata.obs['Age'] = pd.Series([pd.NA] * adata.obs.shape[0], dtype='Int64')
    else:
         adata.obs.Age=adata.obs.Age.round().astype('Int64') #Original BD did not include age column only SZ
    if 'Individual' not in adata.obs.columns:
        if 'id' in adata.obs.columns:
            adata.obs.rename(columns={'id': 'Individual'}, inplace=True)
        elif 'SubID' in adata.obs.columns:
            adata.obs.rename(columns={'SubID': 'Individual'}, inplace=True)
        else:
            print('error individual column not found to be standardized')
    if 'Gender' in adata.obs.columns:
        adata.obs.rename(columns={'Gender': 'Sex'}, inplace=True)
    elif 'sex' in adata.obs.columns:
        adata.obs.rename(columns={'sex': 'Sex'}, inplace=True)
    else:
        adata.obs['Sex'] = np.full(adata.obs.shape[0], 'NA') #BD does not include sex column only SZ
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata)

BD.obs['data_source'] = np.full(BD.obs.shape[0], "Han")
SZ.obs['data_source'] = np.full(SZ.obs.shape[0], "Ruzicka")

try:
    BD_barcodes=BD.obs_names
except:
    try:
        BD_barcodes=BD.obs['Unnamed: 0'].unique()
    except:
        print('missing BD barcodes')
SZ_barcodes=SZ.obs['Barcodes'].unique()
overlap_cells = set(SZ_barcodes) & set(BD_barcodes)
print('Number of Overlapping Cells: ', len(overlap_cells))

BD = BD[~BD.obs_names.isin(overlap_cells)]
print('BD shape post removal of cells also in SZ: ', BD.shape)

BD_SZ = an.concat([BD, SZ], join='inner', keys=["Han", "Ruzicka"]) #only choosing overlapping genes (only same column names in obs will be kept!)
print('Combined adata shape:', BD_SZ.shape)

adatas_list=sh.check_adata_format(adatas=[BD_SZ], batches=['Batch'],
                                      data_sources=['Han_Ruzicka'],
                                      celltype_colnames=['Celltype'], 
                                      subtype_colnames=["Subtype"], 
                                      condition_colnames=['Condition'],
                                      age_colnames=['Age'],
                                      sex_colnames=['Sex'])

BD_SZ_form=adatas_list[0] #since check adata format returns a list

BD_SZ_form.obs.drop(columns='datasource', inplace=True) #dropping the datasource column with Han_Ruzicka as entry (will be needed separated in future analysis)
BD_SZ_form.obs.rename(columns={'data_source': 'datasource'}, inplace=True) #replacing datasource column with 'Han', 'Ruzicka' as values

BD_SZ_form.write(args.output)