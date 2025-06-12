#Mapping Annotations to AD data
#must use TACCO_env
import pandas as pd
import scanpy as sc
import numpy as np
import anndata as an
import tacco as tc # type: ignore
import scanpy_helper_functions as sh
import sklearn
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--AD", required=True)
parser.add_argument("--BD_SZ", required=True)
parser.add_argument("--MS", required=True)
parser.add_argument("--mapping_annotation_key", required=True)
parser.add_argument("--AD_outfile", required=True)
parser.add_argument("--BD_SZ_outfile", required=True)
parser.add_argument("--MS_outfile", required=True)
parser.add_argument("--report", required=True, help="Relative path to report to be saved as txt.")
args = parser.parse_args()

adata_mapping=sc.read_h5ad(args.AD)
if adata_mapping.X.max() < 100: #then log normalized
    adata_mapping.layers['normalized_counts']=adata_mapping.X #saving normalized counts before resetting
    adata_mapping.X=adata_mapping.layers['counts'] #resetting to raw counts since TACCO requires it

if 'Celltype_Class' not in adata_mapping.obs:
        print('adding unified celltype broad to AD')
        adata_mapping.obs['Celltype_Class'] = adata_mapping.obs.apply(
            lambda row: sh.assign_celltype_fullnames(row['Celltype']), axis=1
        )
else:
    print('Mapping adata already has Celltype_Class')

subtype_fullname_key=adata_mapping.obs[['Celltype_Class', 'Subtype']].drop_duplicates()
subtype_to_class = dict(subtype_fullname_key[['Subtype', 'Celltype_Class']].values)

datasets = {'BD_SZ': [args.BD_SZ, args.BD_SZ_outfile], 'MS': [args.MS, args.MS_outfile]}

with open(args.report, "w") as f:
    f.write("TACCO Annotation Report\n")
    f.write("=======================\n\n")

for adata_name in datasets:
    print('starting tacco for', adata_name)
    adata_path=datasets[adata_name][0]
    outfile_path=datasets[adata_name][1]
    adata=sc.read_h5ad(adata_path)
    if adata.X.max() < 100: #then log normalized
        adata.layers['normalized_counts']=adata.X #saving normalized counts before resetting
        adata.X=adata.layers['counts'] #resetting to raw counts since TACCO requires it

    if args.mapping_annotation_key in adata.obs:
        adata.obs[f"{args.mapping_annotation_key}_old"]=adata.obs[args.mapping_annotation_key]
        adata.obs.drop(args.mapping_annotation_key, axis=1, inplace=True) #dropping because otherwise tacco will not provide obsm

    tc.tl.annotate(adata=adata, reference=adata_mapping, annotation_key=args.mapping_annotation_key, result_key='unified_subtype')

    try:
        adata.obs['unified_subtype']
    except:
        try:
            print('manually adding unified celltype to obs')
            unified_celltype_matrix = adata.obsm['unified_subtype']
            max_indices = np.argmax(unified_celltype_matrix, axis=1)
            column_names = unified_celltype_matrix.columns
            adata.obs['unified_subtype'] = column_names[max_indices]
        except:
            print('issue adding obs from obsm... trying again\n')
            try:
                unified_celltype_matrix = adata.obsm['unified_subtype']
                adata.obs['unified_subtype'] = unified_celltype_matrix.idxmax(axis=1)
            except:
                 print('second failure to add obs from obsm\n')

    ##################### Adding Categories for Celltype Columns:
    #unifying mapped annotations into a general broad celltype
    adata.obs['Celltype_Class'] = adata.obs['unified_subtype'].map(lambda x: subtype_to_class.get(x, np.nan)) #adds Celltype_Class for mapped annotations
    print('number of nans in celltype_class for mapped:',adata.obs.Celltype_Class.isna().sum(), 'note this should be 0\n')

    #unifying original annotations into a general broad celltype (for use in comparison metrics)
    adata.obs['celltype_broad'] = adata.obs.apply(
        lambda row: sh.assign_celltype_class(row['Celltype']), axis=1
    )

    ###################### Cluster Uniformity Metrics:
    uniformity_adata=adata[~(adata.obs['celltype_broad'] == 'nan')]
    uniformity_adata=uniformity_adata[~(uniformity_adata.obs['celltype_broad'].isna())]
    cells_not_uniform=adata.shape[0]-uniformity_adata.shape[0] #number of cells with NA in the original dataset's celltype column (therefore should not be included in uniformity score)

    labels_true = uniformity_adata.obs['celltype_broad'].values
    labels_pred = uniformity_adata.obs['Celltype_Class'].values

    ARI_score=sklearn.metrics.adjusted_rand_score(labels_true, labels_pred)
    NMI_score=sklearn.metrics.normalized_mutual_info_score(labels_true, labels_pred)

    ###################### Writing adata with new celltype columns:
    adata.write(outfile_path)

    ##################### Creating the Report:
    with open(args.report, "a") as f:
        f.write(f"{adata_name}\n")
        f.write(f"Results from using tacco to unify celltype annotations\n")
        f.write(f"ARI score for adata with mapped annotations: {ARI_score}\n")
        f.write(f"NMI score for adata with mapped annotations: {NMI_score}\n")
        f.write(f"For cluster uniformity scores, cells not included in ARI or NMI score due to NA celltype in original dataset: {cells_not_uniform}\n")

adata_mapping.obs['unified_subtype'] = adata_mapping.obs['Subtype'] #adding unified_celltype to AD as well so that when joining later on this column is not lost in others
adata_mapping.write(args.AD_outfile)