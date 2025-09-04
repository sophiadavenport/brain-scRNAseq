import scanpy as sc
import pandas as pd
import argparse
import random
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--adata", required=True)
parser.add_argument("--dataset", required=True)
parser.add_argument("--current_celltype", required=True)
parser.add_argument("--celltype_col", required=True)
parser.add_argument("--donor_col", required=True)
parser.add_argument("--batch_col", required=True)
parser.add_argument("--counts_path", required=True)
parser.add_argument("--metadata_path", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.adata)

print('adata format:', adata)

#Clean names
adata.obs["Celltype_clean"] = adata.obs[args.celltype_col].astype(str).str.replace(r" \+ ", "_", regex=True).str.replace(r" ", "_", regex=True).str.replace(r"\(", "", regex=True).str.replace(r"\)", "", regex=True).str.replace("/", "-")

#Select only current celltype
cell_subset = adata[adata.obs['Celltype_clean'] == args.current_celltype,:]
if cell_subset.shape[0] == 0: #writing blank file if no cells present of this cell type
    print(f"No {args.current_celltype} present in {args.dataset}")
    empty_counts = pd.DataFrame()
    empty_metadata = pd.DataFrame()
    
    empty_counts.to_csv(args.counts_path)
    empty_metadata.to_csv(args.metadata_path)
    exit(0)

#Category Type
cell_subset.obs["Celltype_clean"] = cell_subset.obs["Celltype_clean"].astype("category")
cell_subset.obs["Condition"] = cell_subset.obs["Condition"].astype("category") #all adatas should already have Condition standardized
cell_subset.obs["Individual"] = cell_subset.obs[args.donor_col].astype("category") #what will be used to pseudobulk
cell_subset.obs["batch"] = cell_subset.obs[args.batch_col].astype("category") #what will be used to pseudobulk

#Use raw counts
if "counts" in cell_subset.layers:
    cell_subset.X = cell_subset.layers["counts"]
else:
    print(args.dataset, args.current_celltype, 'Counts not available in layers. Current maximum counts is:', cell_subset.X.max(), '\n')

#Metadata safety checks
if 'pct_counts_mt' not in cell_subset.obs.columns:
    if 'percent.mt' in cell_subset.obs.columns:
        cell_subset.obs.rename(columns={'percent.mt': 'pct_counts_mt'}, inplace=True)
    elif 'percent_mt' in cell_subset.obs.columns:
        cell_subset.obs.rename(columns={'percent_mt': 'pct_counts_mt'}, inplace=True)
    elif 'pct_mito' in cell_subset.obs.columns:
        cell_subset.obs.rename(columns={'pct_mito': 'pct_counts_mt'}, inplace=True)
    else:
        print(f'Warning: Could not find mitochondrial column in {args.dataset} for {args.celltype}\n')

if 'n_genes_by_counts' not in cell_subset.obs.columns:
    if 'nFeature_RNA' in cell_subset.obs.columns:
        cell_subset.obs.rename(columns={'nFeature_RNA': 'n_genes_by_counts'}, inplace=True)
    else:
        print(f'Warning: Could not find feature count column in {args.dataset} for {args.celltype}\n')

#Function for spliting data into pseudobulks:
def pseudobulk(celltype_adata, donor_col, keep_obs=[], replicates=2):
    pbs=[]
    for cur_donor in celltype_adata.obs[donor_col].unique():
        indv_cell_subset = celltype_adata[celltype_adata.obs[donor_col]==cur_donor].copy()
        n_cells = indv_cell_subset.n_obs
        if n_cells < 10:
            continue #skip individual if has less than 10 cells contributing to the pseudobulk
        indices = list(indv_cell_subset.obs_names)
        random.shuffle(indices)
        indices = np.array_split(np.array(indices), replicates)
        for i, pseudo_rep in enumerate(indices):
            X = np.array(indv_cell_subset[pseudo_rep].X.sum(axis=0)).reshape(1, -1)
            rep_adata = sc.AnnData(X=X, var=indv_cell_subset.var[[]].copy())
            rep_adata.obs_names = [f"{cur_donor}_{i}"] #new barcode is the donor ID with replicate number
            if keep_obs != []:
                rep_adata.obs[keep_obs] = pd.DataFrame({col: [indv_cell_subset.obs[col].iloc[0]] * rep_adata.n_obs for col in keep_obs}, index=rep_adata.obs_names) #saving the observations in keep_obs
            rep_adata.obs['replicate'] = i
            rep_adata.obs['n_cells_per_indv'] = n_cells
            pbs.append(rep_adata)
    if len(pbs)==0:
        print(args.current_celltype, "does not have enough cells per individual. Saving empty file./n")
        empty_counts = pd.DataFrame()
        empty_metadata = pd.DataFrame()
        empty_counts.to_csv(args.counts_path)
        empty_metadata.to_csv(args.metadata_path)
        exit(0)
    pb = sc.concat(pbs)
    #Writing pseudobulk 
    counts_df = pd.DataFrame(pb.X, columns = pb.var_names, index=pb.obs_names) #counts df for celltype 
    counts_df.T.to_csv(args.counts_path) #genes by cells for edgeR
    pb.obs.to_csv(args.metadata_path)

pseudobulk(celltype_adata=cell_subset, donor_col='Individual', keep_obs=['pct_counts_mt', 'n_genes_by_counts', 'Condition', 'batch', 'Individual', 'Age', 'Sex'], replicates=2)