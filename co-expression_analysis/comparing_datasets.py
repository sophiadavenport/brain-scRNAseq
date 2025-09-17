import pandas as pd
import argparse
import numpy as np
import os

parser=argparse.ArgumentParser()
parser.add_argument("--annotated_csvs", nargs="+", required=True)
parser.add_argument("--outfile", required=True)
parser.add_argument("--excel_outfile", required=True)
args=parser.parse_args()

def celltypes_w_most_interactions(results):
    results=results.copy()
    results['celltype1']=results['module1'].str.split('_').str[0]
    results['celltype2']=results['module2'].str.split('_').str[0]

    new_df=results[['celltype1', 'celltype2']].drop_duplicates().copy()
    new_df['ligand_receptor_count']=None
    new_df['receptor_ligand_count']=None
    new_df['total']=None

    for i, row in new_df.iterrows():
        c1=row['celltype1']
        c2=row['celltype2']

        sub=results[(results['celltype1'] == c1) & (results['celltype2'] == c2)]
        lr_count=sub[sub['m1_ligand'].notna() & sub['m2_receptor'].notna()].shape[0]
        rl_count=sub[sub['m1_receptor'].notna() & sub['m2_ligand'].notna()].shape[0]
        total=lr_count + rl_count
        
        new_df.at[i, 'ligand_receptor_count']=lr_count if lr_count > 0 else None
        new_df.at[i, 'receptor_ligand_count']=rl_count if rl_count > 0 else None
        new_df.at[i, 'total']=total if total > 0 else None

    return new_df.sort_values(by=['total','ligand_receptor_count', 'receptor_ligand_count'], ascending=False)

def modules_w_most_interactions(results):
    new_df=results[['module1', 'module2']].drop_duplicates()
    new_df['ligand_receptor_count']=None
    new_df['receptor_ligand_count']=None
    new_df['total']=None
    for i,row in new_df.iterrows():
        m1=row['module1']
        m2=row['module2']

        sub=results[(results['module1']==m1) & (results['module2']==m2)]
        lr_count=sub[sub['m1_ligand'].notna() & sub['m2_receptor'].notna()].shape[0]
        rl_count=sub[sub['m1_receptor'].notna() & sub['m2_ligand'].notna()].shape[0]
        total=lr_count + rl_count
        
        new_df.at[i, 'ligand_receptor_count']=lr_count if lr_count > 0 else None
        new_df.at[i, 'receptor_ligand_count']=rl_count if rl_count > 0 else None
        new_df.at[i, 'total']=total if total > 0 else None

    return new_df.sort_values(by=['total','ligand_receptor_count', 'receptor_ligand_count'], ascending=False)

data=[]

writer=pd.ExcelWriter(args.excel_outfile, engine="xlsxwriter")
for path in args.annotated_csvs:
    dataset=path.split(os.sep)[2]
    cur_df=pd.read_csv(path)
    cur_df['datasource']=str(dataset)
    data.append(cur_df)

    celltype_df=celltypes_w_most_interactions(cur_df)
    module_df=modules_w_most_interactions(cur_df)
    celltype_df.to_excel(writer, sheet_name=f"{dataset}_celltype", index=False)
    module_df.to_excel(writer, sheet_name=f"{dataset}_module", index=False)
writer.close()

combined_data=pd.concat(data, ignore_index=True)
combined_data['celltype_m1']=combined_data['module1'].str.split('_').str[0]
combined_data['celltype_m2']=combined_data['module2'].str.split('_').str[0]

combined_data.to_csv(args.outfile, index=False)