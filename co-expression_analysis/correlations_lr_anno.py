import liana as li
import pandas as pd
import argparse
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("--all_modules", required=True)
parser.add_argument("--correlation_df", required=True)
parser.add_argument("--liana_csv", required=True)
parser.add_argument("--dataset", required=True)
parser.add_argument("--outfile", required=True)
parser.add_argument("--current_results_folder", required=True)
args=parser.parse_args()

correlations=pd.read_csv(args.correlation_df)
modules=pd.read_csv(args.all_modules)

def get_sig_df(cor_df, modules_df):
    sig_modules_df=modules_df[(modules_df['median_ks_pval']<0.05) & (modules_df['fisher_p_adj']<0.05)]
    sig_cof_df=cor_df[cor_df.padj<0.05]
    both_sigmods_sig_cof_df=sig_cof_df[(sig_cof_df.module1.isin(list(sig_modules_df.module))) & (sig_cof_df.module2.isin(list(sig_modules_df.module)))]
    print('significant modules shape:',sig_modules_df.shape, 'significiant correlations shape:', sig_cof_df.shape, 'both significant modules and correlations shape:',both_sigmods_sig_cof_df.shape)
    return(both_sigmods_sig_cof_df)

sig_correlations=get_sig_df(correlations, modules)

if sig_correlations.empty:
    results_cols=['m1_ligand','m2_receptor','m1_receptor','m2_ligand','module1','module2', 'liana_magnitude_rank', 'liana_lrscore', 'liana_cellphone_pvals']
    merged_cols=results_cols + list(correlations.columns)
    empty_out=pd.DataFrame(columns=merged_cols)
    empty_out.to_csv(args.outfile, index=False)
    print("No significant modules/correlations for {dataset}")
    exit()

liana_LRs=li.rs.select_resource('consensus')
lr_pairs=set(zip(liana_LRs['ligand'], liana_LRs['receptor']))
ligands=set(liana_LRs['ligand'])
receptors=set(liana_LRs['receptor'])

module2genes={}
all_modules=pd.concat([sig_correlations['module1'],sig_correlations['module2']])
modcelltypes=sorted(set([x.rsplit("_", 1)[0] for x in all_modules.unique()]))
for celltype in modcelltypes:
    module_genes=pd.read_csv(f'{args.current_results_folder}{args.dataset}__{celltype}/modules/modules.csv')
    module_genes_dict=(module_genes.groupby("module")["gene"].apply(list).to_dict())
    module_genes_dict={f"{celltype}_{k}": v for k, v in module_genes_dict.items()}
    module2genes=module2genes | module_genes_dict
    # module2genes.update(module_genes_dict)

results=[]

for row in sig_correlations.itertuples(index=False):
    module1=row.module1
    module2=row.module2
    m1_genes=module2genes[module1]
    m2_genes=module2genes[module2]
    
    for lig in set(m1_genes) & ligands: #Ligands in m1 to receptors in m2
        for rec in set(m2_genes) & receptors:
            if (lig, rec) in lr_pairs:
                results.append({'m1_ligand': lig,'m2_receptor': rec,'m1_receptor': 'NA','m2_ligand': 'NA','module1': module1,'module2': module2})
    for lig in set(m2_genes) & ligands: #Ligands in m2 to receptors in m1
        for rec in set(m1_genes) & receptors:
            if (lig, rec) in lr_pairs:
                results.append({'m1_ligand': 'NA','m2_receptor': 'NA','m1_receptor': rec,'m2_ligand': lig,'module1': module1,'module2': module2})

results_df=pd.DataFrame(results)

if results_df.empty:
    results_cols=['m1_ligand','m2_receptor','m1_receptor','m2_ligand','module1','module2', 'liana_magnitude_rank', 'liana_lrscore', 'liana_cellphone_pvals']
    merged_cols=results_cols + list(sig_correlations.columns)
    empty_out=pd.DataFrame(columns=merged_cols)
    empty_out.to_csv(args.outfile, index=False)
    print("No ligand–receptor pairs found for significant modules/correlations in {dataset}")
    exit()

sig_correlations_anno=pd.merge(results_df, sig_correlations, how='outer').sort_values(by='padj')

liana_resdf=pd.read_csv(args.liana_csv)

#comparing with liana results:
new_cols = ['liana_magnitude_rank', 'liana_lrscore', 'liana_cellphone_pvals']
for col in new_cols:
    if col not in sig_correlations_anno.columns:
        sig_correlations_anno[col] = np.nan
for i,row in sig_correlations_anno.iterrows():
    if (row['m1_ligand']!='NA') and pd.notna(row['m1_ligand']):
        ligand=row['m1_ligand']
        receptor=row['m2_receptor']
        result_source=row['module1'].rsplit("_", 1)[0].replace(' + ', ' ').replace(' ', '_')
        result_target=row['module2'].rsplit("_", 1)[0].replace(' + ', ' ').replace(' ', '_')
        liana_row=liana_resdf[(liana_resdf['source']==result_source) & (liana_resdf['target']==result_target) & (liana_resdf['ligand_complex']==ligand) & (liana_resdf['receptor_complex']==receptor)]
        if not liana_row.empty:
            sig_correlations_anno.at[i, 'liana_magnitude_rank'] = liana_row['magnitude_rank'].values[0]
            sig_correlations_anno.at[i, 'liana_lrscore'] = liana_row['lrscore'].values[0]
            sig_correlations_anno.at[i, 'liana_cellphone_pvals'] = liana_row['cellphone_pvals'].values[0]
    elif (row['m2_ligand']!='NA') and pd.notna(row['m2_ligand']):
        ligand=row['m2_ligand']
        receptor=row['m1_receptor']
        result_source=row['module2'].rsplit("_", 1)[0].replace(' + ', ' ').replace(' ', '_')
        result_target=row['module1'].rsplit("_", 1)[0].replace(' + ', ' ').replace(' ', '_')
        liana_row=liana_resdf[(liana_resdf['source']==result_source) & (liana_resdf['target']==result_target) & (liana_resdf['ligand_complex']==ligand) & (liana_resdf['receptor_complex']==receptor)]
        if not liana_row.empty:
            sig_correlations_anno.at[i, 'liana_magnitude_rank'] = liana_row['magnitude_rank'].values[0]
            sig_correlations_anno.at[i, 'liana_lrscore'] = liana_row['lrscore'].values[0]
            sig_correlations_anno.at[i, 'liana_cellphone_pvals'] = liana_row['cellphone_pvals'].values[0]

sig_correlations_anno.to_csv(args.outfile, index=False)