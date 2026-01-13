import scanpy as sc
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import argparse
import sys
import os
from statsmodels.stats.multitest import multipletests

#Multivariable linear regression of module activity on condition with covariate adjustment or covariate adjusted linear model
parser=argparse.ArgumentParser()
parser.add_argument("--adata", required=True)
parser.add_argument("--celltype", required=True)
parser.add_argument("--celltype_col", required=True)
parser.add_argument("--indv_matrix", required=True)
parser.add_argument("--modules_csv_path", required=True)
parser.add_argument("--output_file", required=True)
args = parser.parse_args()

adata=sc.read_h5ad(args.adata, backed='r')
modules=pd.read_csv(args.modules_csv_path)
if modules.shape[0]==0:
    print('no modules for ', args.celltype)
    results_df=pd.DataFrame(data=None, columns=["module","coef","std_err","t","p_value","ci_low","ci_high","n"])
    results_df.to_csv(args.output_file) #writing empty output if no modules
    sys.exit(0)

indv_matrix=pd.read_csv(args.indv_matrix)
indv_matrix=indv_matrix.set_index('Unnamed: 0')

if 'individual_id_anon' in list(adata.obs.columns):
    indv_col='individual_id_anon'
elif 'Individual_ID' in list(adata.obs.columns):
    indv_col='Individual_ID'
elif 'individual_id' in list(adata.obs.columns):
    indv_col='individual_id'
elif 'Individual' in list(adata.obs.columns):
    indv_col='Individual'
elif 'individualID' in list(adata.obs.columns):
    indv_col='individualID'
else:
    print('individual column not found')
    sys.exit(1)

if 'pmi_cat' in list(adata.obs.columns):
    adata.obs.rename(columns={'pmi_cat': 'PMI'}, inplace=True)
    pmi_col='PMI'
elif 'PMI' in list(adata.obs.columns):
    pmi_col='PMI'
else:
    pmi_col=None

if 'batch' in list(adata.obs.columns):
    adata.obs.rename(columns={'batch': 'Batch'}, inplace=True)
    batch_col='Batch'
elif 'Batch' in list(adata.obs.columns):
    batch_col='Batch'
else:
    batch_col=None

meta_cols=[indv_col, 'Sex', 'Age', 'Condition']
default_formula="module_activity ~ Condition + Age + sex + avg_pct_mito + num_cells"
if pmi_col!=None:
    meta_cols.append(pmi_col)
    default_formula=default_formula+" + PMI"
if batch_col!=None:
    meta_cols.append(batch_col)
    default_formula=default_formula+" + C(Batch)" #C() keeps batch as category

meta=adata.obs[meta_cols].drop_duplicates().reset_index()[meta_cols] #indv_col, 'Sex', 'Age', 'Condition'
if meta[indv_col].duplicated().any():
    print(f"Warning: {meta[indv_col].duplicated().sum()} duplicate rows found for {indv_col}.") #drops repeats by keeping first
    meta=meta.drop_duplicates(subset=[indv_col], keep='first').reset_index(drop=True)

if 'pct_mito' not in list(adata.obs.columns):
    if 'percent.mt' in list(adata.obs.columns):
        adata.obs.rename(columns={'percent.mt': 'pct_mito'}, inplace=True)
    elif 'percent_mt' in list(adata.obs.columns):
        adata.obs.rename(columns={'percent_mt': 'pct_mito'}, inplace=True)
    elif 'pct_counts_mt' in list(adata.obs.columns):
        adata.obs.rename(columns={'pct_counts_mt': 'pct_mito'}, inplace=True)
    else:
        print(f'Warning: Could not find mitochondrial column in {args.dataset} for {args.celltype}\n')

avg_pct_mito=(
    adata.obs
      .groupby(indv_col)['pct_mito']
      .mean()
      .rename("avg_pct_mito")
)
sub_cell=adata[adata.obs[args.celltype_col]==args.celltype]
num_cells=sub_cell.obs[indv_col].value_counts().rename('num_cells')

meta=meta.merge(
    avg_pct_mito,
    left_on=indv_col,
    right_index=True,
    how="left"
)
meta = meta.merge(
    num_cells,
    left_on=indv_col,
    right_index=True,
    how="left"
)
meta['num_cells']=meta['num_cells'].fillna(0)

def run_module_regression(module_id, module_matrix, meta_df, indv_expr, indv_col, formula="module_activity ~ Condition + Age + sex + avg_pct_mito + num_cells"):
    indv_z=indv_expr.sub(indv_expr.mean(axis=1), axis=0)
    indv_z=indv_z.div(indv_expr.std(axis=1), axis=0)

    meta_df["sex"]=meta_df["Sex"].map({"F": 0, "M": 1})

    df=(indv_z.reset_index().rename(columns={"Unnamed: 0": "module"}).melt(
        id_vars="module",
        var_name="individual_id",
        value_name="module_activity"
    )
       )
    meta=meta_df.rename(columns={indv_col: "individual_id"})
    df=df.merge(meta, on="individual_id", how="left")
    df.groupby("module")["individual_id"].nunique().head()
    df.groupby("module")["module_activity"].std().describe()

    subdf=df[df["module"] == module_id]
    
    model=smf.ols(formula, data=subdf)
    res=model.fit()

    cond_term=[c for c in res.params.index if c.startswith("Condition")][0]

    return {
        "module": module_id,
        "coef": res.params[cond_term],
        "std_err": res.bse[cond_term],
        "t": res.tvalues[cond_term],
        "p_value": res.pvalues[cond_term],
        "ci_low": res.conf_int().loc[cond_term, 0],
        "ci_high": res.conf_int().loc[cond_term, 1],
        "n": int(res.nobs)
    }

print(f'Covariate formula used: {default_formula}')
results=[]

for module_id in indv_matrix.index:
    out=run_module_regression(module_id=module_id, module_matrix=modules, meta_df=meta, indv_expr=indv_matrix, indv_col=indv_col, formula=default_formula)
    if out is not None:
        results.append(out)

res_df=pd.DataFrame(results)
if not res_df.empty and "p_value" in res_df.columns: #FDR correction across modules
    res_df["p_adj"]=multipletests(res_df["p_value"], method="fdr_bh")[1]
else:
    print('fdr adjustment for across modules failed')
    res_df["p_adj"]=[]

res_df.to_csv(args.output_file, index=False)