import scanpy as sc
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import argparse
import sys
import os

parser=argparse.ArgumentParser()
parser.add_argument("--adata", required=True)
parser.add_argument("--celltype", required=True)
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

if 'individual_id_anon' in adata.obs:
    indv_col='individual_id_anon'
elif 'Individual_ID' in adata.obs:
    indv_col='Individual_ID'
elif 'individual_id' in adata.obs:
    indv_col='individual_id'

meta=adata.obs[[indv_col, 'Sex', 'Age', 'Condition']].drop_duplicates().reset_index()[[indv_col, 'Sex', 'Age', 'Condition']]

avg_pct_mito=(
    adata.obs
      .groupby(indv_col)['pct_mito']
      .mean()
      .rename("avg_pct_mito")
)
sub_cell=adata[adata.obs['unified_subtype']==args.celltype]
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

def run_module_regression(module_id, module_matrix, meta_df, indv_expr, indv_col,
                          formula="module_activity ~ Condition + Age + sex + avg_pct_mito + num_cells"):

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

results=[]

for module_id in indv_matrix.index:
    out=run_module_regression(module_id=module_id, module_matrix=modules, meta_df=meta, indv_expr=indv_matrix, indv_col=indv_col)
    if out is not None:
        results.append(out)

res_df=pd.DataFrame(results)

res_df.to_csv(args.output_file)