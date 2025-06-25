import pandas as pd


def get_sig_genes(file_path, logFC_thres=1, pval_sig=0.05):
    '''
    Purpose is to take an input csv from run_edgeR's outputs and return a list of significant genes from that file.
    Requires: pandas as pd
    Input:
        - file_path: full path to .csv file (each row represents a gene; should have header with: "","logFC","logCPM","F","PValue","FDR","comparison","gene")
        - logFC_thres: threshold for fold-change (will automatically be set to abs-value greater than or equal to)
        - pval_sig: significance threshold for p-value (will be set to less than)
    Output:
        - list of significant genes
    '''
    df = pd.read_csv(file_path)

    required_cols = {"logFC", "PValue", "gene"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required_cols - set(df.columns)}")

    sig_df = df[(df["PValue"] < pval_sig) & (df["logFC"].abs() >= logFC_thres)]

    return sig_df["gene"].tolist()