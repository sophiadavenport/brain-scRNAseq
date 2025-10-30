import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--input_files', nargs='+', required=True)
parser.add_argument('--dataset', required=True)
parser.add_argument('--report', required=True)
parser.add_argument('--figures_dir', required=True)
args=parser.parse_args()

os.makedirs(args.figures_dir, exist_ok=True)

def plot_logFC_heatmap(matrix, tag, output_path, cluster_rows=True, cluster_cols=True):
    num_genes=matrix.shape[0]
    if num_genes <= 5:
        fig_height=8
    elif num_genes <= 10:
        fig_height=6
    else:
        fig_height=min(0.25 * num_genes, 30)
    font_size=max(4, min(8, 12 - 0.03 * num_genes))

    g=sns.clustermap(matrix, cmap="bwr", center=0, row_cluster=cluster_rows, col_cluster=cluster_cols, cbar_kws={'shrink': 0.75}, linecolor='lightgray')
    
    g.ax_row_dendrogram.set_visible(True)
    g.ax_col_dendrogram.set_visible(True)
    g.ax_heatmap.set_xlabel('Cell Type', fontsize=10)
    g.ax_heatmap.set_ylabel('Gene', fontsize=10)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=65, fontsize=8)
    g.ax_heatmap.set_yticks(np.arange(matrix.shape[0]))
    g.ax_heatmap.set_yticklabels(matrix.index, fontsize=font_size)
    g.figure.set_size_inches(10, fig_height)
    g.figure.suptitle(f"{tag} Significant Genes", fontsize=14, y=1.05, ha='center')
    cbar_ax=g.cax
    cbar_ax.set_title('logFC', fontsize=10, pad=10, loc='center')
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

report_lines=[]
report_lines=[f"Summary Report for Dataset: {args.dataset}\n"]
report_lines.append("=" * 60 + "\n")

significant_genes=set()
celltype_logfc={}

for filepath in args.input_files:
    df=pd.read_csv(filepath)
    filename=os.path.basename(filepath)
    celltype=filename.replace(f"{args.dataset}_", "").replace("_edger_results.csv", "")
    if df.empty:
        report_lines.append(f"Cell Type: {celltype}")
        report_lines.append("  No significant genes (file was empty)\n")
        continue
    sig_df=df[(df["FDR"] < 0.05) & (df["logFC"].abs() >= 1)]
    sig_genes=sig_df["gene"].tolist()
    significant_genes.update(sig_genes)
    celltype_logfc[celltype]=sig_df.set_index("gene")["logFC"].to_dict()

    up=sig_df[sig_df["logFC"] >= 1].shape[0]
    down=sig_df[sig_df["logFC"] <= -1].shape[0]
    total=sig_df.shape[0]

    report_lines.append(f"Cell Type: {celltype}")
    report_lines.append(f"  Total significant genes (|logFC| >= 1 & FDR < 0.05): {total}")
    report_lines.append(f"  Upregulated: {up}")
    report_lines.append(f"  Downregulated: {down}\n")

#Create logFC heatmap matrix
sorted_genes=sorted(significant_genes)
sorted_celltypes=sorted(celltype_logfc.keys())
heatmap_matrix=pd.DataFrame(index=sorted_genes, columns=sorted_celltypes)

for celltype in sorted_celltypes:
    for gene in sorted_genes:
        logfc=celltype_logfc[celltype].get(gene, 0)  #use 0 if not significant in this celltype
        heatmap_matrix.loc[gene, celltype]=logfc

heatmap_matrix=heatmap_matrix.astype(float)

heatmap_path=os.path.join(args.figures_dir, f"{args.dataset}_logFC_heatmap.png")
plot_logFC_heatmap(heatmap_matrix, tag=args.dataset, output_path=heatmap_path)

with open(args.report, 'w') as f:
    f.write("\n".join(report_lines))