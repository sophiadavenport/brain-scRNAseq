import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--celltype", required=True)
parser.add_argument("--dataset", required=True)
parser.add_argument("--modules_csv_path", required=True)
parser.add_argument("--de_folder", required=True)
parser.add_argument("--ks_csv_path", required=True)
parser.add_argument("--avg_logfc_csv_path", required=True)
parser.add_argument("--number_dist", type=int, required=True)
args = parser.parse_args()

de_csv_path=f"{args.de_folder}/{args.dataset}_DEG/results/{args.dataset}_celltypeclass__{args.celltype}_edger_results.csv"

def logfc_distribution_per_module(celltype, dataset, modules_csv_path, de_csv_path, ks_csv_path, avg_logfc_csv_path, number_dist=50):
    #requires: scipy.stats as stats, os, matplotlib as plt, seaborn as sns, numpy as np, and pandas as pd
    modules_df = pd.read_csv(modules_csv_path)
    de_df = pd.read_csv(de_csv_path)

    if 'gene' not in modules_df.columns or 'module' not in modules_df.columns:
        print("modules.csv must contain 'gene' and 'module' columns.")
        return
    if 'gene' not in de_df.columns or 'logFC' not in de_df.columns:
        print("DE file must contain 'gene' and 'logFC' columns.")
        return

    de_df['abs_logFC'] = de_df['logFC'].abs()
    merged_df = pd.merge(modules_df, de_df[['gene', 'abs_logFC']], on='gene', how='inner')
    background_logfc = de_df['abs_logFC'].dropna()
    avg_logfc_df = merged_df.groupby('module')['abs_logFC'].mean().reset_index()
    avg_logfc_df.rename(columns={'abs_logFC': 'average_abslogFC'}, inplace=True)

    #Collect GO_BP terms with p.adjust < 0.05
    go_base_folder=modules_csv_path.split("/modules/")[0] + "/modules"
    go_terms = []
    for module_id in avg_logfc_df['module']:
        go_file = os.path.join(go_base_folder, f"module{module_id}_GO.csv")
        go_bp_list = []

        if os.path.isfile(go_file):
            try:
                go_df = pd.read_csv(go_file)
                if {'ONTOLOGY', 'Description', 'p.adjust'}.issubset(go_df.columns):
                    bp_df = go_df[go_df['ONTOLOGY'] == 'BP']
                    bp_df = bp_df[bp_df['p.adjust'] < 0.05]
                    bp_df_sorted = bp_df.sort_values(by='p.adjust', ascending=True)
                    go_bp_list = bp_df_sorted['Description'].dropna().tolist()
            except Exception as e:
                print(f"Error reading GO file for module {module_id}: {e}")
        else:
            print(f"Missing GO file for module {module_id}: {go_file}")

        go_terms.append(go_bp_list)

    avg_logfc_df['GO_BP'] = go_terms

    os.makedirs(os.path.dirname(avg_logfc_csv_path), exist_ok=True)
    avg_logfc_df.to_csv(avg_logfc_csv_path, index=False)

    module_groups = list(merged_df.groupby('module'))

    ks_results = []
    plot_data = {}

    for module_id, group in module_groups:
        module_genes = group['gene'].dropna().unique()
        module_logfc = group['abs_logFC'].dropna()  
        if module_logfc.empty:
            print(f"Module {module_id}: No logFC values found.")
            continue
        background_pool = de_df[~de_df['gene'].isin(module_genes)]['abs_logFC'].dropna() #Pool of background genes (excluding module)
        if len(background_pool) < len(module_logfc):
            print(f"Module {module_id}: Not enough background genes to sample.")
            continue
    
        sampled_dists = []
        p_values = []
        for _ in range(int(number_dist)): #Sample background distributions
            sample = background_pool.sample(n=len(module_logfc), replace=False)
            sampled_dists.append(sample)
            _, p = stats.ks_2samp(module_logfc, sample, alternative='two-sided', method='auto')
            p_values.append(p)
        #For plotting: Choose the background distribution whose median is closest to the median of medians
        medians = [s.median() for s in sampled_dists]
        median_of_medians = np.median(medians)
        median_index = np.argmin([abs(m - median_of_medians) for m in medians])
        chosen_background = sampled_dists[median_index]

        median_p = np.median(p_values)
        confidence = sum(p < 0.05 for p in p_values) / len(p_values)
        ks_results.append({'module': module_id,'median_ks_pval': median_p,'confidence': confidence,'ks_p_values': p_values})
        plot_data[module_id] = {'module_logfc': module_logfc,'background_logfc': chosen_background}

    ks_df=pd.DataFrame(ks_results)
    
    if not ks_df.empty:
        ks_df["ks_p_values"] = ks_df["ks_p_values"].apply(lambda x: [round(float(p), 5) for p in x])
    else:
        print(f"No valid KS-test results for {dataset} {celltype}")
    os.makedirs(os.path.dirname(ks_csv_path), exist_ok=True)
    ks_df.to_csv(ks_csv_path, index=False)

    num_modules = len(plot_data)
    if num_modules > 0:
        ncols = 3
        nrows = (num_modules + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 4, nrows * 3))
        axes = axes.flatten()

        for idx, (module_id, data) in enumerate(plot_data.items()):
            ax = axes[idx]
            sns.histplot(data['background_logfc'], bins=50, color='lightgray', alpha=0.5, kde=False, label='Random Genes', ax=ax)
            sns.histplot(data['module_logfc'], bins=50, color='steelblue', alpha=0.8, kde=False, label='Module Genes', ax=ax)
            ax.set_title(f"Module {module_id}")
            ax.set_xlabel('|logFC|')
            ax.set_ylabel('Number of Genes')
            ax.legend()

        for ax in axes[num_modules:]:
            ax.axis('off')

        plt.suptitle(f"{dataset} {celltype}: logFC Distribution per Module", fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plot_dir = os.path.dirname(ks_csv_path)
        panel_path = os.path.join(plot_dir, f"{dataset}_{celltype}_logFC_panel.png")
        plt.savefig(panel_path, dpi=300)
        plt.close()
        print(f"Saved panel figure: {panel_path}")

        # Save per-module plots
        for module_id, data in plot_data.items():
            plt.figure(figsize=(8, 5))
            sns.histplot(data['background_logfc'], bins=50, color='lightgray', alpha=0.5, kde=False, label='Random Genes')
            sns.histplot(data['module_logfc'], bins=50, color='steelblue', alpha=0.8, kde=False, label='Module Genes')
            plt.title(f"{dataset} {celltype}: Module {module_id}")
            plt.xlabel('|logFC|')
            plt.ylabel('Number of Genes')
            plt.legend()
            plt.tight_layout()
            filename = f"module{module_id}_{dataset}_{celltype}_logFCdist.png"
            save_path = os.path.join(plot_dir, filename)
            plt.savefig(save_path, dpi=300)
            plt.close()
    else:
        print(f"No valid modules for plotting in {dataset} {celltype}")


    if not ks_df.empty and "median_ks_pval" in ks_df.columns:
        significant_modules = ks_df[ks_df['median_ks_pval'] < 0.05]['module'].tolist()
        print(f"Significant modules ({len(significant_modules)}): {significant_modules}")
    else:
        significant_modules = []
        print(f"No significant modules for {dataset} {celltype}")

if not os.path.isfile(de_csv_path):
    print(f"DE file not found: {de_csv_path}. Skipping computation for {args.dataset} - {args.celltype}.")
    #creating blank files so rule all completes successfully
    pd.DataFrame(columns=["module", "KS_stat", "p_value"]).to_csv(args.ks_csv_path, index=False)
    pd.DataFrame(columns=["module", "avg_abs_logFC"]).to_csv(args.avg_logfc_csv_path, index=False)
    sys.exit(0)

logfc_distribution_per_module(celltype=args.celltype, dataset=args.dataset, modules_csv_path=args.modules_csv_path,de_csv_path=de_csv_path,
                              ks_csv_path=args.ks_csv_path, avg_logfc_csv_path=args.avg_logfc_csv_path,number_dist=args.number_dist)