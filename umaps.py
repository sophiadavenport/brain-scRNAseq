import scanpy as sc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--joined_adata", required=True)
parser.add_argument("--celltype_umap_tracked", required=True)
args = parser.parse_args()

adata_combined=sc.read_h5ad(args.joined_adata, backed='r')
tracked=args.celltype_umap_tracked
save_celltype=tracked.split("/", 1)[1]

###UMAPs
sc.settings.figdir = "results/"

sc.pl.umap(adata_combined, color='unified_subtype', size=0.1, save=save_celltype, show=False, title='AD, BD, MS, & SZ Datasets: Subtype')
sc.pl.umap(adata_combined, color='Celltype_Class', size=0.1, save="ADBDMSSZ_full_unified_celltype__broad_umap.png",show=False, title='AD, BD, MS, & SZ Datasets: Cell Class')
sc.pl.umap(adata_combined, color='datasource', size=0.1, save="ADBDMSSZ_full_datasource_umap.png",show=False, title='AD, BD, MS, & SZ Datasets: Cell Class')

sc.pl.umap(adata_combined, color='unified_subtype', size=0.5, save="ADBDMSSZ_full_unified_celltype_umap_sizechange.png", show=False, title='AD, BD, MS, & SZ Datasets: Subtype')
sc.pl.umap(adata_combined, color='Celltype_Class', size=0.5, save="ADBDMSSZ_full_unified_celltype__broad_umap_sizechange.png",show=False, title='AD, BD, MS, & SZ Datasets: Cell Class')
sc.pl.umap(adata_combined, color='datasource', size=0.5, save="ADBDMSSZ_full_datasource_umap_sizechange.png",show=False, title='AD, BD, MS, & SZ Datasets: Cell Class')
