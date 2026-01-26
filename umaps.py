import scanpy as sc
import argparse
import matplotlib as plt
plt.rcParams['pdf.fonttype']=42
parser = argparse.ArgumentParser()
parser.add_argument("--joined_adata", required=True)
parser.add_argument("--celltype_umap_tracked", required=True)
args = parser.parse_args()

adata_combined=sc.read_h5ad(args.joined_adata, backed='r')
tracked=args.celltype_umap_tracked
save_celltype=tracked.split("/", 1)[1]

###UMAPs
sc.settings.figdir="results/"
adata_combined_nona=adata_combined[adata_combined.obs.unified_subtype!='NA']
# sc.pl.umap(adata_combined_nona, color='unified_subtype', size=0.1, save=save_celltype, show=False, title='AD, ASD, BD, MS, & SZ Datasets: Subtype')
# sc.pl.umap(adata_combined, color='Celltype_Class', size=0.1, save="ADASDBDMSSZ_full_unified_celltype__broad_umap.png",show=False, title='AD, ASD, BD, MS, & SZ Datasets: Cell Class')
# sc.pl.umap(adata_combined, color='datasource', size=0.1, save="ADASDBDMSSZ_full_datasource_umap.png",show=False, title='AD, ASD, BD, MS, & SZ Datasets')

# sc.pl.umap(adata_combined_nona, color='unified_subtype', size=0.5, save="ADASDBDMSSZ_full_unified_celltype_umap_sizechange.png", show=False, title='AD, ASD, BD, MS, & SZ Datasets: Subtype')
# sc.pl.umap(adata_combined, color='Celltype_Class', size=0.5, save="ADASDBDMSSZ_full_unified_celltype__broad_umap_sizechange.png",show=False, title='AD, ASD, BD, MS, & SZ Datasets: Cell Class')
# sc.pl.umap(adata_combined, color='datasource', size=0.5, save="ADASDBDMSSZ_full_datasource_umap_sizechange.png",show=False, title='AD, ASD, BD, MS, & SZ Datasets')

adata_combined_nona=