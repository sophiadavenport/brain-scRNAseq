use_conda = True
configfile: "config.yaml"
OUTPUT_PATH=config["output_path"]
AD_H5AD=config["AD_h5ad"]
AD_METADATA=config["AD_metadata"]
BD_H5AD=config["BD_h5ad"]
SZ_COUNTS=config["SZ_counts"]
SZ_GENENAMES=config["SZ_genenames"]
SZ_METADATA=config["SZ_metadata"]
MS_COUNTS=config["MS_counts"]
MS_GENENAMES=config["MS_genenames"]
MS_METADATA=config["MS_metadata"]
ASD_H5AD=config["ASD_h5ad"]

rule all:
    input:
        "results/ADASDBDMSSZ_full_unified_celltype_umap.pdf", "results/common_gene_list.csv", "results/annotations_qc_report.txt", f"{OUTPUT_PATH}/ADASDBDMSSZ_processed.h5ad"

rule filter_BD_SZ: #BD and SZ data must be filtered first since they have overlapping cells
    input:
       BD_adata=BD_H5AD, SZ_counts=SZ_COUNTS, SZ_metadata=SZ_METADATA, SZ_genenames=SZ_GENENAMES
    output:
        BDSZ_intermediate=temp("intermediate/BDSZ_preqccheck.h5ad")
    conda:
        "envs/scanpy.yaml"
    shell:
        """
        python filter_BD_SZ.py --BD_adata {input.BD_adata} --SZ_counts {input.SZ_counts} --SZ_metadata {input.SZ_metadata} --SZ_genenames {input.SZ_genenames} --output {output.BDSZ_intermediate}
        """

rule QC_check: #runs a QC check but also will output a list of genes that will be kept in each dataset
    input:
        BD_SZ_h5ad=f"intermediate/BDSZ_preqccheck.h5ad", AD_h5ad=AD_H5AD, MS_counts=MS_COUNTS, MS_genenames=MS_GENENAMES, MS_metadata=MS_METADATA, AD_metadata=AD_METADATA, ASD_h5ad=ASD_H5AD
    output:
        keep_genes="results/common_gene_list.csv", ASD_outfile="intermediate/ASD_postqc.h5ad", MS_outfile="intermediate/MS_postqc.h5ad", AD_outfile="intermediate/AD_postqc.h5ad", BD_SZ_outfile="intermediate/BD_SZ_postqc.h5ad"
    conda:
        'envs/scanpy.yaml'
    shell:
        """
        python QC_check.py --BD_SZ {input.BD_SZ_h5ad} --AD_h5ad {input.AD_h5ad} --AD_metadata {input.AD_metadata} --MS_counts {input.MS_counts} --MS_genenames {input.MS_genenames} --MS_metadata {input.MS_metadata} --MS_outfile {output.MS_outfile} --AD_outfile {output.AD_outfile} --BD_SZ_outfile {output.BD_SZ_outfile} --ASD_h5ad {input.ASD_h5ad} --ASD_outfile {output.ASD_outfile} --keep_genes_file {output.keep_genes}
        """

rule tacco_annotations: #creating column unified_celltype and unified_celltype_broad with the mapped AD annotations
    input:
        MS=f"intermediate/MS_postqc.h5ad", AD=f"intermediate/AD_postqc.h5ad", BD_SZ=f"intermediate/BD_SZ_postqc.h5ad", ASD="intermediate/ASD_postqc.h5ad"
    output:
        MS_outfile=f"{OUTPUT_PATH}/MS_formatted.h5ad", AD_outfile=f"{OUTPUT_PATH}/AD_formatted.h5ad", BD_outfile=f"{OUTPUT_PATH}/BD_formatted.h5ad", SZ_outfile=f"{OUTPUT_PATH}/SZ_formatted.h5ad", ASD_outfile=f"{OUTPUT_PATH}/ASD_formatted.h5ad", annotations_report="results/annotations_qc_report.txt"
    params:
        mapping_annotation_key='Subtype'
    conda:
        'envs/TACCO_env.yaml'
    shell:
        """
        python unify_annotations.py --AD {input.AD} --BD_SZ {input.BD_SZ} --MS {input.MS} --ASD {input.ASD} --mapping_annotation_key {params.mapping_annotation_key} --MS_outfile {output.MS_outfile} --AD_outfile {output.AD_outfile} --BD_outfile {output.BD_outfile} --SZ_outfile {output.SZ_outfile} --ASD_outfile {output.ASD_outfile} --report {output.annotations_report}
        """

rule join_adatas:
    input:
         MS=f"{OUTPUT_PATH}/MS_formatted.h5ad", AD=f"{OUTPUT_PATH}/AD_formatted.h5ad", BD=f"{OUTPUT_PATH}/BD_formatted.h5ad", SZ=f"{OUTPUT_PATH}/SZ_formatted.h5ad", ASD=f"{OUTPUT_PATH}/ASD_formatted.h5ad"
    output:
        joined_adata=temp(f"{OUTPUT_PATH}/ADASDBDMSSZ_full.h5ad")
    conda:
         "envs/scanpy.yaml"
    shell:
         """
         python join_adatas.py --MS {input.MS} --AD {input.AD} --BD {input.BD} --SZ {input.SZ} --ASD {input.ASD} --joined_outfile {output.joined_adata}
         """

rule redo_clustering:
    input:
        joined_adata=f"{OUTPUT_PATH}/ADASDBDMSSZ_full.h5ad"
    output:
        joined_outfile=f"{OUTPUT_PATH}/ADASDBDMSSZ_processed.h5ad"
    conda:
        "envs/scanpy.yaml"
    shell:
        """
        python data_cluster.py --joined_adata {input.joined_adata} --joined_outfile {output.joined_outfile}
        """

rule umaps:
    input:
        joined_adata=f"{OUTPUT_PATH}/ADASDBDMSSZ_processed.h5ad"
    output:
        celltype_umap_tracked="results/ADASDBDMSSZ_full_unified_celltype_umap.pdf"
    conda:
        "envs/scanpy.yaml"
    shell:
        """
        python umaps.py --joined_adata {input.joined_adata} --celltype_umap_tracked {output.celltype_umap_tracked}
        """