use_conda = True
configfile: "config.yaml"
OUTPUT_PATH=config["output_path"]
AD_H5AD=config["AD_h5ad"]
AD_METADATA=config["AD_metadata"]
BD_COUNTS=config["BD_counts"]
BD_GENENAMES=config["BD_genenames"]
BD_METADATA=config["BD_metadata"]
SZ_COUNTS=config["SZ_counts"]
SZ_GENENAMES=config["SZ_genenames"]
SZ_METADATA=config["SZ_metadata"]
MS_COUNTS=config["MS_counts"]
MS_GENENAMES=config["MS_genenames"]
MS_METADATA=config["MS_metadata"]

rule all:
    input:
        "intermediate/MS_postqc.h5ad", "results/common_gene_list.csv", "intermediate/AD_postqc.h5ad", "intermediate/BD_SZ_postqc.h5ad", "results/common_gene_list.csv", "results/annotations_qc_report.txt", f"{OUTPUT_PATH}/ADBDMSSZ_full.h5ad", "results/ADBDMSSZ_full_unified_celltype_umap.png"

rule filter_BD_SZ: #BD and SZ data must be filtered first since they have overlapping cells
    input:
       BD_counts=BD_COUNTS, BD_metadata=BD_METADATA, BD_genenames=BD_GENENAMES, SZ_counts=SZ_COUNTS, SZ_metadata=SZ_METADATA, SZ_genenames=SZ_GENENAMES
    output:
        BDSZ_intermediate=temp(f"intermediate/BDSZ_preqccheck.h5ad")
    conda:
        "envs/scanpy.yaml"
    shell:
        """
        python filter_BD_SZ.py --BD_counts {input.BD_counts} --BD_metadata {input.BD_metadata} --BD_genenames {input.BD_genenames} --SZ_counts {input.SZ_counts} --SZ_metadata {input.SZ_metadata} --SZ_genenames {input.SZ_genenames} --output {output.BDSZ_intermediate}
        """

rule QC_check: #runs a QC check but also will output a list of genes that will be kept in each dataset
    input:
        BD_SZ_h5ad=f"intermediate/BDSZ_preqccheck.h5ad", AD_h5ad=AD_H5AD, MS_counts=MS_COUNTS, MS_genenames=MS_GENENAMES, MS_metadata=MS_METADATA, AD_metadata=AD_METADATA
    output:
        keep_genes="results/common_gene_list.csv", MS_outfile="intermediate/MS_postqc.h5ad", AD_outfile="intermediate/AD_postqc.h5ad", BD_SZ_outfile="intermediate/BD_SZ_postqc.h5ad"
    conda:
        'envs/scanpy.yaml'
    shell:
        """
        python QC_check.py --BD_SZ {input.BD_SZ_h5ad} --AD_h5ad {input.AD_h5ad} --AD_metadata {input.AD_metadata} --MS_counts {input.MS_counts} --MS_genenames {input.MS_genenames} --MS_metadata {input.MS_metadata} --MS_outfile {output.MS_outfile} --AD_outfile {output.AD_outfile} --BD_SZ_outfile {output.BD_SZ_outfile} --keep_genes_file {output.keep_genes}
        """

rule tacco_annotations: #creating column unified_celltype and unified_celltype_broad with the mapped AD annotations
    input:
        MS=f"intermediate/MS_postqc.h5ad", AD=f"intermediate/AD_postqc.h5ad", BD_SZ=f"intermediate/BD_SZ_postqc.h5ad"
    output:
        MS_outfile=f"{OUTPUT_PATH}/MS_formatted.h5ad", AD_outfile=f"{OUTPUT_PATH}/AD_formatted.h5ad", BD_SZ_outfile=f"{OUTPUT_PATH}/BD_SZ_formatted.h5ad", annotations_report="results/annotations_qc_report.txt"
    params:
        mapping_annotation_key='Subtype'
    conda:
        'envs/TACCO_env.yaml'
    shell:
        """
        python unify_annotations.py --AD {input.AD} --BD_SZ {input.BD_SZ} --MS {input.MS} --mapping_annotation_key {params.mapping_annotation_key} --MS_outfile {output.MS_outfile} --AD_outfile {output.AD_outfile} --BD_SZ_outfile {output.BD_SZ_outfile} --report {output.annotations_report}
        """

rule join_adatas:
    input:
        MS=f"{OUTPUT_PATH}/MS_formatted.h5ad", AD=f"{OUTPUT_PATH}/AD_formatted.h5ad", BD_SZ=f"{OUTPUT_PATH}/BD_SZ_formatted.h5ad"
    output:
        joined_adata=f"{OUTPUT_PATH}/ADBDMSSZ_full.h5ad", umap="results/ADBDMSSZ_full_unified_celltype_umap.png"
    conda:
        "envs/scanpy.yaml"
    shell:
        """
        python --MS {input.MS} --AD {input.AD} --BD_SZ {input.BD_SZ} --joined_outfile {output.joined_adata} --umap_path {output.umap}
        """