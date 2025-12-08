library(devtools)

if (!requireNamespace("CSCORE", quietly = TRUE)) {
  remotes::install_github("ChangSuBiostats/CS-CORE")
}

suppressPackageStartupMessages({
    library(CSCORE)
    library(Seurat)
    library(dplyr)
    library(Matrix)
})

args <- commandArgs(trailingOnly=TRUE)
dataset <- args[1]
genes <- args[2]
barcodes <- args[3]
counts <- args[4]
metadata <- args[5]
coexpr_output <- args[6]
teststats_output <- args[7]
p_vals_output <- args[8]
selected_genes_output <- args [9]

counts_mat <- readMM(counts)
counts_mat <- t(counts_mat)
if (sum(counts_mat)==0){
  cat("counts matrix is empty (all 0). Write empty output...")
  outdir <- dirname(coexpr_output)
    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }
  write.csv(data.frame(), coexpr_output, row.names = FALSE)
    write.csv(data.frame(), teststats_output, row.names = FALSE)
    write.csv(data.frame(), p_vals_output, row.names = FALSE)
    write.table(character(), selected_genes_output, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
    quit(save="no")
}
gene_names <- read.table(genes, header=FALSE, stringsAsFactors=FALSE)[,1]
barcode_names <- read.table(barcodes, header=FALSE, stringsAsFactors=FALSE)[,1]
rownames(counts_mat) <- gene_names
tryCatch({
    colnames(counts_mat) <- barcode_names
}, error = function(e) {
    barcode_df <- read.delim(barcodes, header = FALSE, stringsAsFactors = FALSE)
    barcode_df <- as.character(barcode_df[[1]])
    barcode_df <- trimws(barcode_df)

    colnames(counts_mat) <<- barcode_df
})

seurat_obj <- CreateSeuratObject(counts=counts_mat, project=dataset, min.cells=0, min.features=0)
assay_counts <- GetAssayData(seurat_obj, assay="RNA", layer="counts")


if (file.exists(metadata) && file.info(metadata)$size > 0) {
  cat("Reading metadata\n")
  meta <- read.csv(metadata, row.names=1, stringsAsFactors=FALSE)
  
  common_cells <- intersect(rownames(seurat_obj@meta.data), rownames(meta))
  seurat_obj <- subset(seurat_obj, cells=common_cells)
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, meta[common_cells, , drop=FALSE])
} else {
  cat("No metadata provided or file is empty\n")
}

seurat_obj$nCount_RNA <- colSums(Seurat::GetAssayData(seurat_obj, assay="RNA", layer="counts"))

frac_expressed <- rowSums(assay_counts > 0) / ncol(assay_counts)
expressed_genes <- names(frac_expressed[frac_expressed >= 0.10])

get_valid_covariates <- function(seurat_obj) {
    meta <- seurat_obj@meta.data
    cov_list <- c()

  if ("Sex" %in% colnames(meta)) {
    if (length(unique(meta$Sex)) > 1) {
      cov_list <- c(cov_list, "Sex")
    } 
  }

  if ("Age" %in% colnames(meta)) {
    if (length(unique(meta$Age)) > 1) {
      cov_list <- c(cov_list, "Age")
    } 
  }

  mito_col <- NULL
  if ("pct_mt" %in% colnames(meta)) {
    mito_col <- "pct_mt"
  } else if ("pct_mito" %in% colnames(meta)) {
    mito_col <- "pct_mito"
  } else if ("pct_counts_mt" %in% colnames(meta)) {
    mito_col <- "pct_counts_mt"
  }

  if (!is.null(mito_col)) {
    if (length(unique(meta[[mito_col]])) > 1) {
      cov_list <- c(cov_list, mito_col)
    } 
  }

  if (length(cov_list)==0) {
    message("No valid covariates found")
  }

  message("Valid Covariates: ", cov_list)
  return(cov_list)
}

valid_covariates <- get_valid_covariates(seurat_obj)

res <- CSCORE(seurat_obj, genes=expressed_genes, covariate_names=valid_covariates)

coexpr_mat <- res$est
pval_mat <- res$p_value
teststat_mat <- res$test_stat

outdir <- dirname(coexpr_output)
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

write.csv(as.matrix(coexpr_mat), coexpr_output, row.names=TRUE)
write.csv(as.matrix(pval_mat), p_vals_output, row.names=TRUE)
write.csv(as.matrix(teststat_mat), teststats_output, row.names=TRUE)
write.table(expressed_genes, selected_genes_output, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")