library(devtools)
library(dplyr)
if (!requireNamespace("TRADEtools", quietly=TRUE)) {
  message("TRADEtools not found, trying install...")
  devtools::install_github("immunogenomics/TRADEtools", dependencies=TRUE)
}

library(TRADEtools)

args <- commandArgs(trailingOnly=TRUE)
edger_res <- args[1]
output_rds <- args[2]
output_sig_genes <- args[3]

results <- read.csv(edger_res)

if (nrow(results)==0) {
  message("No rows in results")
  saveRDS(NULL, output_rds)
  write.csv(data.frame(), output_sig_genes, row.names = FALSE)
  quit(save="no", status=0)
}

results <- results |>
  dplyr::rename(
    log2FoldChange=logFC,
    pvalue=PValue
  )
results$F <- suppressWarnings(as.numeric(results$F))
results$F[is.na(results$F) | results$F < 0] <- 0 #ensure no NANs produced in lfcSE
results$lfcSE <- 1 / sqrt(results$F + 1e-8)  #crude proxy from F statistic since edgeR doesn't give lfcSE


TRADE_output <- TRADE(mode="univariate", results1=results, annot_table=NULL, genes_exclude=NULL, n_sample=NULL)
sig_genes <- data.frame(TRADE_output$significant_genes_FDR)

saveRDS(TRADE_output, file=output_rds)
write.csv(sig_genes, output_sig_genes, row.names=FALSE)
