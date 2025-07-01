library(WGCNA)
library(dynamicTreeCut)
library(clusterProfiler)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)
coexpr_path <- args[1]
pval_path <- args[2]
genes_path <- args[3]
celltype_dir <- args[4]

if (!dir.exists(celltype_dir)) dir.create(celltype_dir, recursive = TRUE)

coexpr <- as.matrix(read.csv(coexpr_path, header = FALSE))
pvals  <- as.matrix(read.csv(pval_path,  header = FALSE))
genes_selected <- readLines(genes_path)

#Assign gene names to matrices
rownames(coexpr) <- colnames(coexpr) <- genes_selected
rownames(pvals)  <- colnames(pvals)  <- genes_selected

#Adjust p-values (BH correction)
p_matrix_BH <- matrix(0, nrow = nrow(pvals), ncol = ncol(pvals))
p_matrix_BH[upper.tri(p_matrix_BH)] <- p.adjust(pvals[upper.tri(pvals)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

#Threshold co-expression matrix (non-significant coexpr to 0)
coexpr[p_matrix_BH > 0.05] <- 0
coexpr <- (coexpr + t(coexpr)) / 2 #ensure matrix consistency

#Compute adjacency and TOM
adj <- adjacency.fromSimilarity(abs(coexpr), power = 1)
TOM <- TOMsimilarity(adj)
dissTOM <- 1 - TOM
rownames(dissTOM) <- colnames(dissTOM) <- genes_selected

hclust_dist <- hclust(as.dist(dissTOM), method = "average")
memb <- cutreeDynamic(dendro = hclust_dist, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 10)
names(memb) <- genes_selected

module_df <- data.frame(gene = genes_selected, module = memb)
write.csv(module_df, file.path(output_dir, "modules.csv"), row.names = FALSE)
#write.csv(module_df, sub("_coexpr.csv", "_modules.csv", coexpr_path), row.names = FALSE)

module_list <- lapply(sort(unique(memb[memb > 0])), function(k) names(which(memb == k)))

universe <- genes_selected

for (i in seq_along(module_list)) {
  eg <- tryCatch({
    enrichGO(
      gene          = module_list[[i]],
      OrgDb         = "org.Hs.eg.db",
      keyType       = "SYMBOL",
      ont           = "ALL",
      universe      = universe,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
  }, error = function(e) NULL)

  if (!is.null(eg) && nrow(eg@result) > 0) {
    #out_file <- sub("_coexpr.csv", paste0("_module", i, "_GO.csv"), coexpr_path)
    #write.csv(eg@result, out_file, row.names = FALSE)
    out_file <- file.path(output_dir, paste0("module", i, "_GO.csv"))
    write.csv(eg@result, out_file, row.names = FALSE)
  }
}