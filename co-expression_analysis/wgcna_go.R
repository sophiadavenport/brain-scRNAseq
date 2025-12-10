suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly=TRUE)
coexpr_path <- args[1]
pval_path <- args[2]
genes_path <- args[3]
output_dir <- args[4]

if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

coexpr<-tryCatch({
  as.matrix(read.csv(coexpr_path, header=FALSE))
  }, error = function(e) {
    warning(paste("Error reading coexpr file:", coexpr_path, "->", e$message))
    matrix(nrow = 0, ncol = 0)
    })
pvals<-tryCatch({
  as.matrix(read.csv(pval_path, header = FALSE))
  }, error = function(e) {
    warning(paste("Error reading pvals file:", pval_path, "->", e$message))
    matrix(nrow = 0, ncol = 0)
    })
genes_selected<-tryCatch({
  readLines(genes_path)
  }, error = function(e) {
    warning(paste("Error reading genes file:", genes_path, "->", e$message))
    character(0)
    })
    
if (length(genes_selected)==0 || nrow(coexpr)==0 || nrow(pvals)==0) {
  warning("Empty CSCORE output detected. Writing empty placeholder files.")
  write.csv(data.frame(gene=character(0), module=integer(0)),
            file.path(output_dir, "modules.csv"), row.names=FALSE)
  go_headers <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue",
                  "p.adjust", "qvalue", "geneID", "Count")
  for (i in 1:1) {  #at least one placeholder
    write.csv(as.data.frame(matrix(ncol=length(go_headers), nrow=0,
                                   dimnames=list(NULL, go_headers))),
              file.path(output_dir, paste0("module", i, "_GO.csv")),
              row.names=FALSE)
  }
  quit(save="no")
}

cat("Genes included...", length(genes_selected))

if (dim(coexpr)[1]==length(genes_selected)) { 
  rownames(coexpr) <- colnames(coexpr) <- genes_selected 
  } else { 
    first_row <- coexpr[1,-1] #skip corner 
    first_col <- coexpr[-1, 1] 
    if (all(is.na(suppressWarnings(as.numeric(first_row)))) && all(is.na(suppressWarnings(as.numeric(first_col))))) { 
      col_genes <- coexpr[1,-1] 
      row_genes <- coexpr[-1, 1] 
      mat <- coexpr[-1, -1, drop = FALSE] 
      mat <- apply(mat, 2, function(x) suppressWarnings(as.numeric(x))) 
      mat <- as.matrix(mat) 
      rownames(mat) <- row_genes 
      colnames(mat) <- col_genes 
      coexpr <- mat 
      } else { 
        cat("cannot assign gene names for coexpr") 
        quit(save="no") 
        } 
} 
if (dim(pvals)[1]==length(genes_selected)) { 
  rownames(pvals) <- colnames(pvals) <- genes_selected 
  } else { 
    first_row <- pvals[1,-1] #skip corner 
    first_col <- pvals[-1, 1] 
    if (all(is.na(suppressWarnings(as.numeric(first_row)))) && all(is.na(suppressWarnings(as.numeric(first_col))))) { 
      col_genes <- pvals[1,-1] #skip corner since is "" 
      row_genes <- pvals[-1, 1] 
      mat <- pvals[-1, -1, drop = FALSE] 
      mat <- apply(mat, 2, function(x) suppressWarnings(as.numeric(x))) 
      mat <- as.matrix(mat) 
      rownames(mat) <- row_genes 
      colnames(mat) <- col_genes 
      pvals <- mat 
      } else { 
        cat("cannot assign gene names for pvals") 
        quit(save="no") 
        } 
}

#Adjust p-values (BH correction)
p_matrix_BH <- matrix(0, nrow=nrow(pvals), ncol=ncol(pvals))
p_matrix_BH[upper.tri(p_matrix_BH)] <- p.adjust(pvals[upper.tri(pvals)], method="BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

#Threshold co-expression matrix (non-significant coexpr to 0)
coexpr[p_matrix_BH > 0.05] <- 0
coexpr <- (coexpr + t(coexpr)) / 2 #ensure matrix consistency

#Compute adjacency and TOM
adj <- adjacency.fromSimilarity(abs(coexpr), power=1)
TOM <- TOMsimilarity(adj)
dissTOM <- 1 - TOM
rownames(dissTOM) <- colnames(dissTOM) <- genes_selected

hclust_dist <- hclust(as.dist(dissTOM), method="average")
memb <- cutreeDynamic(dendro=hclust_dist, distM=dissTOM, deepSplit=2, pamRespectsDendro=FALSE, minClusterSize=10)
names(memb) <- genes_selected

module_df <- data.frame(gene=genes_selected, module=memb)
write.csv(module_df, file.path(output_dir, "modules.csv"), row.names=FALSE)

module_list <- lapply(sort(unique(memb[memb > 0])), function(k) names(which(memb==k)))

universe <- genes_selected

for (i in seq_along(module_list)) {
  eg <- tryCatch({
    enrichGO(
      gene         =module_list[[i]],
      OrgDb        ="org.Hs.eg.db",
      keyType      ="SYMBOL",
      ont          ="ALL",
      universe     =universe,
      pAdjustMethod="BH",
      pvalueCutoff =0.1
    )
  }, error=function(e) NULL)

  if (!is.null(eg) && nrow(eg@result) > 0) {
    out_file <- file.path(output_dir, paste0("module", i, "_GO.csv"))
    write.csv(eg@result, out_file, row.names=FALSE)
  }
}