library(tidyverse)
library(limma)
library(edgeR)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
barcodes_file <- args[2]
features_file <- args[3]
metadata_file <- args[4]
condition_col <- args[5]
output_file <- args[6]

counts <- readMM(counts_file)
print('dimensions of matrix:')
print(dim(counts))
print('Maximum counts:')
print(max(counts))
rownames(counts) <- readLines(features_file)
colnames(counts) <- readLines(barcodes_file)
metadata <- read.csv(metadata_file, row.names = 1)
metadata[[condition_col]] <- gsub(" ", ".", metadata[[condition_col]]) #ensure no spaces in conditions

group <- as.factor(metadata[[condition_col]])

print(levels(group))
print(table(group))

covariates <- metadata %>% group_by(.data[[condition_col]]) %>%
  summarise(
    num_cells = n(),
    avg_mt = mean(pct_counts_mt, na.rm = TRUE),
    avg_features = mean(n_genes_by_counts, na.rm = TRUE)
  )

print(covariates)

metadata <- metadata %>%
  left_join(covariates, by = condition_col)

dge <- DGEList(counts = counts, group = group)
keep <- filterByExpr(dge)
print('Sum of kept dge from filterbyexpr:')
sum(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]

#Recalculate library sizes and check
lib_sizes_filtered <- colSums(dge$counts)
print('summary of library sizes after filtering:')
summary(lib_sizes_filtered)
print('Any cells have total count of 0 after filtering:')
any(lib_sizes_filtered == 0)
print('Any cells with NA values after filtering:')
any(is.na(lib_sizes_filtered))

dge <- calcNormFactors(dge) 

initial_covariates <- c("avg_features", "avg_mt")
covariates <- initial_covariates
removed_covariates <- c()
design_success <- FALSE

while (!design_success && length(covariates) >= 0) {
  formula_str <- paste("~ 0 + group", if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else "")
  message("Trying design matrix with formula: ", formula_str)

  design <- model.matrix(as.formula(formula_str), data = metadata)
  colnames(design)[seq_along(levels(group))] <- levels(group)

  fit_try <- try({
    dge <- estimateDisp(dge, design)
    glmQLFit(dge, design)
  }, silent = TRUE)

  if (inherits(fit_try, "try-error")) {
    message("Model fitting failed with covariates: ", paste(covariates, collapse = ", "))
    if (length(covariates) == 0) {
      stop("All covariates removed and model still not estimable. Aborting.")
    }
    removed_covariate <- tail(covariates, 1)
    removed_covariates <- c(removed_covariates, removed_covariate)
    covariates <- head(covariates, -1)
    message("Removed covariate: ", removed_covariate)
  } else {
    fit <- fit_try
    design_success <- TRUE
    message("Model fitting succeeded with covariates: ", paste(covariates, collapse = ", "))
  }
}

if (length(removed_covariates) > 0) {
  message("Final model excludes these covariates due to rank deficiency: ", paste(removed_covariates, collapse = ", "))
}

group_levels <- levels(group)
pairwise_comparisons <- combn(group_levels, 2, simplify = FALSE)
all_results <- list()

for (pair in pairwise_comparisons) {
  contrast_string <- paste0(pair[1], " - ", pair[2])
  contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design)
  fit2 <- glmQLFTest(fit, contrast = contrast_matrix)
  table <- topTags(fit2, n = Inf)$table
  table$comparison <- paste0(pair[1], "_vs_", pair[2])
  table$gene <- rownames(table)
  all_results[[paste(pair, collapse = "_vs_")]] <- table
}

final_results <- bind_rows(all_results)
write.csv(final_results, output_file, row.names = TRUE)