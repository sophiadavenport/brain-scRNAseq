library(tidyverse, quiet=TRUE)
library(limma, quiet=TRUE)
library(edgeR, quiet=TRUE)

args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]
metadata_file <- args[2]
condition_col <- args[3]
output_file <- args[4]

#Safe reading for empty csvs from splitting adata
safe_read_csv <- function(path) {
  result <- try(read.csv(path, row.names=1, check.names=FALSE), silent=TRUE)
  if (inherits(result, "try-error") || nrow(result)==0 || ncol(result)==0 || length(colnames(result))==0) {
    message("File is unreadable or has no valid content: ", path)
    return(NULL)
  }
  return(result)
}

counts <- safe_read_csv(counts_file)
#Function for empty file if no cells:
write_empty_output <- function(out_path) {
  empty_df <- data.frame(logFC=numeric(0), logCPM=numeric(0), F=numeric(0), PValue=numeric(0), FDR=numeric(0), comparison=character(0), gene=character(0))
  write.csv(empty_df, out_path, row.names=TRUE)
  message("Empty input detected. Wrote blank results with header to: ", out_path)
  quit(save="no", status=0)
}

#Check if inputs are empty or invalid: exit with empty file if so
if (is.null(counts) || nrow(counts)==0 || ncol(counts)==0 ) {
  write_empty_output(output_file)
}

metadata <- read.csv(metadata_file, row.names=1)
metadata[[condition_col]] <- gsub(" ", ".", metadata[[condition_col]]) #ensure no spaces in conditions

#Checking number of individuals contributing to each group
indiv_groups_check <- function(metadata, group_col=condition_col, min_n=3) {
  counts_per_group <- table(metadata[[group_col]])
  all(counts_per_group >= min_n)
}
if (indiv_groups_check(metadata)) {
  #Proceed with DEG
} else {
  message("Skipping DEG: insufficient individuals per group")
  write_empty_output(output_file)
}

sample <- factor(metadata[["replicate"]]) #avoid pseudoreplication
group <- as.factor(metadata[[condition_col]])

print(levels(group))
print(table(group))

#Blank file save and no edgeR run if only one Condition found in cells
if (nlevels(group) < 2) {
  message("Only one group level found in condition column (", condition_col, "): ", levels(group))
  write_empty_output(output_file)
}

dge <- DGEList(counts=counts, group=group)
keep <- filterByExpr(dge)
print('Sum of kept dge from filterbyexpr:')
sum(keep)
dge <- dge[keep, , keep.lib.sizes=FALSE]

#Recalculate library sizes and check
lib_sizes_filtered <- colSums(dge$counts)
print('summary of library sizes after filtering:')
summary(lib_sizes_filtered)
print('Any cells have total count of 0 after filtering:')
any(lib_sizes_filtered==0)
print('Any cells with NA values after filtering:')
any(is.na(lib_sizes_filtered))

dge <- calcNormFactors(dge)

initial_covariates <- c("n_cells_per_indv", "n_genes_by_counts", 'batch')
#Check for "Sex"
if ("Sex" %in% colnames(metadata) && all(!is.na(metadata$Sex))) {
  sex_factor <- as.factor(metadata$Sex)
  if (nlevels(sex_factor) > 1) {
    initial_covariates <- c(initial_covariates, "Sex")
    message('"Sex" found in metadata, has more than one level, and added to covariates.')
  } else {
    message('"Sex" found in metadata but only one level present; not adding to covariates.')
  }
} else {
  message('"Sex" not found or contains NA values; not adding to covariates.')
}

#Check for "Age"
if ("Age" %in% colnames(metadata) && all(!is.na(metadata$Age))) {
  age_factor <- as.factor(metadata$Age)
  if (nlevels(age_factor) > 1) {
    initial_covariates <- c(initial_covariates, "Age")
    message('"Age" found in metadata, has more than one level, and added to covariates.')
  } else {
    message('"Age" found in metadata but only one level present; not adding to covariates.')
  }
} else {
  message('"Age" not found or contains NA values; not adding to covariates.')
}

covariates <- initial_covariates
removed_covariates <- c()
design_success <- FALSE

while (!design_success && length(covariates) >= 0) {
  formula_str <- paste("~ 0 + group", if (length(covariates) > 0) paste("+", paste(covariates, collapse=" + ")) else "")
  message("Trying design matrix with formula: ", formula_str)

  design <- model.matrix(as.formula(formula_str), data=metadata)
  colnames(design)[seq_along(levels(group))] <- levels(group)

  fit_try <- try({
    dge <- estimateDisp(dge, design)
    glmQLFit(dge, design)
  }, silent=TRUE)

  if (inherits(fit_try, "try-error")) {
    message("Model fitting failed with covariates: ", paste(covariates, collapse=", "))
    if (length(covariates)==0) {
      stop("All covariates removed and model still not estimable. Aborting.")
    }
    removed_covariate <- tail(covariates, 1)
    removed_covariates <- c(removed_covariates, removed_covariate)
    covariates <- head(covariates, -1)
    message("Removed covariate: ", removed_covariate)
  } else {
    fit <- fit_try
    design_success <- TRUE
    message("Model fitting succeeded with covariates: ", paste(covariates, collapse=", "))
  }
}

if (length(removed_covariates) > 0) {
  message("Final model excludes these covariates due to rank deficiency: ", paste(removed_covariates, collapse=", "))
}

group_levels <- levels(group)
pairwise_comparisons <- combn(group_levels, 2, simplify=FALSE)
all_results <- list()

for (pair in pairwise_comparisons) {
  contrast_string <- paste0(pair[1], " - ", pair[2])
  contrast_matrix <- makeContrasts(contrasts=contrast_string, levels=design)
  fit2 <- glmQLFTest(fit, contrast=contrast_matrix)
  table <- topTags(fit2, n=Inf)$table
  table$comparison <- paste0(pair[1], "_vs_", pair[2])
  table$gene <- rownames(table)
  all_results[[paste(pair, collapse="_vs_")]] <- table
}

final_results <- bind_rows(all_results)
write.csv(final_results, output_file, row.names=TRUE)


#Volcano Plot (not tracked):
volcano_plot_path <- sub("([^/]+)$", "Volcano_\\1", output_file)
volcano_plot_path <- sub("\\.csv$", ".png", volcano_plot_path)

filename <- basename(output_file)
celltype_name <- sub(".*__(.*)_edger_results.*", "\\1", filename)
celltype_name <- gsub("_", " ", celltype_name)

fdr_threshold <- 0.05
logfc_threshold <- 1

if (length(all_results) > 0) {
  volcano_data <- all_results[[1]]
  
  p <- ggplot(volcano_data, aes(x=logFC, y=- log10(FDR))) +
  geom_point(aes(color=FDR < 0.05), alpha=0.6, size=0.5) +
  geom_vline(xintercept=c(-logfc_threshold, logfc_threshold), linetype="dashed", color="#484747") +
    geom_hline(yintercept=-log10(fdr_threshold), linetype="longdash", color="#2f2fc4") +
    labs(
      title=paste(celltype_name, ":", names(all_results)[1]),
      x="Log Fold Change",
      y="-Log10(FDR p-value)"
    ) +
    theme_light()
  
  ggsave(volcano_plot_path, plot=p, width=7, height=6)
}