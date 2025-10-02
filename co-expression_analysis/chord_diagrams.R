library(CCPlotR)
library(dplyr)
library(stringr)
library(circlize)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
sig_correlations_annotated_path <- args[1]
dataset_name <- args[2]
celltype_column <- args[3]
pdf_plot_path <- args[4]
png_plot_path <- args[5]

correlations <- read.csv(sig_correlations_annotated_path)

prep_cc_data <- function(df) {
  clean_module <- function(x) {
    x %>%
    str_replace_all("_", " ") %>%   
    str_remove_all("[0-9]") %>%     
    str_trim()
    }

  filtered_df <- df %>% filter (!(m1_ligand=="" & m2_ligand=="" & m1_receptor =="" & m2_receptor=="")) %>% 
    filter(!(is.na(m1_ligand) & is.na(m2_ligand) & is.na(m1_receptor) & is.na(m2_receptor)))

  return_df <- filtered_df %>%
    mutate(
      case_a = !(is.na(m1_ligand)) & !(is.na(m2_receptor)) & is.na(m2_ligand) & is.na(m1_receptor),
      case_b = !(is.na(m2_ligand)) & !(is.na(m1_receptor)) & is.na(m1_ligand) & is.na(m2_receptor),
    ligand = case_when(
      case_a ~ m1_ligand,
      case_b ~ m2_ligand,
      TRUE ~ NA_character_),
    receptor = case_when(
      case_a ~ m2_receptor,
      case_b ~ m1_receptor,
      TRUE ~ NA_character_),
    source_celltype = case_when(
      case_a ~ clean_module(module1),
      case_b ~ clean_module(module2),
      TRUE ~ NA_character_),
    receiver_celltype = case_when(
      case_a ~ clean_module(module2),
      case_b ~ clean_module(module1),
      TRUE ~ NA_character_)
    ) %>% 
    mutate(
      source   = ifelse(source_celltype   == "OPCs COPs", "OPCs + COPs", source_celltype),
      target = ifelse(receiver_celltype == "OPCs COPs", "OPCs + COPs", receiver_celltype),
      score=correlation
      ) %>%
    filter(padj < 0.05) %>%
    select(source, target, ligand, receptor, score)
  return(return_df)
}

if (celltype_column=='Celltype_Class'){
    cell_cols<-c(`Excitatory neurons`="#006d2c", `OPCs + COPs`="#8c510a", `Oligodendrocytes`='#f781bf', `Microglia`='#e41a1c', `Inhibitory neurons`='#377eb8', `Astrocytes`='#FF7518',`Vascular Cells`= "#bebada")

} else if (celltype_column=='unified_subtype'){
    cell_cols <- c(`Exc L2-3 CBLN2 LINC02306`="#00441b", `Exc L3-4 RORB CUX2`="#006d2c", `Exc L3-5 RORB PLCH1`="#238b45",
  `Exc L4-5 RORB GABRG1`="#41ab5d", `Exc L4-5 RORB IL1RAPL2`="#74c476", `Exc L5 ET`="#a1d99b",
  `Exc L5/6 IT Car3`="#c7e9c0", `Exc L5/6 NP`="#2ca25f", `Exc L5-6 RORB LINC02196`="#66c2a4",
  `Exc L6 CT`="#99d8c9", `Exc L6 THEMIS NFIA`="#ccece6", `Exc L6b`="#5aae61",
  `Exc NRGN`="#1b7837", `Exc RELN CHD7`="#d9f0d3",
  `Inh ALCAM TRPM3`="#08306b", `Inh CUX2 MSR1`="#08519c", `Inh ENOX2 SPHKAP`="#2171b5",
  `Inh FBN2 EPB41L4A`="#4292c6", `Inh GPC5 RIT2`="#6baed6", `Inh L1 PAX6 CA4`="#9ecae1",
  `Inh L1-2 PAX6 SCGN`="#c6dbef", `Inh L1-6 LAMP5 CA13`="#9fb3c8", `Inh L3-5 SST MAFB`="#252b6c",
  `Inh L5-6 PVALB STON2`="#2c7fb8", `Inh L5-6 SST TH`="#3690c0", `Inh L6 SST NPY`="#74a9cf",
  `Inh LAMP5 NRG1 (Rosehip)`="#a6bddb", `Inh LAMP5 RELN`="#d0d1e6", `Inh PTPRK FAM19A1`="#1d3557",
  `Inh PVALB CA8 (Chandelier)`="#0570b0", `Inh PVALB HTR4`="#3690c0", `Inh PVALB SULF1`="#74a9cf",
  `Inh RYR3 TSHZ2`="#a6bddb", `Inh SGCD PDE3A`="#045a8d", `Inh SORCS1 TTN`="#2b8cbe",
  `Inh VIP ABI3BP`="#627d98", `Inh VIP CLSTN2`="#005082", `Inh VIP THSD7B`="#506680", `Inh VIP TSHZ2`="#41526a",
  `Ast DPP10`="#d94801", `Ast GRM3`="#ff7f2a", `Ast TPST1`="#ffae66",
  `Mic DUSP1`="#660000", `Mic MKI67`="#990000", `Mic P2RY12`="#cc0000", `Mic TPT1`="#e60000", `T cells`="#ff1a1a",
  `Oli OPALIN`="#fa91c9", `Oli RASGRF1`="#ff4da6",
  `OPC DOCK5`="#5c4033", `Opc GRIA4`="#8b5a2b", `Opc TPST1`="#d2b48c",
  `CAMs`="#3f007d", `Fib1`="#54278f", `Fib2`="#6a51a3", `Fib3`="#7b3294", `Per1`="#807dba", `Per2`="#9e9ac8",
  `aEndo`="#9e4db3", `aSMC`="#b084c1", `capEndo`="#cbc9e2", `vEndo`="#dadaeb", `vSMC`="#decbe4",
  `NA`="#000000")
  } else {
  unique_vals <- unique(cc_df[[celltype_column]])
  cell_cols <- setNames(rainbow(length(unique_vals)), unique_vals) #generate distinct colors
}

cc_df <- prep_cc_data(correlations)

#PNG plot
png(png_plot_path, width = 4000, height = 4000, res = 300)
cc_circos(cc_df, option = "A", cell_cols = cell_cols, cex = 2.5)
title(paste0(dataset_name, ": Cell–Cell Interaction Network"), line = -0.5)
dev.off()
circos.clear()

#PDF plot
pdf(pdf_plot_path, width = 10, height = 10, useDingbats = FALSE)
cc_circos(cc_df, option = "A", cell_cols = cell_cols, cex = 2.5)
title(paste0(dataset_name, ": Cell–Cell Interaction Network"), line = -0.5)
dev.off()
circos.clear()