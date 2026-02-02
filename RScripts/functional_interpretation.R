# Load libraries
# Note that here I use edgeR, just like my paper does. 
library(dplyr)
library(edgeR)
library(tidyverse)
library(purrr)
library(tibble)

# Load and clean metadata
# Since I am only interested in root data, I filtered for root
meta <- read.csv("Downloads/SraRunTable_Final_Project.csv") %>%
  mutate(
    genotype = recode(genotype,
                      "chl1.5" = "CH",
                      "hy5-215" = "HY",
                      "Col-0" = "COL"),
    treatment = recode(treatment,
                       "21C" = "control",
                       "28C" = "exposed")
  ) %>%
  dplyr::filter(tissue == "Root") %>%
  dplyr::select(Run, genotype, treatment)

# Load count files
directory <- "Downloads/Final_Project_Counts"
sampleFiles <- paste0(meta$Run, ".counts")

# Check that files match metadata
# This part successfully returns "True"
all(str_remove(sampleFiles, ".counts") == meta$Run) 
print(all(str_remove(sampleFiles, ".counts") == meta$Run))

# Merge all counts by gene_id into a single table
# When I run edgeR, it will grab sub data from my counts_matrix object
counts_matrix <- purrr::map(file.path(directory, sampleFiles),
                            ~ read.delim(.x, header = FALSE, col.names = c("gene_id", basename(.x)))) %>%
  purrr::reduce(full_join, by = "gene_id") %>%   
  column_to_rownames("gene_id") %>%
  as.matrix()

# This function takes in the genotype. 
# It uses an LFC threshold with an absolute value of at least 1
# It uses an FDR threshold of 0.05
# These two thresholds are consistent with my paper
# The use of edgeR is also from my paper 

run_edgeR <- function(geno, lfc_threshold = 1, fdr_threshold = 0.05) {
  
  #This only looks at the subsections of meta and my counts_matrix that are relevant to my genotype
  meta_sub <- meta %>% filter(genotype == geno)
  counts_sub <- counts_matrix[, colnames(counts_matrix) %in% paste0(meta_sub$Run, ".counts")]
  
  # Factor the data
  sample_groups <- factor(meta_sub$treatment)
  dge <- DGEList(counts = counts_sub, group = sample_groups)
  
  # Filter out lowly expressed genes
  keep <- rowSums(cpm(dge) > 1) >= 2
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Fit the data to a model
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ sample_groups)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  
  # Here, I make the contrast vector. 
  # The intercept (0) is my control group, which are roots exposed to a control temperature. 
  # 1 represents my treatment group, which are roots exposed to a heightened temperature. 
  contrast <- c(0, 1)  
  
  # I now apply the contrast vector to the model
  res <- glmQLFTest(fit, contrast = contrast)
  res_df <- topTags(res, n = nrow(dge), adjust.method = "BH", sort.by = "none")$table
  res_df$gene_id <- rownames(res_df)
  
  # Apply thresholds from the paper: FDR < 0.05 & |log2FC| > 1
  res_df <- res_df %>%
    filter(FDR < fdr_threshold & abs(logFC) > lfc_threshold)
  
  return(res_df)
}


# Run my edgeR function on each genotype
res_COL <- run_edgeR("COL")
res_CH  <- run_edgeR("CH")
res_HY  <- run_edgeR("HY")

# Make a table of DEGs for each genotype
summarize_DEGs <- function(res_df, geno_name) {
  low  <- sum(res_df$logFC < -1)
  high <- sum(res_df$logFC > 1)
  total <- low + high
  
  # Return as a tibble
  tibble(genotype = geno_name, low = low, high = high, total = total)
}

# Make tables 
table_COL <- summarize_DEGs(res_COL, "COL")
table_CH  <- summarize_DEGs(res_CH, "CH")
table_HY  <- summarize_DEGs(res_HY, "HY")

# Functional interpretation

#Load libraries
library(biomaRt)
library(clusterProfiler)
library(org.At.tair.db)

#Over-representation analysis, treatment groups vs. control groups
# Background = all genes in counts matrix
background <- rownames(counts_matrix)

#Get arabidopsis functional annotation
mart <- useMart(
  biomart = "plants_mart",
  dataset = "athaliana_eg_gene",
  host = "https://plants.ensembl.org"
)

run_GO_for_foreground <- function(foreground, mart, background) {
  
  #Map foreground DEGs to TAIR 
  mapping_foreground <- getBM(
    attributes = c("ensembl_gene_id", "tair_locus"),
    filters = "ensembl_gene_id",
    values = foreground,
    mart = mart
  )
  foreground_tair <- mapping_foreground$tair_locus
  
  # Map all genes, which is my background, to TAIR 
  mapping_background <- getBM(
    attributes = c("ensembl_gene_id", "tair_locus"),
    values = background,
    filters = "ensembl_gene_id",
    mart = mart
  )
  background_tair <- mapping_background$tair_locus
  
  #GO enrichment
  ego <- enrichGO(
    gene          = as.character(foreground_tair),
    universe      = as.character(background_tair),
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.1,
    readable      = TRUE
  )
  
  return(ego)
}

# Genes in Col-0 but NOT in chl1-5 (NRT1.1-dependent)
NRT1.1_dependent <- setdiff(res_COL$gene_id, res_CH$gene_id)

# Genes in Col-0 but NOT in hy5-215 (HY5-dependent)
HY5_dependent <- setdiff(res_COL$gene_id, res_HY$gene_id)

# Genes in Col-0 AND chl1-5, but NOT in hy5-215 (HY5-NRT1.1-dependent)
HY5_NRT1.1_dependent <- setdiff(intersect(res_COL$gene_id, res_CH$gene_id), res_HY$gene_id)


GO_NRT1.1 <- run_GO_for_foreground(NRT1.1_dependent, mart, background)
GO_HY5    <- run_GO_for_foreground(HY5_dependent, mart, background)
GO_HY5_NRT1.1 <- run_GO_for_foreground(HY5_NRT1.1_dependent, mart, background)

print(GO_NRT1.1)
print(GO_HY5)
print(GO_HY5_NRT1.1)

plot1 <- dotplot(GO_HY5, showCategory = 15) +
  ggtitle("HY5 dependent GO enrichment")

plot2 <- dotplot(GO_HY5_NRT1.1, showCategory = 15) +
  ggtitle("HY5-NRT1.1 dependent GO enrichment")

print(plot1)
print(plot2)
