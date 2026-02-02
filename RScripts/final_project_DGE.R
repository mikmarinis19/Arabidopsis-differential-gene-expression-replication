# Load libraries
# Note that here I use edgeR, just like my paper does. 
library(edgeR)
library(tidyverse)
library(purrr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ggVennDiagram)

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
  
  # Return as a one-row tibble
  tibble(genotype = geno_name, low = low, high = high, total = total)
}

# Make tables 
table_COL <- summarize_DEGs(res_COL, "COL")
table_CH  <- summarize_DEGs(res_CH, "CH")
table_HY  <- summarize_DEGs(res_HY, "HY")

# Print tables
print(table_COL)
print(table_CH)
print(table_HY)

# Make venn diagram
venn_data <- list(
  "Col-0" = res_COL$gene_id,
  "chl1-5"  = res_CH$gene_id,
  "hy5-215"  = res_HY$gene_id
)

# Adjust and print venn diagram showing overlapping DEGs between the 3 genotypes.
venn_plot <- ggVennDiagram(venn_data, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "white") +  # makes fill white / invisible
  theme(legend.position = "none") +
  theme(text = element_text(size = 14))

print(venn_plot)

#Make heatmap

# Combine all DEGs into one table
all_res <- bind_rows(res_COL, res_CH, res_HY)

# Select top 50 by absolute logFC (unique genes)
top50_genes <- all_res %>%
  arrange(desc(abs(logFC))) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  slice_head(n = 50) %>%
  pull(gene_id)

# Calculate logCPM for all genes, for better comparison
dge <- DGEList(counts = counts_matrix)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE, prior.count = 1)

# Subset to top 50 DEGs
logCPM_top50 <- logCPM[top50_genes, ]

# Subset meta to only the samples in top50
# This way, we can add control vs treatment and genotype labels to the heatmap
annotation_col <- meta %>%
  dplyr::filter(paste0(Run, ".counts") %in% colnames(logCPM_top50)) %>%
  dplyr::mutate(Run = paste0(Run, ".counts")) %>%
  column_to_rownames("Run") %>%
  dplyr::select(genotype, treatment)

# Make the heatmap with labels 
pheatmap(
  logCPM_top50,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  scale = "row",
  annotation_col = annotation_col,
  main = "Top 50 DEGs Heatmap"
)

#Make PCA plot

# Select top 500 by absolute logFC (unique genes)
top500_genes <- all_res %>%
  arrange(desc(abs(logFC))) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  slice_head(n = 500) %>%
  pull(gene_id)

# Subset logCPM to top 500 DEGs
logCPM_top500 <- logCPM[top500_genes, ]

# PCA on top 500 DEGs

# Fix column names to match meta, get rid of "counts"
colnames(logCPM) <- str_remove(colnames(logCPM), "\\.counts$")

# Subset logCPM to top 500 DEGs that exist in counts
top500_genes <- top500_genes[top500_genes %in% rownames(logCPM)]
logCPM_top500 <- logCPM[top500_genes, ]

# Run PCA and variance calculation
pca <- prcomp(t(logCPM_top500))
percentVar <- (pca$sdev^2) / sum(pca$sdev^2)

# Build PCA dataframe
meta_sub <- meta %>% filter(Run %in% colnames(logCPM_top500)) %>%
  slice(match(colnames(logCPM_top500), Run))

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  genotype = meta_sub$genotype,
  treatment = meta_sub$treatment,
  Run = meta_sub$Run
)

# Plot PCA and print
p <-ggplot(pca_df, aes(x = PC1, y = PC2, color = genotype, shape = treatment, label = Run)) +
  geom_point(size = 3) +
  geom_label_repel(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1]*100, 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2]*100, 2), "% variance")) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

print(p)

#Print tables for the genes with the highest Log2FC for each genotype
top_ten_res_COL <- res_COL %>% arrange(desc(abs(logFC))) %>% slice(1:10) %>% select(gene_id, logFC)
top_ten_res_CH <- res_CH %>% arrange(desc(abs(logFC))) %>% slice(1:10) %>% select(gene_id, logFC)
top_ten_res_HY <- res_HY %>% arrange(desc(abs(logFC))) %>% slice(1:10) %>% select(gene_id, logFC)

print("Top DEGs for Col-0")
print(top_ten_res_COL)
print("Top DEGs for chl1-5")
print(top_ten_res_CH)
print("Top DEGs for hy5-215")
print(top_ten_res_HY)