#load libraries
library(tidyverse)

# Load and clean metadata
# Since I am only interested in root data, I filtered for root
meta <- read.csv("Downloads/SraRunTable_Final_Project.csv") %>%
  mutate(
    genotype = str_replace(genotype, "chl.*", "CH"),
    genotype = str_replace(genotype, "hy5.*", "HY"),
    genotype = str_replace(genotype, "Co.*", "COL"),
    treatment = case_when(
      treatment == "21C" ~ "control",
      treatment == "28C" ~ "exposed"
    )
  ) %>%
  filter(tissue == "Root") %>%
  select(Run, genotype, treatment)


# Load count files
directory <- "Downloads/Final_Project_Counts"
sampleFiles <- paste0(meta$Run, ".counts")

#Make sampletable for genotype and treatment
sampleTable <- data.frame(
  sampleName = meta$Run,
  fileName = sampleFiles,
  genotype = meta$genotype,
  treatment = meta$treatment
)

# write out sampleTable
sampleTable

# Check that files match metadata
# This part successfully returns "True"
all(str_remove(sampleFiles, ".counts") == meta$Run) 
print(all(str_remove(sampleFiles, ".counts") == meta$Run))

# Read SN.txt
sn_raw <- read.delim("Downloads/SN.txt", header = TRUE, check.names = FALSE)

# Get the row for raw total sequences
raw_row <- sn_raw %>% filter(grepl("raw total sequences", rownames(sn_raw) %||% .[[1]]))

# Make it so the sample numbers are rows
# I'm not sure why, but probably due to a weird formatting problem, I had to clean it
raw_seq <- data.frame(
  sampleID = str_remove(colnames(raw_row), "^\\.\\./"),        
  raw_total_sequences = as.numeric(raw_row[1, ])
)

# Clean it in case of formatting issue
raw_seq$sampleID <- str_remove(raw_seq$sampleID, "^\\.\\./")

#Check to make sure it prints properly, it does (18 samples)
print(raw_seq$sampleID)

# Join with metadata
sn <- left_join(sampleTable, raw_seq, by = c("sampleName" = "sampleID"))

# Now plot
plot <- ggplot(sn, aes(x = treatment, y = raw_total_sequences)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~ genotype)

# Print plot
print(plot)