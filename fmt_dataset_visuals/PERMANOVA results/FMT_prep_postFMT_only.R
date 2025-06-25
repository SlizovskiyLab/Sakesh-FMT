# Load necessary packages
library(vegan)     # for adonis2
library(dplyr)     # for data manipulation
library(readr)     # for reading CSVs

# === Load metadata (shared for both datasets) ===
metadata <- read_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/metadata_for_PostFMTonly_prep.csv")

# === MOBILOME PERMANOVA ===

# Load Aitchison distance matrix for mobilome
mobilome_df <- read_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/aitchison_dist_matrix_prep_PostFMTonly_mobilome.csv")
mobilome_mat <- as.matrix(mobilome_df[,-1])
rownames(mobilome_mat) <- mobilome_df[[1]]
mobilome_dist <- as.dist(mobilome_mat)

# Align metadata with mobilome
mobilome_meta <- metadata %>% filter(ID %in% rownames(mobilome_mat))
mobilome_meta <- mobilome_meta[match(rownames(mobilome_mat), mobilome_meta$ID), ]
stopifnot(all(mobilome_meta$ID == rownames(mobilome_mat)))

# Run PERMANOVA for mobilome
mobilome_permanova <- adonis2(
  mobilome_dist ~ fmt_prep,
  data = mobilome_meta,
  permutations = 999,
  strata = mobilome_meta$Patient
)

# Print mobilome result
cat("\n=== PERMANOVA: Mobilome ===\n")
print(mobilome_permanova)


# === RESISTOME PERMANOVA ===

# Load Aitchison distance matrix for resistome
resistome_df <- read_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/aitchison_dist_matrix_prep_resistome_PostFMTonly.csv")
resistome_mat <- as.matrix(resistome_df[,-1])
rownames(resistome_mat) <- resistome_df[[1]]
resistome_dist <- as.dist(resistome_mat)

# Align metadata with resistome
resistome_meta <- metadata %>% filter(ID %in% rownames(resistome_mat))
resistome_meta <- resistome_meta[match(rownames(resistome_mat), resistome_meta$ID), ]
stopifnot(all(resistome_meta$ID == rownames(resistome_mat)))

# Run PERMANOVA for resistome
resistome_permanova <- adonis2(
  resistome_dist ~ fmt_prep,
  data = resistome_meta,
  permutations = 999,
  strata = resistome_meta$Patient
)

# Print resistome result
cat("\n=== PERMANOVA: Resistome ===\n")
print(resistome_permanova)


# === FDR-adjusted p-values (Benjamini-Hochberg) ===

# Extract raw p-values
raw_pvals <- c(
  mobilome = mobilome_permanova$`Pr(>F)`[1],
  resistome = resistome_permanova$`Pr(>F)`[1]
)

# Apply BH correction
fdr_adjusted <- p.adjust(raw_pvals, method = "BH")

# Report FDR-corrected p-values
cat("\n=== FDR-adjusted p-values ===\n")
fdr_result <- data.frame(
  Test = names(raw_pvals),
  Raw_P = raw_pvals,
  FDR_Adjusted_P = fdr_adjusted
)
print(fdr_result)
