# Load necessary packages
library(vegan)     # for adonis2
library(dplyr)     # for data manipulation
library(readr)     # for reading CSVs

# === Load the Aitchison distance matrix ===
aitchison_df <- read_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/aitchison_dist_matrix_prep_mobilome.csv")

# Convert to distance matrix
aitchison_mat <- as.matrix(aitchison_df[,-1])
rownames(aitchison_mat) <- aitchison_df[[1]]
distance_matrix <- as.dist(aitchison_mat)

# === Load metadata ===
metadata <- read_csv("C:/Users/asake/OneDrive/Desktop/Homework/FMT/metadata_for_prep.csv")

# Ensure row order matches distance matrix
metadata <- metadata %>% filter(ID %in% rownames(aitchison_mat))
metadata <- metadata[match(rownames(aitchison_mat), metadata$ID), ]

# Check alignment
stopifnot(all(metadata$ID == rownames(aitchison_mat)))

# === Run PERMANOVA ===
permanova_result <- adonis2(
  distance_matrix ~ fmt_prep,
  data = metadata,
  permutations = 999,
  strata = metadata$Patient # Controls for patient as a blocking factor
)

# === View results ===
print(permanova_result)


# library(pairwiseAdonis)
# 
# # Convert the distance object to a square matrix
# distance_matrix_mat <- as.matrix(distance_matrix)
# 
# # Subset to only rows/columns matching metadata
# distance_matrix_mat <- distance_matrix_mat[metadata$ID, metadata$ID]
# 
# # Run pairwise.adonis2 on the full matrix
# pairwise_results <- pairwise.adonis2(
#   distance_matrix_mat,
#   factors = metadata$fmt_prep,
#   strata = metadata$Patient,
#   perm = 999,
#   p.adjust.m = "BH"
# )
# 
# # View the result
# print(pairwise_results)
# 
# # View pairwise comparison results
# print(pairwise_results)