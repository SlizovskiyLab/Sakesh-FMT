run_pairwise_permanova_multi_groups <- function(distance_matrix_path, metadata_path, group_column, permutations = 999) {
  library(vegan)
  library(dplyr)
  library(readr)
  
  # Load metadata and distance matrix
  metadata <- read_csv(metadata_path)
  distance_df <- read_csv(distance_matrix_path)
  
  # fix rows
  distance_df <- distance_df[!duplicated(distance_df[[1]]), ]
  
  # Make row and column names unique and consistent
  row_ids <- make.unique(as.character(distance_df[[1]]))
  distance_mat_raw <- distance_df[, -1]
  colnames(distance_mat_raw) <- row_ids 
  
  # Convert to matrix
  distance_mat <- as.matrix(distance_mat_raw)
  rownames(distance_mat) <- row_ids
  
  # Make metadata IDs unique
  metadata$ID <- make.unique(as.character(metadata$ID))
  
  # Align metadata with matrix
  common_ids <- intersect(metadata$ID, rownames(distance_mat))
  metadata <- metadata %>% filter(ID %in% common_ids)
  distance_mat <- distance_mat[common_ids, common_ids]
  metadata <- metadata[match(common_ids, metadata$ID), ]
  stopifnot(all(metadata$ID == rownames(distance_mat)))
  
  # Prepare grouping and strata
  metadata$GroupFactor <- as.factor(metadata[[group_column]])
  metadata$Patient <- as.factor(metadata$Patient)
  
  # Generate pairwise group combinations
  group_levels <- levels(metadata$GroupFactor)
  pairwise_combos <- combn(group_levels, 2, simplify = FALSE)
  
  # Store results
  pairwise_results <- data.frame(
    Group1 = character(),
    Group2 = character(),
    R2 = numeric(),
    p.value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (pair in pairwise_combos) {
    g1 <- pair[1]
    g2 <- pair[2]
    
    sub_meta <- metadata %>% filter(GroupFactor %in% c(g1, g2))
    sub_mat <- distance_mat[sub_meta$ID, sub_meta$ID]
    
    result <- adonis2(
      as.dist(sub_mat) ~ GroupFactor,
      data = sub_meta,
      permutations = permutations,
      strata = sub_meta$Patient,
      method = "euclidean"
    )
    
    pairwise_results <- rbind(pairwise_results, data.frame(
      Group1 = g1,
      Group2 = g2,
      R2 = result$R2[1],
      p.value = result$`Pr(>F)`[1]
    ))
  }
  
  pairwise_results$padj <- p.adjust(pairwise_results$p.value, method = "fdr")
  cat("\n=== Pairwise PERMANOVA (3+ groups, stratified by Patient, FDR-corrected) ===\n")
  print(pairwise_results)
  
  return(pairwise_results)
}




extrac_mdrb_mob <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/aitchison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/metadata_mdrb.csv",
  group_column = "DNA_extraction_kit"
)

extrac_rcdi_mob <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/metadata_rcdi.csv",
  group_column = "DNA_extraction_kit"
)

# 
# route_mdrb_mob <- run_pairwise_permanova_multi_groups(
#   distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Route/aitchison_mdrb.csv",
#   metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Route/metadata_mdrb.csv",
#   group_column = "fmt_route"
# )
# route_rcdi_mob <- run_pairwise_permanova_multi_groups(
#   distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Route/aitchison_rcdi.csv",
#   metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Route/metadata_rcdi.csv",
#   group_column = "fmt_route"
# )

sequencer_rcdi_mob <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Sequencer/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Sequencer/metadata_rcdi.csv",
  group_column = "sequencer"
)

study_mdrb_mob <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/aitchison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/metadata_mdrb.csv",
  group_column = "study_data"
)

study_rcdi_mob <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/metadata_rcdi.csv",
  group_column = "study_data"
)


donrec_rcdi_mob <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/metadata_rcdi.csv",
  group_column = "donor_pre_post"
)
donrec_mdrb_mob <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/aitchison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/metadata_mdrb.csv",
  group_column = "donor_pre_post"
)

# donrec_mdrb_melanoma <- run_pairwise_permanova_multi_groups(
#   distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/aitchison_melanoma.csv",
#   metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/PrePostDonor/metadata_melanoma.csv",
#   group_column = "donor_pre_post"
# )






extrac_mdrb_res <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Extraction/aitchison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Extraction/metadata_mdrb.csv",
  group_column = "DNA_extraction_kit"
)

extrac_rcdi_res <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Extraction/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Extraction/metadata_rcdi.csv",
  group_column = "DNA_extraction_kit"
)


# route_mdrb_res <- run_pairwise_permanova_multi_groups(
#   distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Route/aitchison_mdrb.csv",
#   metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Route/metadata_mdrb.csv",
#   group_column = "fmt_route"
# )
# route_rcdi_res <- run_pairwise_permanova_multi_groups(
#   distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Route/aitchison_rcdi.csv",
#   metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Route/metadata_rcdi.csv",
#   group_column = "fmt_route"
# )

sequencer_rcdi_res <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Sequencer/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Sequencer/metadata_rcdi.csv",
  group_column = "sequencer"
)

study_mdrb_res <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/aitchison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/metadata_mdrb.csv",
  group_column = "study_data"
)

study_rcdi_res <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/metadata_rcdi.csv",
  group_column = "study_data"
)


donrec_rcdi_res <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/metadata_rcdi.csv",
  group_column = "donor_pre_post"
)
donrec_mdrb_res <- run_pairwise_permanova_multi_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/aitchison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/metadata_mdrb.csv",
  group_column = "donor_pre_post"
)

# donrec_mdrb_res <- run_pairwise_permanova_multi_groups(
#   distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/aitchison_melanoma.csv",
#   metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/PrePostDonor/metadata_melanoma.csv",
#   group_column = "donor_pre_post"
# )







