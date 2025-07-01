run_permanova_two_groups <- function(distance_matrix_path, metadata_path, group_column, permutations = 999) {
  library(vegan)
  library(dplyr)
  library(readr)
  
  # Load metadata
  metadata <- read_csv(metadata_path)
  
  # Load Aitchison distance matrix
  distance_df <- read_csv(distance_matrix_path)
  distance_mat <- as.matrix(distance_df[,-1])
  rownames(distance_mat) <- distance_df[[1]]
  distance_dist <- as.dist(distance_mat)
  
  # Align metadata
  metadata <- metadata %>% filter(ID %in% rownames(distance_mat))
  metadata <- metadata[match(rownames(distance_mat), metadata$ID), ]
  stopifnot(all(metadata$ID == rownames(distance_mat)))
  
  # TEMP: Copy the grouping column into a fixed name
  metadata$GroupFactor <- metadata[[group_column]]
  
  # Run PERMANOVA using the fixed column name
  result <- adonis2(
    distance_dist ~ GroupFactor,
    data = metadata,
    permutations = permutations,
    strata = metadata$Patient,
    method = "euclidean"
  )
  
  cat("\n=== PERMANOVA: Two Groups ===\n")
  print(result)
  return(result)
}




prep_mdrb_mob <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Prep/aitichison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Prep/metadata_mdrb.csv",
  group_column = "fmt_prep"
)

prep_rcdi_mob <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Prep/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Prep/metadata_rcdi.csv",
  group_column = "fmt_prep"
)

extrac_melanoma_mob <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/aitchison_melanoma.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Extraction/metadata_melanoma.csv",
  group_column = "DNA_extraction_kit"
)

route_melanoma_mob <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Route/aitchison_melanoma.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Route/metadata_melanoma.csv",
  group_column = "fmt_route"
)

sequencer_melanoma_mob <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Sequencer/aitchison_melanoma.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Sequencer/metadata_melanoma.csv",
  group_column = "sequencer"
)

sequencer_mdrb_mob <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Sequencer/aitchison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Sequencer/metadata_mdrb.csv",
  group_column = "sequencer"
)

study_melanoma_mob <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/aitchison_melanoma.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Mobilome_PCA/Study/metadata_melanoma.csv",
  group_column = "study_data"
)


prep_mdrb_res <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Prep/aitchison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Prep/metadata_mdrb.csv",
  group_column = "fmt_prep"
)

prep_rcdi_res <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Prep/aitchison_rcdi.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Prep/metadata_rcdi.csv",
  group_column = "fmt_prep"
)

extrac_melanoma_res <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Extraction/aitchison_melanoma.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Extraction/metadata_melanoma.csv",
  group_column = "DNA_extraction_kit"
)

route_melanoma_res <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Route/aitchison_melanoma.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Route/metadata_melanoma.csv",
  group_column = "fmt_route"
)

sequencer_melanoma_res <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Sequencer/aitchison_melanoma.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Sequencer/metadata_melanoma.csv",
  group_column = "sequencer"
)

sequencer_mdrb_res <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Sequencer/aitchison_mdrb.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Sequencer/metadata_mdrb.csv",
  group_column = "sequencer"
)

study_melanoma_res <- run_permanova_two_groups(
  distance_matrix_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/aitchison_melanoma.csv",
  metadata_path = "C:/Users/asake/OneDrive/Desktop/Homework/FMT/Resistome_PCA/Study/metadata_melanoma.csv",
  group_column = "study_data"
)

