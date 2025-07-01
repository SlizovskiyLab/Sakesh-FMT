run_pairwise_permanova_multi_group <- function(distance_matrix_path, metadata_path, group_column, permutations = 999) {
  library(vegan)
  library(dplyr)
  library(readr)
  
  # Read input
  metadata <- read_csv(metadata_path)
  distance_df <- read_csv(distance_matrix_path)
  
  # Format distance matrix
  distance_mat <- as.matrix(distance_df[,-1])
  rownames(distance_mat) <- distance_df[[1]]
  
  # Align metadata
  metadata <- metadata %>% filter(ID %in% rownames(distance_mat))
  metadata <- metadata[match(rownames(distance_mat), metadata$ID), ]
  metadata[[group_column]] <- as.factor(metadata[[group_column]])
  
  group_levels <- levels(metadata[[group_column]])
  combos <- combn(group_levels, 2, simplify = FALSE)
  
  results <- data.frame(
    Group1 = character(),
    Group2 = character(),
    R2 = numeric(),
    p.value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (pair in combos) {
    g1 <- pair[1]
    g2 <- pair[2]
    
    sub_meta <- metadata %>%
      filter(.data[[group_column]] %in% c(g1, g2))
    
    if (any(table(sub_meta[[group_column]]) < 2)) {
      warning(paste("Skipping", g1, "vs", g2, "— not enough samples in one or both groups"))
      next
    }
    
    sub_mat <- distance_mat[sub_meta$ID, sub_meta$ID]
    
    res <- adonis2(
      as.dist(sub_mat) ~ .data[[group_column]],
      data = sub_meta,
      permutations = permutations
    )
    
    results <- rbind(results, data.frame(
      Group1 = g1,
      Group2 = g2,
      R2 = res$R2[1],
      p.value = res$`Pr(>F)`[1]
    ))
  }
  
  results$padj <- p.adjust(results$p.value, method = "fdr")
  return(results)
}
