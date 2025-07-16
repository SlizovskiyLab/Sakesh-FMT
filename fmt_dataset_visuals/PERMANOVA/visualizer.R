library(ggplot2)
library(dplyr)
library(tidyr)

plot_pairwise_permanova_heatmap <- function(pairwise_df, title = NULL, fdr_threshold = 0.05) {
  # Check required columns
  stopifnot(all(c("Group1", "Group2", "R2", "padj") %in% names(pairwise_df)))
  
  # Filter to significant results only
  sig_df <- pairwise_df %>%
    filter(padj < fdr_threshold) %>%
    mutate(Signif = "*")
  
  # If no significant results, return message
  if (nrow(sig_df) == 0) {
    message("No significant pairwise comparisons (padj < ", fdr_threshold, ")")
    return(NULL)
  }
  
  # Create symmetric pairs
  sig_df_sym <- sig_df %>%
    bind_rows(sig_df %>%
                rename(Group1 = Group2, Group2 = Group1)) %>%
    mutate(Signif = "*")
  
  # Ensure full factor levels for axis consistency
  all_groups <- unique(c(sig_df$Group1, sig_df$Group2))
  
  # Plot
  ggplot(sig_df_sym, aes(x = Group1, y = Group2, fill = R2)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("R²=%.2f\n%s", R2, Signif)), size = 3.5) +
    scale_fill_gradient(low = "white", high = "firebrick", na.value = "grey90") +
    scale_x_discrete(limits = all_groups) +
    scale_y_discrete(limits = all_groups) +
    coord_fixed() +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = "PERMANOVA R²"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}