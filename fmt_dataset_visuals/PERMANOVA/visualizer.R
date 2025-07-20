plot_pairwise_permanova_results <- function(pairwise_results) {
  library(ggplot2)
  library(dplyr)
  
  # Prepare data
  plot_data <- pairwise_results %>%
    mutate(
      Variable = paste(Group1, "vs", Group2),
      Significant = padj < 0.05
    )
  
  # Order by R2
  plot_data$Variable <- factor(
    plot_data$Variable,
    levels = plot_data$Variable[order(plot_data$R2, decreasing = TRUE)]
  )
  
  # Plot
  ggplot(plot_data, aes(x = Variable, y = R2)) +
    geom_bar(stat = "identity", fill = "goldenrod") +
    geom_tile(
      data = plot_data %>% filter(Significant),
      aes(y = -0.005, height = 0.005, width = 0.9),
      fill = "red"
    ) +
    coord_flip() +
    theme_minimal() +
    labs(
      y = "Adonis R²",
      x = "",
      title = "Pairwise PERMANOVA Results"
    ) +
    theme(axis.text.y = element_text(size = 10))
}

plot_pairwise_permanova_results(donrec_rcdi_mob)



