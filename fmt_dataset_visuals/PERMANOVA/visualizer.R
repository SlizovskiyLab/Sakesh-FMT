library(ggplot2)
library(dplyr)
library(tibble)

plot_permanova_r2 <- function(permanova_results_named_list, title = "PERMANOVA R² by Term") {
  # Collect R² values for each term from each result
  r2_df <- lapply(names(permanova_results_named_list), function(name) {
    result <- permanova_results_named_list[[name]]
    
    # Convert to data frame and remove "Residual" and "Total" rows
    result_df <- as.data.frame(result)
    result_df$Term <- rownames(result_df)
    result_df <- result_df %>% 
      filter(!(Term %in% c("Residual", "Total"))) %>%
      select(Term, R2)
    
    result_df$Comparison <- name
    return(result_df)
  }) %>% bind_rows()
  
  # Plot
  ggplot(r2_df, aes(x = R2, y = Term, fill = Comparison)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_minimal() +
    labs(
      title = title,
      x = expression(R^2),
      y = "PERMANOVA Term"
    ) +
    theme(
      axis.text.y = element_text(size = 10),
      legend.position = "bottom"
    )
}
