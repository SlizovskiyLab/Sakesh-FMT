plot_r2_permanova_individual <- function(permanova_results_named_list, output_dir = NULL, prefix = "permanova_r2") {
  library(ggplot2)
  library(dplyr)
  
  for (name in names(permanova_results_named_list)) {
    result <- permanova_results_named_list[[name]]
    
    # Convert result to data frame and extract terms
    result_df <- as.data.frame(result)
    result_df$Term <- rownames(result_df)
    result_df <- result_df %>%
      filter(!(Term %in% c("Residual", "Total"))) %>%
      select(Term, R2)
    
    # Create plot
    p <- ggplot(result_df, aes(x = reorder(Term, -R2), y = R2)) +
      geom_bar(stat = "identity", fill = "darkorange") +
      theme_minimal() +
      labs(
        title = paste("PERMANOVA R² -", name),
        x = "Term",
        y = expression(R^2)
      ) +
      ylim(0, 1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p)  # Show plot
    
    # Optionally save plot
    if (!is.null(output_dir)) {
      ggsave(filename = file.path(output_dir, paste0(prefix, "_", name, ".png")), plot = p, width = 6, height = 4)
    }
  }
}


all_permanova_results <- list(
  prep_mdrb_mob = prep_mdrb_mob,
  prep_rcdi_mob = prep_rcdi_mob,
  extrac_melanoma_mob = extrac_melanoma_mob,
  route_melanoma_mob = route_melanoma_mob,
  sequencer_melanoma_mob = sequencer_melanoma_mob,
  sequencer_mdrb_mob = sequencer_mdrb_mob,
  study_melanoma_mob = study_melanoma_mob,
  prep_mdrb_res = prep_mdrb_res,
  prep_rcdi_res = prep_rcdi_res,
  extrac_melanoma_res = extrac_melanoma_res,
  route_melanoma_res = route_melanoma_res,
  sequencer_melanoma_res = sequencer_melanoma_res,
  sequencer_mdrb_res = sequencer_mdrb_res,
  study_melanoma_res = study_melanoma_res
)

plot_r2_permanova_individual(all_permanova_results, output_dir = NULL)
