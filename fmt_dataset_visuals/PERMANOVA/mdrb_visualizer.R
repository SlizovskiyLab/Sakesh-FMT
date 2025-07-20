# Helper function to extract R2 and p-value
extract_adonis_result <- function(adonis_result, factor_name, data_type) {
  data.frame(
    Factor = factor_name,
    DataType = data_type,  # "Mobilome" or "Resistome"
    R2 = adonis_result$R2[1],
    P = adonis_result$`Pr(>F)`[1]
  )
}

# Collect all results
results <- bind_rows(
  extract_adonis_result(prep_mdrb_mob, "FMT Prep", "Mobilome"),
  extract_adonis_result(prep_mdrb_res, "FMT Prep", "Resistome"),
  extract_adonis_result(extrac_mdrb_mob, "DNA Extraction", "Mobilome"),
  extract_adonis_result(extrac_mdrb_res, "DNA Extraction", "Resistome"),
  extract_adonis_result(route_mdrb_mob, "FMT Route", "Mobilome"),
  extract_adonis_result(route_mdrb_res, "FMT Route", "Resistome"),
  extract_adonis_result(sequencer_mdrb_mob, "Sequencer", "Mobilome"),
  extract_adonis_result(sequencer_mdrb_res, "Sequencer", "Resistome"),
  extract_adonis_result(study_mdrb_mob, "Study", "Mobilome"),
  extract_adonis_result(study_mdrb_res, "Study", "Resistome"),
  extract_adonis_result(donrec_mdrb_mob, "Donor/Recipient", "Mobilome"),
  extract_adonis_result(donrec_mdrb_res, "Donor/Recipient", "Resistome")
)



library(ggplot2)
library(dplyr)
library(patchwork)

# Assign row for tile matrix
results <- results %>%
  mutate(TileRow = DataType)

# Order Factors by max R²
factor_order <- results %>%
  group_by(Factor) %>%
  summarise(max_r2 = max(R2)) %>%
  arrange(desc(max_r2)) %>%
  pull(Factor)

results$Factor <- factor(results$Factor, levels = factor_order)

# Tile rows will be explicitly ordered
results$TileRow <- factor(results$TileRow, levels = c("Mobilome", "Resistome"))


bar_plot <- ggplot(results, aes(x = Factor, y = R2, fill = DataType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values = c("Mobilome" = "#FBB04B", "Resistome" = "#80B1D3")) +
  labs(y = "Adonis R²", x = NULL, fill = "MDRB") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),  # hide to avoid duplication
    axis.title.x = element_blank(),
    plot.margin = margin(b = 0),   # reduce spacing
    legend.position = "top"
  )


tile_grid <- expand.grid(
  Factor = levels(results$Factor),
  TileRow = c("Mobilome", "Resistome"),
  stringsAsFactors = FALSE
)

tile_data <- left_join(tile_grid, results, by = c("Factor", "TileRow")) %>%
  mutate(
    Significant = ifelse(!is.na(P) & P < 0.05, "p < 0.05", "Not Significant")
  )

# Plot with red legend
tile_plot <- ggplot(tile_data, aes(x = Factor, y = TileRow)) +
  geom_tile(aes(fill = Significant), width = 0.6, height = 0.6, color = "gray80") +
  scale_fill_manual(
    values = c("p < 0.05" = "red", "Not Significant" = "white"),
    breaks = c("p < 0.05"),
    name = NULL
  ) +
  scale_y_discrete(position = "right") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    plot.margin = margin(t = 0),
    axis.text.y = element_text(size = 12),
    legend.position = "top"
  )



final_plot <- bar_plot / tile_plot + plot_layout(heights = c(3, 1))
print(final_plot)
