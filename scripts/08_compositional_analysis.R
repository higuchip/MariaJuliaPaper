# ==============================================================================
# COMPOSITIONAL ANALYSIS: PERMANOVA, BETA PARTITIONING, AND PCoA
# ==============================================================================
#
# Description: Analyzes compositional differences between regenerating
#              communities in abandoned pastures and reference forests using
#              PERMANOVA, beta diversity partitioning (turnover vs nestedness),
#              and PCoA visualization. Demonstrates that shrubland communities
#              are distinct assemblages, not impoverished forest subsets.
#
# Project:     Vegetation dynamics in abandoned Atlantic Forest highland pastures
#              São Joaquim National Park, Brazil (2014-2025)
#
# Author:      Pedro Higuchi
# Contact:     higuchip@gmail.com
# Repository:  https://github.com/higuchip/MariaJuliaPaper
#
# Created:     2025
# Last update: 2025-12-16
#
# Citation:    Cruz MJC et al. (in review) Forest recovery or shrubland assembly?
#              Vegetation dynamics in abandoned Atlantic Forest highland pastures.
#              Biotropica.
#
# License:     MIT License (see LICENSE file)
#
# Input:       ref_dados_floresta from 03_floristic_similarity.R
#              df_processed from 01_data_processing.R
#
# Output:      ../artigo/biotropica/review/Figure_S5_PCoA_composition.tiff
#              outputs/tables/Compositional_Analysis_Results.csv
#
# Note:        Requires objects from scripts 01 and 03
#              R project located in MariaJulia/dissertacao
#
# ==============================================================================

# --- Setup --------------------------------------------------------------------

set.seed(123)

# Required packages
suppressPackageStartupMessages({
  library(vegan)
  library(betapart)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Output directories
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)

# Output path for publication figures
output_path <- "../artigo/biotropica/review/"

cat("================================================================================\n")
cat("COMPOSITIONAL ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

cat("1. DATA PREPARATION\n")
cat("-------------------\n\n")

# Check dependencies
if(!exists("ref_dados_floresta")) {
  stop("Object 'ref_dados_floresta' not found. Run 03_floristic_similarity.R first.")
}
if(!exists("df_processed")) {
  stop("Object 'df_processed' not found. Run 01_data_processing.R first.")
}

# --- A. Reference Forest data ---
cat("Processing reference forest data...\n")

comp_forest <- ref_dados_floresta %>%
  mutate(
    Plot_ID = paste0("Forest_A", Area, "_", Parcela),
    Ecosystem = "Reference forest"
  ) %>%
  select(Plot_ID, Ecosystem, Species = Especie) %>%
  mutate(Presence = 1) %>%
  distinct()

n_forest_plots <- n_distinct(comp_forest$Plot_ID)
n_forest_species <- n_distinct(comp_forest$Species)
cat(sprintf("  Reference forest: %d plots, %d species\n", n_forest_plots, n_forest_species))

# --- B. Regeneration data (most recent year) ---
cat("Processing regeneration data...\n")

comp_height_cols <- grep("^h\\d{4}$", names(df_processed), value = TRUE)
comp_recent_year <- max(as.numeric(gsub("h", "", comp_height_cols)))
comp_year_col <- paste0("h", comp_recent_year)

cat(sprintf("  Using regeneration census: %d\n", comp_recent_year))

# Abandoned pastures (Areas 1, 2, 4)
comp_pastures <- df_processed %>%
  filter(area %in% c(1, 2, 4)) %>%
  filter(!is.na(!!sym(comp_year_col)) & !!sym(comp_year_col) > 0) %>%
  mutate(
    Plot_ID = paste0("Pasture_A", area, "_", parc.1),
    Ecosystem = "Abandoned pastures"
  ) %>%
  select(Plot_ID, Ecosystem, Species = especies) %>%
  mutate(Presence = 1) %>%
  distinct()

n_pasture_plots <- n_distinct(comp_pastures$Plot_ID)
n_pasture_species <- n_distinct(comp_pastures$Species)
cat(sprintf("  Abandoned pastures: %d plots, %d species\n", n_pasture_plots, n_pasture_species))

# --- C. Combine datasets ---
cat("\nCombining datasets...\n")

comp_combined <- bind_rows(comp_forest, comp_pastures)

cat(sprintf("  Total: %d plots, %d species\n", 
            n_distinct(comp_combined$Plot_ID),
            n_distinct(comp_combined$Species)))

# --- D. Create presence/absence matrix ---
cat("Creating presence/absence matrix...\n")

comp_matrix_wide <- comp_combined %>%
  pivot_wider(
    id_cols = c(Plot_ID, Ecosystem),
    names_from = Species,
    values_from = Presence,
    values_fill = 0
  )

comp_ecosystem <- factor(comp_matrix_wide$Ecosystem)
comp_plot_ids <- comp_matrix_wide$Plot_ID

comp_pa_matrix <- comp_matrix_wide %>%
  select(-Plot_ID, -Ecosystem) %>%
  as.data.frame()

rownames(comp_pa_matrix) <- comp_plot_ids

cat(sprintf("  Matrix dimensions: %d plots × %d species\n", 
            nrow(comp_pa_matrix), ncol(comp_pa_matrix)))
cat("\n  Plots per ecosystem:\n")
print(table(comp_ecosystem))

# ==============================================================================
# 2. SPECIES POOL ANALYSIS
# ==============================================================================

cat("\n================================================================================\n")
cat("2. SPECIES POOL ANALYSIS\n")
cat("================================================================================\n\n")

# Species lists
spp_forest <- comp_combined %>% 
  filter(Ecosystem == "Reference forest") %>% 
  pull(Species) %>% unique()

spp_pastures <- comp_combined %>% 
  filter(Ecosystem == "Abandoned pastures") %>% 
  pull(Species) %>% unique()

# Set operations
spp_shared <- intersect(spp_forest, spp_pastures)
spp_only_forest <- setdiff(spp_forest, spp_pastures)
spp_only_pastures <- setdiff(spp_pastures, spp_forest)

# Summary
cat("Species pool comparison:\n")
cat(sprintf("  Reference forest: %d species\n", length(spp_forest)))
cat(sprintf("  Abandoned pastures: %d species\n", length(spp_pastures)))
cat(sprintf("  Shared: %d species (%.1f%% of pasture flora)\n", 
            length(spp_shared), 
            length(spp_shared)/length(spp_pastures)*100))
cat(sprintf("  Exclusive to forest: %d species\n", length(spp_only_forest)))
cat(sprintf("  Exclusive to pastures: %d species (%.1f%% of pasture flora)\n", 
            length(spp_only_pastures),
            length(spp_only_pastures)/length(spp_pastures)*100))

cat("\nSpecies EXCLUSIVE to abandoned pastures:\n")
print(sort(spp_only_pastures))

# ==============================================================================
# 3. PERMANOVA ANALYSIS
# ==============================================================================

cat("\n================================================================================\n")
cat("3. PERMANOVA ANALYSIS\n")
cat("================================================================================\n\n")

# Calculate Jaccard distance matrix
comp_dist <- vegdist(comp_pa_matrix, method = "jaccard", binary = TRUE)

# Run PERMANOVA
cat("Running PERMANOVA (999 permutations)...\n\n")

comp_permanova <- adonis2(comp_dist ~ comp_ecosystem, permutations = 999)

print(comp_permanova)

# Extract results
permanova_F <- comp_permanova$F[1]
permanova_R2 <- comp_permanova$R2[1]
permanova_p <- comp_permanova$`Pr(>F)`[1]

cat(sprintf("\nPERMANOVA Results:\n"))
cat(sprintf("  F = %.3f\n", permanova_F))
cat(sprintf("  R² = %.3f (%.1f%% of variance explained)\n", permanova_R2, permanova_R2*100))
cat(sprintf("  p = %.3f\n", permanova_p))

# ==============================================================================
# 4. BETA DIVERSITY PARTITIONING
# ==============================================================================

cat("\n================================================================================\n")
cat("4. BETA DIVERSITY PARTITIONING\n")
cat("================================================================================\n\n")

# Aggregate to ecosystem-level pools
forest_pool <- colSums(comp_pa_matrix[comp_ecosystem == "Reference forest", ]) > 0
pastures_pool <- colSums(comp_pa_matrix[comp_ecosystem == "Abandoned pastures", ]) > 0

ecosystem_pool_matrix <- rbind(
  "Reference forest" = as.numeric(forest_pool),
  "Abandoned pastures" = as.numeric(pastures_pool)
)
colnames(ecosystem_pool_matrix) <- colnames(comp_pa_matrix)

cat("Ecosystem-level species pools:\n")
cat(sprintf("  Reference forest: %d species\n", sum(forest_pool)))
cat(sprintf("  Abandoned pastures: %d species\n", sum(pastures_pool)))

# Beta diversity partitioning using betapart
cat("\nCalculating beta diversity components (Jaccard family)...\n\n")

beta_pair_result <- beta.pair(ecosystem_pool_matrix, index.family = "jaccard")

# Extract components
beta_total <- as.numeric(beta_pair_result$beta.jac)
beta_turnover <- as.numeric(beta_pair_result$beta.jtu)
beta_nestedness <- as.numeric(beta_pair_result$beta.jne)

# Calculate percentages
pct_turnover <- (beta_turnover / beta_total) * 100
pct_nestedness <- (beta_nestedness / beta_total) * 100

cat("Beta diversity partitioning results:\n")
cat(sprintf("  Total dissimilarity (β_jac): %.3f\n", beta_total))
cat(sprintf("  Turnover component (β_jtu): %.3f (%.1f%%)\n", beta_turnover, pct_turnover))
cat(sprintf("  Nestedness component (β_jne): %.3f (%.1f%%)\n", beta_nestedness, pct_nestedness))


# ==============================================================================
# 5. BETADISPER ANALYSIS
# ==============================================================================

cat("\n================================================================================\n")
cat("5. BETADISPER ANALYSIS\n")
cat("================================================================================\n\n")

# Betadisper for PCoA coordinates
comp_betadisper <- betadisper(comp_dist, comp_ecosystem, type = "centroid")

# Extract eigenvalues for variance explained
comp_eig <- comp_betadisper$eig
comp_var_explained <- round(100 * comp_eig / sum(comp_eig), 1)

cat("PCoA variance explained:\n")
cat(sprintf("  Axis 1: %.1f%%\n", comp_var_explained[1]))
cat(sprintf("  Axis 2: %.1f%%\n", comp_var_explained[2]))

# Test for dispersion differences
cat("\nBetadisper test (homogeneity of dispersions):\n")
comp_betadisper_test <- anova(comp_betadisper)
print(comp_betadisper_test)

# Mean distances to centroid
cat("\nMean distance to centroid:\n")
comp_distances_df <- data.frame(
  Ecosystem = comp_ecosystem,
  Distance = comp_betadisper$distances
) %>%
  group_by(Ecosystem) %>%
  summarise(
    Mean = mean(Distance),
    SD = sd(Distance),
    n = n(),
    .groups = 'drop'
  )
print(comp_distances_df)

# ==============================================================================
# 6. PCoA FIGURE
# ==============================================================================

cat("\n================================================================================\n")
cat("6. GENERATING PCoA FIGURE\n")
cat("================================================================================\n\n")

# Extract PCoA scores
pcoa_scores <- data.frame(
  Plot_ID = rownames(comp_pa_matrix),
  PCoA1 = comp_betadisper$vectors[, 1],
  PCoA2 = comp_betadisper$vectors[, 2],
  Ecosystem = comp_ecosystem
)

# Extract centroids
pcoa_centroids <- data.frame(
  Ecosystem = rownames(comp_betadisper$centroids),
  PCoA1 = comp_betadisper$centroids[, 1],
  PCoA2 = comp_betadisper$centroids[, 2]
)

# Define colors
comp_colors <- c(
  "Abandoned pastures" = "#D55E00",
  "Reference forest" = "#0072B2"
)

# Create figure
fig_pcoa <- ggplot() +
  # Confidence ellipses (95%)
  stat_ellipse(data = pcoa_scores, 
               aes(x = PCoA1, y = PCoA2, color = Ecosystem, fill = Ecosystem),
               level = 0.95, 
               geom = "polygon",
               alpha = 0.15,
               linewidth = 1.2) +
  # Points
  geom_point(data = pcoa_scores, 
             aes(x = PCoA1, y = PCoA2, color = Ecosystem, shape = Ecosystem),
             size = 3, 
             alpha = 0.7) +
  # Centroids
  geom_point(data = pcoa_centroids, 
             aes(x = PCoA1, y = PCoA2, color = Ecosystem),
             size = 6, 
             shape = 3, 
             stroke = 2) +
  # Scales
  scale_color_manual(values = comp_colors) +
  scale_fill_manual(values = comp_colors) +
  scale_shape_manual(values = c("Abandoned pastures" = 16, "Reference forest" = 17)) +
  # Labels
  labs(
    x = paste0("PCoA 1 (", comp_var_explained[1], "%)"),
    y = paste0("PCoA 2 (", comp_var_explained[2], "%)"),
    color = NULL,
    fill = NULL,
    shape = NULL
  ) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  # Theme
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    legend.key.size = unit(1.2, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed(ratio = 1)

# Save figure
ggsave(paste0(output_path, "Figure_S4_PCoA_composition.tiff"), 
       fig_pcoa,
       width = 16, height = 14, units = "cm", dpi = 300, compression = "lzw")

cat("Figure saved:", paste0(output_path, "Figure_S5_PCoA_composition.tiff"), "\n")

# ==============================================================================
# 7. SUMMARY RESULTS TABLE
# ==============================================================================

cat("\n================================================================================\n")
cat("7. SUMMARY RESULTS\n")
cat("================================================================================\n\n")

results_summary <- data.frame(
  Analysis = c(
    "--- Species Pool ---",
    "Reference forest species",
    "Abandoned pastures species",
    "Shared species",
    "Exclusive to pastures",
    "% shared (of pasture flora)",
    "% exclusive (of pasture flora)",
    "--- PERMANOVA ---",
    "F statistic",
    "R² (variance explained)",
    "p-value",
    "--- Beta Partitioning ---",
    "Total dissimilarity (Jaccard)",
    "Turnover component",
    "Nestedness component",
    "% Turnover",
    "% Nestedness",
    "--- Interpretation ---",
    "Dominant process"
  ),
  Value = c(
    "",
    length(spp_forest),
    length(spp_pastures),
    length(spp_shared),
    length(spp_only_pastures),
    round(length(spp_shared)/length(spp_pastures)*100, 1),
    round(length(spp_only_pastures)/length(spp_pastures)*100, 1),
    "",
    round(permanova_F, 3),
    round(permanova_R2, 3),
    ifelse(permanova_p < 0.001, "< 0.001", round(permanova_p, 3)),
    "",
    round(beta_total, 3),
    round(beta_turnover, 3),
    round(beta_nestedness, 3),
    round(pct_turnover, 1),
    round(pct_nestedness, 1),
    "",
    ifelse(pct_turnover > pct_nestedness, "Species turnover", "Nestedness")
  )
)

print(results_summary, row.names = FALSE)

write.csv(results_summary, 
          "outputs/tables/Compositional_Analysis_Results.csv", 
          row.names = FALSE)

cat("\nResults saved: outputs/tables/Compositional_Analysis_Results.csv\n")


# ==============================================================================
# END OF SCRIPT
# ==============================================================================
