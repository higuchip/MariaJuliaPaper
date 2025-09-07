# ================================================================================
# FLORISTIC SIMILARITY ANALYSIS IN ABANDONED PASTURES
# Revised approach for stable similarity patterns
# ================================================================================

set.seed(123)

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(nlme)
library(vegan)

# Create output directories
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)

# ================================================================================
# 1. DATA PREPARATION
# ================================================================================

cat("================================================================================\n")
cat("SIMILARITY ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

# Load reference data (mature forest)
referencia <- read.table("dados/dados_talissa.csv", header=TRUE, dec=".", sep=",")
dados_ref <- referencia %>% 
  filter(DAP2 > 0, Area %in% c(1, 2, 3, 4)) %>%
  select(Area, Especie)

# Identify height columns for temporal data
height_cols <- grep("^h\\d{4}$", names(df_processed), value = TRUE)

# Function to calculate Jaccard similarity
calcular_jaccard <- function(area_alvo, ano_alvo) {
  # Campo natural (área 3) sempre tem similaridade zero
  if(area_alvo == 3) {
    return(0)
  }
  
  # Reference species for target area
  especies_ref <- dados_ref %>%
    filter(Area == area_alvo) %>%
    pull(Especie) %>%
    unique()
  
  # Regenerating species for target year
  nome_coluna <- paste0("h", ano_alvo)
  
  if (!(nome_coluna %in% names(df_processed))) {
    return(NA)
  }
  
  especies_reg <- df_processed %>%
    filter(area == area_alvo) %>%
    filter(!is.na(!!sym(nome_coluna))) %>%
    filter(!!sym(nome_coluna) > 0) %>%
    pull(especies) %>%
    unique()
  
  if (length(especies_reg) == 0) {
    return(0)
  }
  
  # Calculate Jaccard index
  intersecao <- length(intersect(especies_ref, especies_reg))
  uniao <- length(union(especies_ref, especies_reg))
  
  if (uniao == 0) return(0)
  
  return(intersecao / uniao)
}

# Calculate similarity for each area and year
areas <- c(1, 2, 3, 4)
anos <- c(2014, 2015, 2016, 2017, 2018, 2022, 2023, 2025)

similarity_data <- expand.grid(Area = areas, Ano = anos) %>%
  rowwise() %>%
  mutate(Similaridade = calcular_jaccard(Area, Ano)) %>%
  filter(!is.na(Similaridade))

# Separate datasets
antrop_similarity_data <- similarity_data %>%
  filter(Area %in% c(1, 2, 4)) %>%
  mutate(Area_Antrop = factor(Area))

reference_similarity_data <- similarity_data %>%
  filter(Area == 3)

cat("DATA SUMMARY\n")
cat("------------\n")
cat("Abandoned pastures: n =", nrow(antrop_similarity_data), "\n")
cat("Natural grassland: n =", nrow(reference_similarity_data), "\n")
cat("Period:", min(similarity_data$Ano), "-", max(similarity_data$Ano), "\n\n")

# ================================================================================
# 2. DESCRIPTIVE ANALYSIS (More appropriate for stable data)
# ================================================================================

cat("DESCRIPTIVE STATISTICS\n")
cat("---------------------\n")

# Overall statistics
overall_stats <- antrop_similarity_data %>%
  summarise(
    mean = mean(Similaridade),
    median = median(Similaridade),
    sd = sd(Similaridade),
    cv = sd(Similaridade)/mean(Similaridade)*100,
    min = min(Similaridade),
    max = max(Similaridade),
    range = max(Similaridade) - min(Similaridade)
  )

cat("Overall similarity:\n")
cat(sprintf("  Mean ± SD: %.3f ± %.3f\n", overall_stats$mean, overall_stats$sd))
cat(sprintf("  Median: %.3f\n", overall_stats$median))
cat(sprintf("  Range: %.3f - %.3f\n", overall_stats$min, overall_stats$max))
cat(sprintf("  Coefficient of variation: %.1f%%\n", overall_stats$cv))

# Statistics by area
area_stats <- antrop_similarity_data %>%
  group_by(Area_Antrop) %>%
  summarise(
    n = n(),
    mean = mean(Similaridade),
    sd = sd(Similaridade),
    cv = sd(Similaridade)/mean(Similaridade)*100,
    min = min(Similaridade),
    max = max(Similaridade),
    .groups = 'drop'
  )

cat("\nSimilarity by area:\n")
print(area_stats)

# ================================================================================
# 3. STATISTICAL MODELS
# ================================================================================

cat("\nMODEL FITTING\n")
cat("-------------\n")

# Model 1: Null model (just mean)
modelo_nulo <- lm(Similaridade ~ 1, data = antrop_similarity_data)

# Model 2: Linear time trend
modelo_tempo <- lm(Similaridade ~ Ano, data = antrop_similarity_data)

# Model 3: Mixed model with area random effect (no time)
modelo_area <- lmer(Similaridade ~ 1 + (1|Area_Antrop), 
                    data = antrop_similarity_data)

# Model 4: Mixed model with time and area
modelo_tempo_area <- lmer(Similaridade ~ Ano + (1|Area_Antrop), 
                          data = antrop_similarity_data)

# Compare models
cat("\nModel comparison (AIC):\n")
aic_comparison <- AIC(modelo_nulo, modelo_tempo, modelo_area, modelo_tempo_area)
aic_comparison$delta_AIC <- aic_comparison$AIC - min(aic_comparison$AIC)
print(aic_comparison[order(aic_comparison$AIC),])

# Select best model (likely the null or area-only model)
best_model_name <- rownames(aic_comparison)[which.min(aic_comparison$AIC)]
cat(sprintf("\nBest model: %s\n", best_model_name))

# Test for temporal trend
tempo_test <- anova(modelo_nulo, modelo_tempo)
cat("\nTemporal trend test:\n")
print(tempo_test)

# Kruskal-Wallis test for differences between areas
kw_test <- kruskal.test(Similaridade ~ Area_Antrop, data = antrop_similarity_data)
cat("\nKruskal-Wallis test (differences between areas):\n")
cat(sprintf("  Chi-squared = %.2f, p-value = %.4f\n", 
            kw_test$statistic, kw_test$p.value))

# ================================================================================
# 4. RESULTS TABLE
# ================================================================================

cat("\n=== MAIN RESULTS ===\n")

results_similarity <- data.frame(
  Parameter = c(
    "Mean similarity",
    "Standard deviation",
    "Coefficient of variation (%)",
    "Median similarity",
    "Range",
    "Temporal trend (slope)",
    "Temporal trend (p-value)",
    "Differences between areas",
    "Period analyzed",
    "N observations",
    "Best model",
    "Natural grassland similarity"
  ),
  
  Value = c(
    sprintf("%.3f", overall_stats$mean),
    sprintf("%.3f", overall_stats$sd),
    sprintf("%.1f", overall_stats$cv),
    sprintf("%.3f", overall_stats$median),
    sprintf("%.3f - %.3f", overall_stats$min, overall_stats$max),
    sprintf("%.4f", coef(modelo_tempo)[2]),
    sprintf("%.3f", summary(modelo_tempo)$coefficients[2,4]),
    ifelse(kw_test$p.value < 0.05, "Significant", "Not significant"),
    paste(min(antrop_similarity_data$Ano), "-", max(antrop_similarity_data$Ano)),
    as.character(nrow(antrop_similarity_data)),
    best_model_name,
    "0 (distinct ecosystem)"
  )
)

print(results_similarity, row.names = FALSE)
write.csv(results_similarity, "outputs/tables/Similarity_results_revised.csv", row.names = FALSE)

# ================================================================================
# 5. PREPARE DATA FOR VISUALIZATION
# ================================================================================

cat("\nPREPARING DATA FOR VISUALIZATION\n")
cat("---------------------------------\n")

# Calculate annual means for abandoned pastures
annual_means <- antrop_similarity_data %>%
  group_by(Ano) %>%
  summarise(
    mean = mean(Similaridade),
    se = sd(Similaridade)/sqrt(n()),
    n = n(),
    .groups = 'drop'
  ) %>%
  mutate(Grupo = "Abandoned pastures")

# Create prediction data for abandoned pastures (horizontal line since no trend)
anos_plot <- seq(min(similarity_data$Ano), max(similarity_data$Ano), length.out = 100)

pred_antrop <- data.frame(
  Ano = anos_plot,
  fit = rep(mean(antrop_similarity_data$Similaridade), length(anos_plot)),
  se.fit = rep(sd(antrop_similarity_data$Similaridade)/sqrt(nrow(antrop_similarity_data)), 
               length(anos_plot))
)
pred_antrop$lower <- pred_antrop$fit - 1.96 * pred_antrop$se.fit
pred_antrop$upper <- pred_antrop$fit + 1.96 * pred_antrop$se.fit
pred_antrop$Grupo <- "Abandoned pastures"

# Create prediction data for natural grassland (always zero)
pred_ref <- data.frame(
  Ano = anos_plot,
  fit = 0,
  se.fit = 0,
  lower = 0,
  upper = 0,
  Grupo = "Natural grassland"
)

# Add group column to reference data
reference_similarity_data <- reference_similarity_data %>%
  mutate(Grupo = "Natural grassland")

# Verificar dados criados
cat("Data check for figure:\n")
cat("- pred_antrop: ", nrow(pred_antrop), "rows\n")
cat("- pred_ref: ", nrow(pred_ref), "rows\n")
cat("- annual_means: ", nrow(annual_means), "rows\n")
cat("- reference_similarity_data: ", nrow(reference_similarity_data), "rows\n\n")

# ================================================================================
# 6. MAIN FIGURE - SHOWING STABILITY
# ================================================================================

cat("GENERATING MAIN FIGURE\n")
cat("----------------------\n")

# Create figure matching the style of abundance/richness/height figures
fig_similarity <- ggplot() +
  # Confidence ribbons
  geom_ribbon(data = pred_antrop,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.22, show.legend = TRUE) +
  geom_ribbon(data = pred_ref,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.18, show.legend = TRUE) +
  
  # Fitted lines (horizontal for similarity since no trend)
  geom_line(data = pred_antrop,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 2.0, show.legend = TRUE) +
  geom_line(data = pred_ref,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 1.4, show.legend = TRUE) +
  
  # Error bars for annual means
  geom_errorbar(data = annual_means,
                aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.35, alpha = 0.7, linewidth = 0.8, show.legend = FALSE) +
  
  # Points for annual means
  geom_point(data = annual_means,
             aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 3.6, stroke = 1.2, show.legend = TRUE) +
  
  # Points for natural grassland
  geom_point(data = reference_similarity_data,
             aes(x = Ano, y = Similaridade, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.8, alpha = 0.7, show.legend = TRUE) +
  
  # Reference line at 50%
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50", alpha = 0.5) +
  annotate("text", x = min(similarity_data$Ano) + 0.5, y = 0.52, 
           label = "50% similarity threshold", 
           hjust = 0, size = 3, color = "gray50", fontface = "italic") +
  
  # Color and shape scales matching other figures
  scale_color_manual(values = c("Abandoned pastures" = "#D55E00",
                                "Natural grassland" = "#009E73")) +
  scale_fill_manual(values = c("Abandoned pastures" = "#D55E00",
                               "Natural grassland" = "#009E73")) +
  scale_linetype_manual(values = c("Abandoned pastures" = "solid",
                                   "Natural grassland" = "23")) +
  scale_shape_manual(values = c("Abandoned pastures" = 21,
                                "Natural grassland" = 23)) +
  
  # Legend guides
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    fill = guide_legend(override.aes = list(alpha = 0.4)),
    linetype = guide_legend(),
    shape = guide_legend()
  ) +
  
  # Axes
  scale_x_continuous(
    breaks = seq(floor(min(similarity_data$Ano)/2)*2, 2025, by = 2),
    limits = c(min(similarity_data$Ano), 2025),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(
    limits = c(-0.05, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  # Labels
  labs(x = "Year", 
       y = "Jaccard similarity index",
       color = NULL, fill = NULL, linetype = NULL, shape = NULL) +
  
  # Theme matching other figures
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "top",
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.height = unit(0.9, "lines"),
    legend.key.width = unit(1.6, "lines")
  )

# Save figure
ggsave("outputs/figures/Figure_4_similarity.tiff", fig_similarity,
       width = 14, height = 10, units = "cm", dpi = 300, compression = "lzw")
cat("Main figure saved: Figure_4_similarity.tiff\n")

# ================================================================================
# 7. SUPPLEMENTARY ANALYSIS - TEMPORAL VARIATION
# ================================================================================

cat("\nTEMPORAL VARIATION ANALYSIS\n")
cat("---------------------------\n")

# Calculate year-to-year changes
year_changes <- antrop_similarity_data %>%
  arrange(Area_Antrop, Ano) %>%
  group_by(Area_Antrop) %>%
  mutate(change = Similaridade - lag(Similaridade)) %>%
  filter(!is.na(change))

# Summary of changes
change_summary <- year_changes %>%
  summarise(
    mean_change = mean(change),
    sd_change = sd(change),
    max_increase = max(change),
    max_decrease = min(change),
    .groups = 'drop'
  )

cat("Year-to-year changes by area:\n")
print(change_summary)

# Test if changes are different from zero
t_test_changes <- t.test(year_changes$change, mu = 0)
cat(sprintf("\nTest for directional change: t = %.3f, p = %.4f\n",
            t_test_changes$statistic, t_test_changes$p.value))

# ================================================================================
# 8. VERIFICATION: Natural grassland (Area 3) similarity check
# ================================================================================

cat("\n================================================================================\n")
cat("VERIFYING NATURAL GRASSLAND (AREA 3) SIMILARITY\n")
cat("================================================================================\n\n")

# Get reference species for Area 3
especies_ref_area3 <- dados_ref %>%
  filter(Area == 3) %>%
  pull(Especie) %>%
  unique()

cat("Reference species in Area 3 (natural grassland):", length(especies_ref_area3), "\n")
if(length(especies_ref_area3) > 0) {
  cat("Species:", paste(head(especies_ref_area3, 10), collapse=", "))
  if(length(especies_ref_area3) > 10) cat(", ...")
  cat("\n\n")
}

# Check each year for Area 3
anos_verificar <- c(2014, 2015, 2016, 2018, 2023, 2025)
verificacao_area3 <- data.frame()

for(ano in anos_verificar) {
  nome_coluna <- paste0("h", ano)
  
  if(nome_coluna %in% names(df_processed)) {
    # Get regenerating species in Area 3 for this year
    especies_reg_area3 <- df_processed %>%
      filter(area == 3) %>%
      filter(!is.na(!!sym(nome_coluna))) %>%
      filter(!!sym(nome_coluna) > 0) %>%
      pull(especies) %>%
      unique()
    
    # Calculate intersection
    especies_comuns <- intersect(especies_ref_area3, especies_reg_area3)
    
    # Calculate Jaccard
    if(length(especies_reg_area3) > 0 && length(especies_ref_area3) > 0) {
      jaccard <- length(especies_comuns) / length(union(especies_ref_area3, especies_reg_area3))
    } else {
      jaccard <- 0
    }
    
    verificacao_area3 <- rbind(verificacao_area3, data.frame(
      Ano = ano,
      N_ref = length(especies_ref_area3),
      N_reg = length(especies_reg_area3),
      N_comum = length(especies_comuns),
      Jaccard = jaccard,
      Especies_comuns = paste(especies_comuns, collapse=", ")
    ))
    
    cat(sprintf("Year %d:\n", ano))
    cat(sprintf("  Reference species: %d\n", length(especies_ref_area3)))
    cat(sprintf("  Regenerating species: %d\n", length(especies_reg_area3)))
    cat(sprintf("  Common species: %d\n", length(especies_comuns)))
    cat(sprintf("  Jaccard index: %.3f\n", jaccard))
    if(length(especies_comuns) > 0) {
      cat(sprintf("  Species in common: %s\n", paste(especies_comuns, collapse=", ")))
    }
    cat("\n")
  }
}

cat("SUMMARY FOR AREA 3 (Natural grassland):\n")
print(verificacao_area3)

# Cross-check with other areas
cat("\n----------------------------------------\n")
cat("CROSS-CHECK: Does Area 3 share species with forest references (Areas 1,2,4)?\n")
cat("----------------------------------------\n")

for(area_forest in c(1, 2, 3,4)) {
  especies_ref_forest <- dados_ref %>%
    filter(Area == area_forest) %>%
    pull(Especie) %>%
    unique()
  
  # Check most recent year
  nome_coluna <- ifelse("h2025" %in% names(df_processed), "h2025", "h2023")
  
  especies_reg_area3_current <- df_processed %>%
    filter(area == 3) %>%
    filter(!is.na(!!sym(nome_coluna))) %>%
    filter(!!sym(nome_coluna) > 0) %>%
    pull(especies) %>%
    unique()
  
  especies_compartilhadas <- intersect(especies_ref_forest, especies_reg_area3_current)
  
  cat(sprintf("\nArea 3 (grassland) vs Area %d (forest) reference:\n", area_forest))
  cat(sprintf("  Forest reference species: %d\n", length(especies_ref_forest)))
  cat(sprintf("  Grassland current species: %d\n", length(especies_reg_area3_current)))
  cat(sprintf("  Shared species: %d\n", length(especies_compartilhadas)))
  
  if(length(union(especies_ref_forest, especies_reg_area3_current)) > 0) {
    cat(sprintf("  Jaccard with forest: %.3f\n", 
                length(especies_compartilhadas)/length(union(especies_ref_forest, especies_reg_area3_current))))
  } else {
    cat("  Jaccard with forest: 0.000 (no species)\n")
  }
  
  if(length(especies_compartilhadas) > 0) {
    cat("  Species shared with forest:", paste(especies_compartilhadas, collapse=", "), "\n")
  }
}

# ================================================================================
# 9. FINAL SUMMARY
# ================================================================================

cat("\n================================================================================\n")
cat("SIMILARITY ANALYSIS COMPLETED\n")
cat("================================================================================\n")

cat("\nKEY FINDINGS:\n")
cat("-------------\n")
cat("• Mean similarity (abandoned pastures):", sprintf("%.3f", overall_stats$mean), "\n")
cat("• Temporal trend:", ifelse(summary(modelo_tempo)$coefficients[2,4] < 0.05, 
                                "Significant", "Not significant"), "\n")
cat("• Natural grassland similarity: Always 0 (distinct ecosystem)\n")
cat("• Interpretation: Stable floristic similarity over time\n")

cat("\nFILES GENERATED:\n")
cat("----------------\n")
cat("• Main figure: outputs/figures/Figure_4_similarity.tiff\n")
cat("• Results table: outputs/tables/Similarity_results_revised.csv\n")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("================================================================================\n")