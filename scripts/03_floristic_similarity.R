# ==============================================================================
# FLORISTIC SIMILARITY ANALYSIS
# ==============================================================================
#
# Description: Analyzes floristic similarity between regenerating communities
#              and adjacent reference forests using Jaccard index. Includes 
#              GAMM modeling with temporal autocorrelation correction (CAR1)
#              and internal similarity calculation within reference forests.
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
# Input:       dados/dados_2025_atualizado.csv (regeneration data, via df_processed)
#              dados/HIGUCHI-SILVA LABDENDRO DATABASE  2025.csv (reference forest)
#
# Output:      ../artigo/biotropica/review/Figure_4_similarity.tiff
#              ../artigo/biotropica/review/Figure_S4_similarity_diagnostics.tiff
#              outputs/tables/Similarity_Results_[INVENTORY].csv
#              outputs/tables/Internal_Similarity_Reference_Forests_[INVENTORY].csv
#
# Note:        Requires df_processed from 01_data_processing.R
#              R project located in MariaJulia/dissertacao
#
# ==============================================================================

# --- Setup --------------------------------------------------------------------

set.seed(123)

# Required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(nlme)
  library(vegan)
  library(tibble)
})

# Output directories
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)

# Output path for publication figures
output_path <- "../artigo/biotropica/review/"

# ==============================================================================
# USER CONFIGURATION - SELECT REFERENCE INVENTORY YEAR
# ==============================================================================

#' Choose which inventory to use for reference forest composition:
#' Options: "CAP1", "CAP2", or "CAP3"
#'   CAP1 = First inventory (ANO1: 2016/2017)
#'   CAP2 = Second inventory (ANO2: 2020)
#'   CAP3 = Third inventory (ANO3: 2024)

REF_INVENTORY_YEAR <- "CAP2"  # <-- CHANGE THIS TO SELECT DIFFERENT YEAR

# ==============================================================================
# 1. DATA PREPARATION - REFERENCE FORESTS
# ==============================================================================

cat("================================================================================\n")
cat("SIMILARITY ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

cat("CONFIGURATION:\n")
cat("--------------\n")
cat("Reference inventory selected:", REF_INVENTORY_YEAR, "\n\n")

# Load reference data from LABDENDRO database
ref_database_path <- "dados/HIGUCHI-SILVA LABDENDRO DATABASE  2025.csv"

if(file.exists(ref_database_path)) {
  ref_raw_data <- read.table(ref_database_path, header = TRUE, dec = ",", sep = ";",
                             fileEncoding = "latin1", stringsAsFactors = FALSE)
  
  cat("Reference database loaded successfully\n")
  cat("Total rows:", nrow(ref_raw_data), "\n")
  cat("Columns:", paste(names(ref_raw_data), collapse = ", "), "\n\n")
  
} else {
  stop(paste("File not found:", ref_database_path))
}

# ==============================================================================
# 2. PROCESS REFERENCE FOREST DATA
# ==============================================================================

cat("PROCESSING REFERENCE FOREST DATA\n")
cat("---------------------------------\n")

# Map selected inventory to column names
ref_cap_column <- REF_INVENTORY_YEAR
ref_year_column <- switch(REF_INVENTORY_YEAR,
                          "CAP1" = "ANO1",
                          "CAP2" = "ANO2",
                          "CAP3" = "ANO3")

cat("Using column:", ref_cap_column, "\n")
cat("Year column:", ref_year_column, "\n\n")

#' Check if CAP value is valid (non-empty and represents measurement)
#' Handles cases with multiple stems (e.g., "15+12+10")
is_valid_cap <- function(x) {
  if(is.na(x) || x == "" || x == " ") return(FALSE)
  values <- strsplit(as.character(x), "\\+")[[1]]
  values <- trimws(values)
  values <- values[values != ""]
  if(length(values) == 0) return(FALSE)
  first_val <- suppressWarnings(as.numeric(gsub(",", ".", values[1])))
  return(!is.na(first_val) && first_val > 0)
}

# Filter reference data based on selected inventory
ref_dados_floresta <- ref_raw_data %>%
  rowwise() %>%
  mutate(cap_valid = is_valid_cap(!!sym(ref_cap_column))) %>%
  ungroup() %>%
  filter(cap_valid == TRUE, T %in% c(1, 2, 3, 4)) %>%
  select(Area = T, 
         Parcela = P, 
         Especie = SPP, 
         CAP = !!sym(ref_cap_column),
         Ano_Inventario = !!sym(ref_year_column),
         Estagio = ESTAGIO_SUCESSIONAL_AREA) %>%
  filter(Especie != "NI")  # Remove unidentified species

ref_inventory_years <- unique(ref_dados_floresta$Ano_Inventario)

cat("Reference forest data after filtering (", ref_cap_column, " > 0):\n", sep = "")
cat("Total individuals:", nrow(ref_dados_floresta), "\n")
cat("Species:", n_distinct(ref_dados_floresta$Especie), "\n")
cat("Inventory year(s):", paste(ref_inventory_years, collapse = ", "), "\n\n")

# Summary by area
ref_summary_by_area <- ref_dados_floresta %>%
  group_by(Area, Estagio) %>%
  summarise(
    n_individuals = n(),
    n_species = n_distinct(Especie),
    n_plots = n_distinct(Parcela),
    .groups = 'drop'
  )

cat("Summary by area:\n")
print(ref_summary_by_area, n = Inf)
cat("\n")

# ==============================================================================
# 3. CALCULATE JACCARD SIMILARITY (Regeneration vs Reference)
# ==============================================================================

cat("CALCULATING FLORISTIC SIMILARITY\n")
cat("---------------------------------\n")

#' Calculate Jaccard similarity index
#' 
#' Compares species composition between regenerating community and 
#' adjacent reference forest for a given area and year
#' 
#' @param area_alvo Target area (1-4)
#' @param ano_alvo Target year
#' @param dados_regeneracao Regeneration data (df_processed)
#' @param dados_referencia Reference forest data
#' @return Jaccard similarity index (0-1)

calcular_jaccard_adapted <- function(area_alvo, ano_alvo, dados_regeneracao, dados_referencia) {
  
  # Reference species for target area
  ref_especies_ref <- dados_referencia %>%
    filter(Area == area_alvo) %>%
    pull(Especie) %>%
    unique()
  
  if(length(ref_especies_ref) == 0) {
    warning(paste("No reference species found for area", area_alvo))
    return(NA)
  }
  
  # Regenerating species for target year
  nome_coluna <- paste0("h", ano_alvo)
  
  if (!(nome_coluna %in% names(dados_regeneracao))) {
    return(NA)
  }
  
  ref_especies_reg <- dados_regeneracao %>%
    filter(area == area_alvo) %>%
    filter(!is.na(!!sym(nome_coluna))) %>%
    filter(!!sym(nome_coluna) > 0) %>%
    pull(especies) %>%
    unique()
  
  if (length(ref_especies_reg) == 0) {
    return(0)
  }
  
  # Calculate Jaccard index
  intersecao <- length(intersect(ref_especies_ref, ref_especies_reg))
  uniao <- length(union(ref_especies_ref, ref_especies_reg))
  
  if (uniao == 0) return(0)
  
  return(intersecao / uniao)
}

# Check dependency
if(!exists("df_processed")) {
  stop("Object 'df_processed' not found. Please run 01_data_processing.R first.")
}

# Extract years from data columns
ref_height_cols <- grep("^h\\d{4}$", names(df_processed), value = TRUE)
ref_anos <- as.numeric(gsub("h", "", ref_height_cols))
ref_areas <- c(1, 2, 3, 4)

# Calculate for all combinations
ref_similarity_data <- expand.grid(Area = ref_areas, Ano = ref_anos) %>%
  rowwise() %>%
  mutate(Similaridade = calcular_jaccard_adapted(Area, Ano, df_processed, ref_dados_floresta)) %>%
  filter(!is.na(Similaridade))

cat("Similarity data calculated:\n")
cat("Total observations:", nrow(ref_similarity_data), "\n\n")

# Separate datasets for analysis
# Areas 1, 2, 4: Abandoned pastures (anthropogenic grasslands)
# Area 3: Natural grassland (distinct ecosystem)
ref_antrop_similarity <- ref_similarity_data %>%
  filter(Area %in% c(1, 2, 4)) %>%
  mutate(Area_Factor = factor(Area))

ref_natural_similarity <- ref_similarity_data %>%
  filter(Area == 3)

ref_primeiro_ano <- min(ref_similarity_data$Ano)
ref_ultimo_ano <- max(ref_similarity_data$Ano)

cat("DATA SUMMARY\n")
cat("------------\n")
cat("Abandoned pastures (Areas 1,2,4): n =", nrow(ref_antrop_similarity), "\n")
cat("Natural grassland (Area 3): n =", nrow(ref_natural_similarity), "\n")
cat("Period:", ref_primeiro_ano, "-", ref_ultimo_ano, "\n")
cat("Mean similarity (abandoned pastures):", round(mean(ref_antrop_similarity$Similaridade), 3), "\n")
cat("Mean similarity (natural grassland):", round(mean(ref_natural_similarity$Similaridade), 3), "\n\n")

# ==============================================================================
# 4. STATISTICAL MODELS (GAMM + AIC Selection)
# ==============================================================================

cat("MODEL FITTING & SELECTION\n")
cat("-------------------------\n")

ref_k_value <- min(length(unique(ref_antrop_similarity$Ano)) - 1, 4)

# --- A. Baseline Model (No Temporal Correlation) ---
cat("1. Fitting Baseline Model (No Temporal Correlation)...\n")
ref_gamm_base <- gamm(
  Similaridade ~ s(Ano, k = ref_k_value),
  random = list(Area_Factor = ~1),
  data = ref_antrop_similarity,
  method = "REML"
)

# --- B. CAR1 Model (With Temporal Correlation) ---
cat("2. Fitting CAR1 Model (With Temporal Correlation)...\n")
ref_cor_struct <- corCAR1(form = ~ Ano | Area_Factor)

ref_gamm_cor <- gamm(
  Similaridade ~ s(Ano, k = ref_k_value),
  random = list(Area_Factor = ~1),
  correlation = ref_cor_struct,
  data = ref_antrop_similarity,
  method = "REML"
)

# --- C. Model Comparison (AIC) ---
ref_aic_base <- AIC(ref_gamm_base$lme)
ref_aic_cor <- AIC(ref_gamm_cor$lme)
ref_delta_aic <- ref_aic_base - ref_aic_cor

cat("\nMODEL COMPARISON (AIC):\n")
cat("Baseline AIC:", round(ref_aic_base, 3), "\n")
cat("CAR1 AIC:    ", round(ref_aic_cor, 3), "\n")
cat("Delta AIC:   ", round(ref_delta_aic, 3), "\n")

# Select best model
if(ref_delta_aic > 2) {
  cat("-> CONCLUSION: CAR1 model justified.\n")
  ref_best_model <- ref_gamm_cor
  ref_used_cor <- TRUE
} else {
  cat("-> CONCLUSION: CAR1 not significantly better. Using baseline.\n")
  ref_best_model <- ref_gamm_base
  ref_used_cor <- FALSE
}

# --- D. Heteroscedasticity Check ---
ref_residuals_temp <- resid(ref_best_model$lme, type = "normalized")
ref_fitted_temp <- fitted(ref_best_model$lme)
ref_bp_model <- lm(ref_residuals_temp^2 ~ ref_fitted_temp)
ref_bp_pval <- pchisq(length(ref_residuals_temp) * summary(ref_bp_model)$r.squared, 
                      df = 1, lower.tail = FALSE)

cat("\nHeteroscedasticity check (p-value):", round(ref_bp_pval, 4), "\n")

ref_final_model <- ref_best_model
ref_heterosced_corrected <- FALSE

if (ref_bp_pval < 0.05) {
  cat("Heteroscedasticity detected. Adding varIdent...\n")
  if(ref_used_cor) {
    ref_final_model <- gamm(
      Similaridade ~ s(Ano, k = ref_k_value),
      random = list(Area_Factor = ~1),
      correlation = ref_cor_struct,
      weights = varIdent(form = ~ 1 | Area_Factor),
      data = ref_antrop_similarity,
      method = "REML"
    )
  } else {
    ref_final_model <- gamm(
      Similaridade ~ s(Ano, k = ref_k_value),
      random = list(Area_Factor = ~1),
      weights = varIdent(form = ~ 1 | Area_Factor),
      data = ref_antrop_similarity,
      method = "REML"
    )
  }
  ref_heterosced_corrected <- TRUE
}

# Linear model for reporting slope
ref_lm_sim <- lm(Similaridade ~ Ano, data = ref_antrop_similarity)

# ==============================================================================
# 5. EXTRACT RESULTS
# ==============================================================================

cat("\nEXTRACTING RESULTS\n")
cat("------------------\n")

ref_final_summary <- summary(ref_final_model$gam)
ref_final_smooth <- ref_final_summary$s.table["s(Ano)",]
ref_residuals_final <- resid(ref_final_model$lme, type = "normalized")

# Phi value
ref_phi_value <- NA
if(ref_used_cor) {
  ref_cs <- ref_final_model$lme$modelStruct$corStruct
  ref_phi_value <- as.numeric(coef(ref_cs, unconstrained = FALSE))
} else {
  ref_phi_value <- 0
}

# Slope & Increment
ref_slope_linear <- coef(ref_lm_sim)[2]
ref_se_slope <- summary(ref_lm_sim)$coefficients[2, 2]

ref_pred_ini <- predict(ref_final_model$gam, newdata = data.frame(Ano = ref_primeiro_ano, Area_Factor = "1"))
ref_pred_fin <- predict(ref_final_model$gam, newdata = data.frame(Ano = ref_ultimo_ano, Area_Factor = "1"))
ref_incremento_pct <- if(ref_pred_ini > 0.001) ((ref_pred_fin - ref_pred_ini) / ref_pred_ini) * 100 else 0

# Diagnostics
ref_shapiro_test <- shapiro.test(ref_residuals_final)

ref_fitted_final <- fitted(ref_final_model$lme)
ref_bp_final_model <- lm(ref_residuals_final^2 ~ ref_fitted_final)
ref_bp_pval_final <- pchisq(length(ref_residuals_final) * summary(ref_bp_final_model)$r.squared, 
                            df = 1, lower.tail = FALSE)

ref_acf_result <- acf(ref_residuals_final, plot = FALSE)
ref_lags_sig <- which(abs(ref_acf_result$acf[-1]) > qnorm(0.975)/sqrt(length(ref_residuals_final)))
ref_status_acf <- ifelse(length(ref_lags_sig) == 0, "OK", paste(length(ref_lags_sig), "lags sig."))

# ==============================================================================
# 6. RESULTS TABLE
# ==============================================================================

cat("\n=== DETAILED RESULTS (SIMILARITY - ", REF_INVENTORY_YEAR, ") ===\n\n", sep = "")

ref_results_table <- data.frame(
  Parameter = c(
    "Reference inventory",
    "Temporal trend (F)",
    "p-value",
    "Trend type (EDF)",
    "Adjusted R²",
    "Rate of change (Jaccard/year)",
    "Total increment (%)",
    "Baseline Model (AIC)",
    "CAR1 Model (AIC)",
    "Delta AIC",
    "Phi (Autocorrelation)",
    "Heteroscedasticity corrected",
    "Normality (Shapiro p)",
    "Homoscedasticity (BP p)",
    "Residual autocorrelation status"
  ),
  
  Value = c(
    paste(REF_INVENTORY_YEAR, "(", paste(ref_inventory_years, collapse = "/"), ")"),
    round(ref_final_smooth["F"], 2),
    formatC(ref_final_smooth["p-value"], format = "e", digits = 2),
    paste0(ifelse(ref_final_smooth["edf"] > 1.5, "Non-linear", "Linear"), " (", round(ref_final_smooth["edf"], 2), ")"),
    round(ref_final_summary$r.sq, 3),
    paste0(format(ref_slope_linear, digits = 3), " ± ", format(ref_se_slope, digits = 3)),
    round(ref_incremento_pct, 1),
    round(ref_aic_base, 1),
    round(ref_aic_cor, 1),
    round(ref_delta_aic, 2),
    round(ref_phi_value, 3),
    ifelse(ref_heterosced_corrected, "Yes", "No"),
    round(ref_shapiro_test$p.value, 3),
    round(ref_bp_pval_final, 3),
    ref_status_acf
  ),
  
  Interpretation = c(
    paste("Year:", paste(ref_inventory_years, collapse = "/")),
    ifelse(ref_final_smooth["p-value"] < 0.05, "Significant", "Stable/No trend"),
    ifelse(ref_final_smooth["p-value"] < 0.01, "**", "NS"),
    "Curve complexity",
    "Model fit",
    "Compositional change",
    "Variation over period",
    "Without temporal correction",
    "With temporal correction",
    ifelse(ref_delta_aic > 2, "CAR1 Justified", "CAR1 Indifferent"),
    ifelse(ref_used_cor, "Dependency strength", "Not applied"),
    "Variance weights",
    ifelse(ref_shapiro_test$p.value > 0.05, "Normal", "Deviation"),
    ifelse(ref_bp_pval_final > 0.05, "Homogeneous", "Heterogeneous"),
    ifelse(length(ref_lags_sig) == 0, "Success", "Alert")
  ),
  check.names = FALSE
)

print(ref_results_table, row.names = FALSE)

# Save with inventory year in filename
ref_output_suffix <- paste0("_", REF_INVENTORY_YEAR)
write.csv(ref_results_table, 
          paste0("outputs/tables/Similarity_Results", ref_output_suffix, ".csv"), 
          row.names = FALSE)

# ==============================================================================
# 7. INTERNAL SIMILARITY WITHIN REFERENCE FORESTS
# ==============================================================================

cat("\n=== INTERNAL SIMILARITY WITHIN REFERENCE FORESTS (", REF_INVENTORY_YEAR, ") ===\n", sep = "")
cat("----------------------------------------------------\n")

#' Calculate mean pairwise Jaccard similarity within an area
#' 
#' @param dados Data frame with columns: Area, Parcela, Especie
#' @param area_alvo Target area
#' @return Data frame with similarity statistics

calcular_similaridade_interna_adapted <- function(dados, area_alvo) {
  
  dados_area <- dados %>% filter(Area == area_alvo)
  
  # Create presence/absence matrix (plots x species)
  matriz_pa <- dados_area %>%
    mutate(presenca = 1) %>%
    select(Parcela, Especie, presenca) %>%
    distinct() %>%
    pivot_wider(
      names_from = Especie,
      values_from = presenca,
      values_fill = 0
    ) %>%
    column_to_rownames("Parcela")
  
  cat(sprintf("\nArea %d: %d plots x %d species\n", 
              area_alvo, nrow(matriz_pa), ncol(matriz_pa)))
  
  if(nrow(matriz_pa) < 2) {
    return(data.frame(
      Area = area_alvo,
      n_plots = nrow(matriz_pa),
      n_species = ncol(matriz_pa),
      mean_jaccard = NA,
      sd_jaccard = NA,
      min_jaccard = NA,
      max_jaccard = NA,
      n_pairs = 0
    ))
  }
  
  # Calculate Jaccard distance and convert to similarity
  dist_jaccard <- vegdist(matriz_pa, method = "jaccard", binary = TRUE)
  sim_jaccard <- 1 - as.matrix(dist_jaccard)
  pairwise_values <- sim_jaccard[upper.tri(sim_jaccard)]
  
  data.frame(
    Area = area_alvo,
    n_plots = nrow(matriz_pa),
    n_species = ncol(matriz_pa),
    mean_jaccard = mean(pairwise_values),
    sd_jaccard = sd(pairwise_values),
    min_jaccard = min(pairwise_values),
    max_jaccard = max(pairwise_values),
    n_pairs = length(pairwise_values)
  )
}

# Calculate for all reference forest areas
ref_areas_list <- c(1, 2, 3, 4)
ref_resultados_internos <- lapply(ref_areas_list, function(a) {
  calcular_similaridade_interna_adapted(ref_dados_floresta, a)
})

ref_tabela_similaridade_interna <- do.call(rbind, ref_resultados_internos)

cat("\n")
print(ref_tabela_similaridade_interna, row.names = FALSE)

# Overall statistics
ref_media_geral <- mean(ref_tabela_similaridade_interna$mean_jaccard, na.rm = TRUE)
ref_sd_geral <- sd(ref_tabela_similaridade_interna$mean_jaccard, na.rm = TRUE)
ref_se_geral <- ref_sd_geral / sqrt(sum(!is.na(ref_tabela_similaridade_interna$mean_jaccard)))

cat(sprintf("\nOverall mean internal similarity: %.3f ± %.3f (SE)\n", ref_media_geral, ref_se_geral))

ref_range_min <- min(ref_tabela_similaridade_interna$mean_jaccard, na.rm = TRUE)
ref_range_max <- max(ref_tabela_similaridade_interna$mean_jaccard, na.rm = TRUE)
cat(sprintf("Range: %.3f - %.3f\n", ref_range_min, ref_range_max))

# ==============================================================================
# 8. COMPARISON SUMMARY
# ==============================================================================

cat("\n=== COMPARISON SUMMARY ===\n")
cat("--------------------------\n")

ref_mean_regen_similarity <- mean(ref_antrop_similarity$Similaridade)

cat(sprintf("Internal similarity (reference forests): %.3f (range: %.3f - %.3f)\n", 
            ref_media_geral, ref_range_min, ref_range_max))
cat(sprintf("Regeneration vs reference: %.3f\n", ref_mean_regen_similarity))
cat(sprintf("Ratio: %.1f%%\n", (ref_mean_regen_similarity / ref_media_geral) * 100))

write.csv(ref_tabela_similaridade_interna, 
          paste0("outputs/tables/Internal_Similarity_Reference_Forests", ref_output_suffix, ".csv"), 
          row.names = FALSE)

# ==============================================================================
# 9. MAIN FIGURE
# ==============================================================================

cat("\nGENERATING MAIN FIGURE\n")
cat("----------------------\n")

# Prediction grid
ref_anos_pred <- seq(ref_primeiro_ano, ref_ultimo_ano, length.out = 100)

# Anthropogenic predictions
ref_pred_antrop <- data.frame(Ano = ref_anos_pred, Area_Factor = "1")
ref_pred_list <- predict(ref_final_model$gam, newdata = ref_pred_antrop, se.fit = TRUE)

ref_pred_antrop$fit <- ref_pred_list$fit
ref_pred_antrop$se.fit <- ref_pred_list$se.fit
ref_pred_antrop$lower <- ref_pred_antrop$fit - 1.96 * ref_pred_antrop$se.fit
ref_pred_antrop$upper <- ref_pred_antrop$fit + 1.96 * ref_pred_antrop$se.fit
ref_pred_antrop$Grupo <- "Abandoned pastures"

# Natural grassland (flat at observed mean)
ref_mean_natural <- mean(ref_natural_similarity$Similaridade, na.rm = TRUE)
ref_pred_natural <- data.frame(
  Ano = ref_anos_pred,
  fit = ref_mean_natural,
  se.fit = 0,
  lower = ref_mean_natural,
  upper = ref_mean_natural,
  Grupo = "Natural grassland"
)

# Annual means for abandoned pastures
ref_annual_means <- ref_antrop_similarity %>%
  group_by(Ano) %>%
  summarise(
    mean = mean(Similaridade),
    se = sd(Similaridade)/sqrt(n()),
    n = n(),
    .groups = 'drop'
  )
ref_annual_means$Grupo <- "Abandoned pastures"

# Natural grassland points
ref_natural_points <- ref_natural_similarity %>%
  mutate(Grupo = "Natural grassland")

# Create plot
ref_fig_similarity <- ggplot() +
  # Ribbons
  geom_ribbon(data = ref_pred_antrop, 
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo), 
              alpha = 0.22) +
  
  # Lines
  geom_line(data = ref_pred_antrop, 
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), 
            linewidth = 2.0) +
  geom_line(data = ref_pred_natural, 
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), 
            linewidth = 1.4) +
  
  # Points and error bars
  geom_errorbar(data = ref_annual_means, 
                aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.35, alpha = 0.7, linewidth = 0.8, show.legend = FALSE) +
  geom_point(data = ref_annual_means, 
             aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 3.6, stroke = 1.2) +
  geom_point(data = ref_natural_points, 
             aes(x = Ano, y = Similaridade, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.8, alpha = 0.7) +
  
  # Reference line at 0.5
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50", alpha = 0.5) +
  
  # Scales
  scale_color_manual(values = c("Abandoned pastures" = "#D55E00", "Natural grassland" = "#009E73")) +
  scale_fill_manual(values = c("Abandoned pastures" = "#D55E00", "Natural grassland" = "#009E73")) +
  scale_linetype_manual(values = c("Abandoned pastures" = "solid", "Natural grassland" = "23")) +
  scale_shape_manual(values = c("Abandoned pastures" = 21, "Natural grassland" = 23)) +
  
  # Layout
  scale_x_continuous(breaks = seq(floor(ref_primeiro_ano/2)*2, 2025, by = 2),
                     limits = c(ref_primeiro_ano, 2025), 
                     expand = expansion(mult = c(0.02, 0.08))) +
  scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, by = 0.2), 
                     expand = expansion(mult = c(0.02, 0.02))) +
  
  labs(x = "Year", 
       y = "Jaccard similarity to reference forest", 
       color = NULL, fill = NULL, linetype = NULL, shape = NULL) +
  
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 10, face = "italic", hjust = 1),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "bottom",
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.height = unit(0.9, "lines"),
    legend.key.width = unit(1.6, "lines")
  )

ggsave(paste0(output_path, "Figure_4_similarity.tiff"), 
       ref_fig_similarity,
       width = 14, height = 10, units = "cm", dpi = 300, compression = "lzw")
cat("Main figure saved:", paste0(output_path, "Figure_4_similarity.tiff"), "\n")

# ==============================================================================
# 10. DIAGNOSTIC FIGURE
# ==============================================================================

cat("GENERATING DIAGNOSTIC FIGURE\n")

tiff(paste0(output_path, "Figure_S4_similarity_diagnostics.tiff"),
     width = 2400, height = 1800, res = 300, compression = "lzw")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2), cex.lab = 1.1, cex.axis = 1)

plot(ref_fitted_final, ref_residuals_final,
     xlab = "Fitted values", ylab = "Standardized residuals",
     main = "A. Residuals vs Fitted", pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
abline(h = 0, col = "red", lty = 2, lwd = 2)
lines(lowess(ref_fitted_final, ref_residuals_final), col = "blue", lwd = 2)

qqnorm(ref_residuals_final, main = "B. Normal Q-Q Plot", pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
qqline(ref_residuals_final, col = "red", lwd = 2)

hist(ref_residuals_final, breaks = 15, probability = TRUE,
     main = "C. Distribution of Residuals",
     xlab = "Standardized residuals", col = "lightgray", border = "darkgray")
curve(dnorm(x, mean = mean(ref_residuals_final), sd = sd(ref_residuals_final)), 
      col = "red", lwd = 2, add = TRUE)

acf(ref_residuals_final, main = "D. Autocorrelation Function", lwd = 2)

dev.off()
cat("Diagnostic figure saved:", paste0(output_path, "Figure_S4_similarity_diagnostics.tiff"), "\n")


cat("\n================================================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("Reference inventory used:", REF_INVENTORY_YEAR, "\n")
cat("================================================================================\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================

