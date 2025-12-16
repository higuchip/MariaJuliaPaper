# ==============================================================================
# TEMPORAL DYNAMICS OF ABUNDANCE
# ==============================================================================
#
# Description: Analyzes temporal changes in woody plant abundance across
#              abandoned pastures and natural grassland using GAMMs with
#              temporal autocorrelation correction (CAR1).
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
# Input:       dados/dados_2025_atualizado.csv (via df_processed)
#
# Output:      ../artigo/biotropica/review/Figure_2_abundance.tiff
#              ../artigo/biotropica/review/Figure_S1_abundance_diagnostics.tiff
#              outputs/tables/Abundance_Results_summary.csv
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
})

# Output directories
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)

# Output path for publication figures
output_path <- "../artigo/biotropica/review/"

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

cat("================================================================================\n")
cat("ABUNDANCE ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

# Check dependency
if(!exists("df_processed")) {
  stop("Object 'df_processed' not found. Please run 01_data_processing.R first.")
}

# Identify height columns
height_cols <- grep("^h\\d{4}$", names(df_processed), value = TRUE)

# Convert to long format
abundance_long <- df_processed %>%
  select(area, parc, ni, especies, all_of(height_cols)) %>%
  pivot_longer(
    cols = all_of(height_cols),
    names_to = "Ano",
    values_to = "Altura",
    names_prefix = "h"
  ) %>%
  mutate(Ano = as.numeric(Ano)) %>%
  filter(!is.na(Altura))

# Calculate abundance (counts per area per year)
abundance_data <- abundance_long %>%
  group_by(Area = as.factor(area), Ano) %>%
  summarise(Abundancia = n(), .groups = 'drop')

# Separate datasets by ecosystem type
# Areas 1, 2, 4: Abandoned pastures (anthropogenic grasslands)
# Area 3: Natural grassland
antrop_abundance_data <- abundance_data %>%
  filter(Area %in% c("1", "2", "4")) %>%
  mutate(Area_Antrop = factor(Area, levels = c("1", "2", "4")))

reference_abundance_data <- abundance_data %>%
  filter(Area == "3")

primeiro_ano <- min(abundance_data$Ano)
ultimo_ano <- max(abundance_data$Ano)

cat("DATA SUMMARY\n")
cat("------------\n")
cat("Anthropogenic sites: n =", nrow(antrop_abundance_data), "\n")
cat("Reference site: n =", nrow(reference_abundance_data), "\n")
cat("Period:", primeiro_ano, "-", ultimo_ano, "\n\n")

# ==============================================================================
# 2. STATISTICAL MODELS (GAMM + AIC Selection)
# ==============================================================================

cat("MODEL FITTING & SELECTION\n")
cat("-------------------------\n")

k_antrop <- min(length(unique(antrop_abundance_data$Ano)) - 1, 4)

# --- A. Baseline Model (No Temporal Correlation) ---
cat("1. Fitting Baseline Model (No Temporal Correlation)...\n")
antrop_gamm_base <- gamm(
  Abundancia ~ s(Ano, k = k_antrop),
  random = list(Area_Antrop = ~1),
  data = antrop_abundance_data,
  method = "REML"
)

# --- B. CAR1 Model (With Temporal Correlation) ---
cat("2. Fitting CAR1 Model (With Temporal Correlation)...\n")
cor_struct <- corCAR1(form = ~ Ano | Area_Antrop)

antrop_gamm_cor <- gamm(
  Abundancia ~ s(Ano, k = k_antrop),
  random = list(Area_Antrop = ~1),
  correlation = cor_struct,
  data = antrop_abundance_data,
  method = "REML"
)

# --- C. Model Comparison (AIC) ---
aic_base <- AIC(antrop_gamm_base$lme)
aic_cor  <- AIC(antrop_gamm_cor$lme)
delta_aic <- aic_base - aic_cor

cat("\nMODEL COMPARISON (AIC):\n")
cat("Baseline AIC:", round(aic_base, 3), "\n")
cat("CAR1 AIC:    ", round(aic_cor, 3), "\n")
cat("Delta AIC:   ", round(delta_aic, 3), "\n")

# Select best model
if(delta_aic > 2) {
  cat("-> CONCLUSION: CAR1 model justified.\n")
  best_model <- antrop_gamm_cor
  used_cor <- TRUE
} else {
  cat("-> CONCLUSION: CAR1 not significantly better. Using baseline.\n")
  best_model <- antrop_gamm_base
  used_cor <- FALSE
}

# --- D. Heteroscedasticity Check ---
res_temp <- resid(best_model$lme, type = "normalized")
fit_temp <- fitted(best_model$lme)
bp_model <- lm(res_temp^2 ~ fit_temp)
bp_stat  <- length(res_temp) * summary(bp_model)$r.squared
bp_pval  <- pchisq(bp_stat, df = 1, lower.tail = FALSE)

cat("\nHeteroscedasticity check (p-value):", round(bp_pval, 4), "\n")

final_model <- best_model
heterosced_corrected <- FALSE

if (bp_pval < 0.05) {
  cat("Heteroscedasticity detected. Adding varIdent...\n")
  
  if(used_cor) {
    final_model <- gamm(
      Abundancia ~ s(Ano, k = k_antrop),
      random = list(Area_Antrop = ~1),
      correlation = cor_struct,
      weights = varIdent(form = ~ 1 | Area_Antrop),
      data = antrop_abundance_data,
      method = "REML"
    )
  } else {
    final_model <- gamm(
      Abundancia ~ s(Ano, k = k_antrop),
      random = list(Area_Antrop = ~1),
      weights = varIdent(form = ~ 1 | Area_Antrop),
      data = antrop_abundance_data,
      method = "REML"
    )
  }
  heterosced_corrected <- TRUE
  cat("Model updated with Variance Structure.\n")
}

# Reference model (simple GAM for natural grassland)
ref_gam <- gam(Abundancia ~ s(Ano, k = 3),
               data = reference_abundance_data,
               method = "REML")

# Linear models for rate reporting
lm_antrop <- lm(Abundancia ~ Ano, data = antrop_abundance_data)
lm_ref <- lm(Abundancia ~ Ano, data = reference_abundance_data)

# ==============================================================================
# 3. EXTRACT RESULTS
# ==============================================================================

cat("\nEXTRACTING RESULTS\n")
cat("------------------\n")

final_summary <- summary(final_model$gam)
final_smooth <- final_summary$s.table["s(Ano)",]
residuos_final <- resid(final_model$lme, type = "normalized")

# Phi parameter (temporal autocorrelation)
phi_value <- NA
if(used_cor) {
  cs <- final_model$lme$modelStruct$corStruct
  phi_value <- as.numeric(coef(cs, unconstrained = FALSE))
  cat("Estimated temporal autocorrelation (Phi):", round(phi_value, 3), "\n")
} else {
  phi_value <- 0
  cat("No correlation structure used (Phi = 0)\n")
}

# Rates & Increment
slope_antrop <- coef(lm_antrop)[2]
se_slope_antrop <- summary(lm_antrop)$coefficients[2, 2]
slope_ref <- coef(lm_ref)[2]
se_slope_ref <- summary(lm_ref)$coefficients[2, 2]

pred_inicial <- predict(final_model$gam, 
                        newdata = data.frame(Ano = primeiro_ano, Area_Antrop = "1"))
pred_final <- predict(final_model$gam, 
                      newdata = data.frame(Ano = ultimo_ano, Area_Antrop = "1"))
incremento_pct <- ((pred_final - pred_inicial) / pred_inicial) * 100

# Diagnostics
shapiro_test <- shapiro.test(residuos_final)

fit_final <- fitted(final_model$lme)
bp_final_model <- lm(residuos_final^2 ~ fit_final)
bp_pval_final <- pchisq(length(residuos_final) * summary(bp_final_model)$r.squared, df=1, lower.tail=FALSE)

acf_result <- acf(residuos_final, plot = FALSE)
lags_sig <- which(abs(acf_result$acf[-1]) > qnorm(0.975)/sqrt(length(residuos_final)))
status_acf <- ifelse(length(lags_sig) == 0, "OK", paste(length(lags_sig), "lags sig."))

# ==============================================================================
# 4. RESULTS TABLE
# ==============================================================================

cat("\n=== DETAILED RESULTS (ABUNDANCE) ===\n")

results_simple <- data.frame(
  Parameter = c(
    "Temporal trend (F)",
    "p-value",
    "Trend type (EDF)",
    "Adjusted R²",
    "Rate of change (ind/year)",
    "Total increment (%)",
    "Baseline Model (AIC)",
    "CAR1 Model (AIC)",
    "Delta AIC",
    "Phi (Autocorrelation)",
    "Heteroscedasticity corrected",
    "Normality (Shapiro p)",
    "Homoscedasticity (BP p)",
    "Residual autocorrelation status",
    "Reference rate (ind/year)"
  ),
  
  Value = c(
    round(final_smooth["F"], 2),
    formatC(final_smooth["p-value"], format = "e", digits = 2),
    paste0(ifelse(final_smooth["edf"] > 1.5, "Non-linear", "Linear"), " (", round(final_smooth["edf"], 2), ")"),
    round(final_summary$r.sq, 3),
    paste0(round(slope_antrop, 1), " ± ", round(se_slope_antrop, 1)),
    round(incremento_pct, 1),
    round(aic_base, 1),
    round(aic_cor, 1),
    round(delta_aic, 2),
    round(phi_value, 3),
    ifelse(heterosced_corrected, "Yes", "No"),
    round(shapiro_test$p.value, 3),
    round(bp_pval_final, 3),
    status_acf,
    round(slope_ref, 1)
  ),
  
  Interpretation = c(
    ifelse(final_smooth["p-value"] < 0.05, "Significant", "NS"),
    ifelse(final_smooth["p-value"] < 0.01, "**", "*"),
    "Curve complexity",
    "Model fit",
    "Demographic aggregation",
    "Structural recovery",
    "Without temporal correction",
    "With temporal correction",
    ifelse(delta_aic > 2, "CAR1 Justified", "CAR1 Indifferent"),
    ifelse(used_cor, "Dependency strength", "Not applied"),
    "Variance weights",
    ifelse(shapiro_test$p.value > 0.05, "Normal", "Deviation"),
    ifelse(bp_pval_final > 0.05, "Homogeneous", "Heterogeneous"),
    ifelse(length(lags_sig) == 0, "Success", "Alert"),
    "Stable natural grassland"
  ),
  check.names = FALSE
)

print(results_simple, row.names = FALSE)
write.csv(results_simple, "outputs/tables/Abundance_Results_summary.csv", row.names = FALSE)

# ==============================================================================
# 5. MAIN FIGURE
# ==============================================================================

cat("\nGENERATING MAIN FIGURE\n")
cat("----------------------\n")

# Predictions for anthropogenic sites
anos_pred <- seq(min(antrop_abundance_data$Ano), max(antrop_abundance_data$Ano), length.out = 100)
pred_antrop <- data.frame(Ano = anos_pred)

pred_antrop_list <- predict(final_model$gam, 
                            newdata = data.frame(Ano = anos_pred, Area_Antrop = "1"), 
                            se.fit = TRUE)
pred_antrop$fit <- pred_antrop_list$fit
pred_antrop$se.fit <- pred_antrop_list$se.fit
pred_antrop$lower <- pred_antrop$fit - 1.96 * pred_antrop$se.fit
pred_antrop$upper <- pred_antrop$fit + 1.96 * pred_antrop$se.fit

# Predictions for reference site
anos_ref <- seq(min(reference_abundance_data$Ano), max(reference_abundance_data$Ano), length.out = 100)
pred_ref <- data.frame(Ano = anos_ref)

pred_ref_list <- predict(ref_gam, newdata = data.frame(Ano = anos_ref), se.fit = TRUE)
pred_ref$fit <- pred_ref_list$fit
pred_ref$se.fit <- pred_ref_list$se.fit
pred_ref$lower <- pred_ref$fit - 1.96 * pred_ref$se.fit
pred_ref$upper <- pred_ref$fit + 1.96 * pred_ref$se.fit

# Annual means for anthropogenic sites
annual_means <- antrop_abundance_data %>%
  group_by(Ano) %>%
  summarise(
    mean = mean(Abundancia),
    se = sd(Abundancia)/sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

# Add group labels
pred_antrop$Grupo <- "Abandoned pastures"
pred_ref$Grupo <- "Natural grassland"
annual_means$Grupo <- "Abandoned pastures"
reference_abundance_data$Grupo <- "Natural grassland"

# Build figure
fig_main <- ggplot() +
  # Ribbons
  geom_ribbon(data = pred_antrop,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.22) +
  geom_ribbon(data = pred_ref,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.18) +
  
  # Lines
  geom_line(data = pred_antrop,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 2.0, lineend = "round") +
  geom_line(data = pred_ref,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 1.4) +
  
  # Points and error bars
  geom_errorbar(data = annual_means,
                aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.35, alpha = 0.7, linewidth = 0.8, show.legend = FALSE) +
  geom_point(data = annual_means,
             aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 3.6, stroke = 1.2) +
  geom_point(data = reference_abundance_data,
             aes(x = Ano, y = Abundancia, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.8, alpha = 0.7) +
  
  # Scales
  scale_color_manual(values = c("Abandoned pastures" = "#D55E00",
                                "Natural grassland" = "#009E73")) +
  scale_fill_manual(values = c("Abandoned pastures" = "#D55E00",
                               "Natural grassland" = "#009E73")) +
  scale_linetype_manual(values = c("Abandoned pastures" = "solid",
                                   "Natural grassland" = "23")) +
  scale_shape_manual(values = c("Abandoned pastures" = 21,
                                "Natural grassland" = 23)) +
  
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    fill = guide_legend(override.aes = list(alpha = 0.4)),
    linetype = guide_legend(),
    shape = guide_legend()
  ) +
  
  # Axes
  scale_x_continuous(
    breaks = seq(floor(primeiro_ano/2)*2, 2025, by = 2),
    limits = c(primeiro_ano, 2025),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  
  # Labels and theme
  labs(x = "Year", 
       y = "Abundance (individuals)", 
       color = NULL, fill = NULL, linetype = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "bottom",
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.height = unit(0.9, "lines"),
    legend.key.width = unit(1.6, "lines")
  )

ggsave(paste0(output_path, "Figure_2_abundance.tiff"), fig_main,
       width = 14, height = 10, units = "cm", dpi = 300, compression = "lzw")
cat("Main figure saved:", paste0(output_path, "Figure_2_abundance.tiff"), "\n")

# ==============================================================================
# 6. DIAGNOSTIC FIGURE
# ==============================================================================

cat("\nGENERATING DIAGNOSTIC FIGURE\n")
cat("-----------------------------\n")

residuos <- resid(final_model$lme, type = "normalized")
valores_ajustados <- fitted(final_model$lme)

tiff(paste0(output_path, "Figure_S1_abundance_diagnostics.tiff"),
     width = 2400, height = 1800, res = 300, compression = "lzw")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2), cex.lab = 1.1, cex.axis = 1)

# A. Residuals vs Fitted
plot(valores_ajustados, residuos,
     xlab = "Fitted values", ylab = "Standardized residuals",
     main = "A. Residuals vs Fitted",
     pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
abline(h = 0, col = "red", lty = 2, lwd = 2)
lines(lowess(valores_ajustados, residuos), col = "blue", lwd = 2)

# B. Q-Q Plot
qqnorm(residuos, main = "B. Normal Q-Q Plot",
       pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
qqline(residuos, col = "red", lwd = 2)

# C. Histogram
hist(residuos, breaks = 15, probability = TRUE,
     main = "C. Distribution of Residuals",
     xlab = "Standardized residuals",
     col = "lightgray", border = "darkgray")
curve(dnorm(x, mean = mean(residuos), sd = sd(residuos)),
      col = "red", lwd = 2, add = TRUE)

# D. ACF
acf(residuos, main = "D. Autocorrelation Function", lwd = 2)

dev.off()
cat("Diagnostic figure saved:", paste0(output_path, "Figure_S1_abundance_diagnostics.tiff"), "\n")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("================================================================================\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
