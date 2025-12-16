# ==============================================================================
# TEMPORAL DYNAMICS OF MEAN HEIGHT
# ==============================================================================
#
# Description: Analyzes temporal changes in mean height of woody plants across
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
# Output:      ../artigo/biotropica/review/Figure_3_mean_height.tiff
#              ../artigo/biotropica/review/Figure_S2_height_diagnostics.tiff
#              outputs/tables/MeanHeight_Results_summary.csv
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
cat("MEAN HEIGHT ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

# Check dependency
if(!exists("df_processed")) {
  stop("Object 'df_processed' not found. Please run 01_data_processing.R first.")
}

# Identify height columns (hYYYY)
height_cols <- grep("^h\\d{4}$", names(df_processed), value = TRUE)

# Long format (individual-level heights)
height_long <- df_processed %>%
  select(area, parc, ni, especies, all_of(height_cols)) %>%
  pivot_longer(
    cols = all_of(height_cols),
    names_to = "Ano",
    values_to = "Altura",
    names_prefix = "h"
  ) %>%
  mutate(Ano = as.numeric(Ano)) %>%
  filter(!is.na(Altura))

# Mean height by Area x Year
mean_height_data <- height_long %>%
  group_by(Area = as.factor(area), Ano) %>%
  summarise(
    AlturaMedia   = mean(Altura, na.rm = TRUE),
    N_Individuos  = n(),
    .groups = "drop"
  )

# Separate datasets by ecosystem type
# Areas 1, 2, 4: Abandoned pastures (anthropogenic grasslands)
# Area 3: Natural grassland
antrop_height_data <- mean_height_data %>%
  filter(Area %in% c("1","2","4")) %>%
  mutate(Area_Antrop = factor(Area, levels = c("1","2","4")))

reference_height_data <- mean_height_data %>%
  filter(Area == "3")

primeiro_ano <- min(mean_height_data$Ano, na.rm = TRUE)
ultimo_ano_obs <- max(mean_height_data$Ano, na.rm = TRUE)

cat("DATA SUMMARY\n")
cat("------------\n")
cat("Anthropogenic sites (mean height): n =", nrow(antrop_height_data), "\n")
cat("Reference site (mean height): n =", nrow(reference_height_data), "\n")
cat("Period (observed):", primeiro_ano, "-", ultimo_ano_obs, "\n\n")

# ==============================================================================
# 2. STATISTICAL MODELS (GAMM + AIC Selection)
# ==============================================================================

cat("MODEL FITTING & SELECTION\n")
cat("-------------------------\n")

k_antrop_h <- max(3, min(4, length(unique(antrop_height_data$Ano)) - 1))

# --- A. Baseline Model (No Temporal Correlation) ---
cat("1. Fitting Baseline Model (No Temporal Correlation)...\n")
gamm_height_base <- gamm(
  AlturaMedia ~ s(Ano, k = k_antrop_h),
  random = list(Area_Antrop = ~ 1),
  data   = antrop_height_data,
  method = "REML"
)

# --- B. CAR1 Model (With Temporal Correlation) ---
cat("2. Fitting CAR1 Model (With Temporal Correlation)...\n")
cor_struct <- corCAR1(form = ~ Ano | Area_Antrop)

gamm_height_cor <- gamm(
  AlturaMedia ~ s(Ano, k = k_antrop_h),
  random = list(Area_Antrop = ~ 1),
  correlation = cor_struct,
  data   = antrop_height_data,
  method = "REML"
)

# --- C. Model Comparison (AIC) ---
aic_base <- AIC(gamm_height_base$lme)
aic_cor  <- AIC(gamm_height_cor$lme)
delta_aic <- aic_base - aic_cor

cat("\nMODEL COMPARISON (AIC):\n")
cat("Baseline AIC:", round(aic_base, 3), "\n")
cat("CAR1 AIC:    ", round(aic_cor, 3), "\n")
cat("Delta AIC:   ", round(delta_aic, 3), "\n")

# Select best model
if(delta_aic > 2) {
  cat("-> CONCLUSION: CAR1 model justified.\n")
  best_model_h <- gamm_height_cor
  used_cor <- TRUE
} else {
  cat("-> CONCLUSION: CAR1 not significantly better. Using baseline.\n")
  best_model_h <- gamm_height_base
  used_cor <- FALSE
}

# --- D. Heteroscedasticity Check ---
res_temp_h <- resid(best_model_h$lme, type = "normalized")
fit_temp_h <- fitted(best_model_h$lme)
bp_model_h <- lm(res_temp_h^2 ~ fit_temp_h)
bp_stat_h  <- length(res_temp_h) * summary(bp_model_h)$r.squared
bp_pval_h  <- pchisq(bp_stat_h, df = 1, lower.tail = FALSE)

cat("\nHeteroscedasticity check (p-value):", round(bp_pval_h, 4), "\n")

final_model_h <- best_model_h
heterosced_corrected_h <- FALSE

if (bp_pval_h < 0.05) {
  cat("Heteroscedasticity detected. Adding varIdent...\n")
  
  if(used_cor) {
    final_model_h <- gamm(
      AlturaMedia ~ s(Ano, k = k_antrop_h),
      random = list(Area_Antrop = ~ 1),
      correlation = cor_struct,
      weights = varIdent(form = ~ 1 | Area_Antrop),
      data   = antrop_height_data,
      method = "REML"
    )
  } else {
    final_model_h <- gamm(
      AlturaMedia ~ s(Ano, k = k_antrop_h),
      random = list(Area_Antrop = ~ 1),
      weights = varIdent(form = ~ 1 | Area_Antrop),
      data   = antrop_height_data,
      method = "REML"
    )
  }
  heterosced_corrected_h <- TRUE
  cat("Model updated with Variance Structure.\n")
}

# Reference model (simple GAM for natural grassland)
k_ref_h <- max(3, min(4, length(unique(reference_height_data$Ano)) - 1))
gam_height_ref <- gam(
  AlturaMedia ~ s(Ano, k = k_ref_h),
  data = reference_height_data,
  method = "REML"
)

# Linear models for rate reporting
lm_height_antrop <- lm(AlturaMedia ~ Ano, data = antrop_height_data)
lm_height_ref    <- lm(AlturaMedia ~ Ano, data = reference_height_data)

# ==============================================================================
# 3. EXTRACT RESULTS
# ==============================================================================

cat("\nEXTRACTING RESULTS (MEAN HEIGHT)\n")
cat("--------------------------------\n")

final_summary_h <- summary(final_model_h$gam)
final_smooth_h  <- final_summary_h$s.table["s(Ano)",]
residuos_final_h <- resid(final_model_h$lme, type = "normalized")

# Phi parameter (temporal autocorrelation)
phi_value <- NA
if(used_cor) {
  cs <- final_model_h$lme$modelStruct$corStruct
  phi_value <- as.numeric(coef(cs, unconstrained = FALSE))
  cat("Estimated temporal autocorrelation (Phi):", round(phi_value, 3), "\n")
} else {
  phi_value <- 0
  cat("No correlation structure used (Phi = 0)\n")
}

# Rates & Increment
slope_antrop_h    <- coef(lm_height_antrop)[2]
se_slope_antrop_h <- summary(lm_height_antrop)$coefficients[2,2]
slope_ref_h       <- coef(lm_height_ref)[2]
se_slope_ref_h    <- summary(lm_height_ref)$coefficients[2,2]

pred_ini_h <- predict(final_model_h$gam,
                      newdata = data.frame(Ano = primeiro_ano, Area_Antrop = "1"))
pred_fin_h <- predict(final_model_h$gam,
                      newdata = data.frame(Ano = ultimo_ano_obs, Area_Antrop = "1"))
incremento_pct_h <- ((pred_fin_h - pred_ini_h) / pred_ini_h) * 100

# Diagnostics
shapiro_test_h <- shapiro.test(residuos_final_h)

fit_final_h <- fitted(final_model_h$lme)
bp_final_model <- lm(residuos_final_h^2 ~ fit_final_h)
bp_pval_final_h <- pchisq(length(residuos_final_h) * summary(bp_final_model)$r.squared, df=1, lower.tail=FALSE)

acf_result_h <- acf(residuos_final_h, plot = FALSE)
lags_sig_h <- which(abs(acf_result_h$acf[-1]) > qnorm(0.975)/sqrt(length(residuos_final_h)))
status_acf <- ifelse(length(lags_sig_h) == 0, "OK", paste(length(lags_sig_h), "lags sig."))

# ==============================================================================
# 4. RESULTS TABLE
# ==============================================================================

cat("\n=== DETAILED RESULTS (MEAN HEIGHT) ===\n")

results_height_simple <- data.frame(
  Parameter = c(
    "Temporal trend (F)",
    "p-value",
    "Trend type (EDF)",
    "Adjusted R²",
    "Rate of change (m/year)",
    "Total increment (%)",
    "Baseline Model (AIC)",
    "CAR1 Model (AIC)",
    "Delta AIC",
    "Phi (Autocorrelation)",
    "Heteroscedasticity corrected",
    "Normality (Shapiro p)",
    "Homoscedasticity (BP p)",
    "Residual autocorrelation status",
    "Reference rate (m/year)"
  ),
  
  Value = c(
    round(final_smooth_h["F"], 2),
    formatC(final_smooth_h["p-value"], format = "e", digits = 2),
    paste0(ifelse(final_smooth_h["edf"] > 1.5, "Non-linear", "Linear"), " (", round(final_smooth_h["edf"], 2), ")"),
    round(final_summary_h$r.sq, 3),
    paste0(round(slope_antrop_h, 3), " ± ", round(se_slope_antrop_h, 3)),
    round(incremento_pct_h, 1),
    round(aic_base, 1),
    round(aic_cor, 1),
    round(delta_aic, 2),
    round(phi_value, 3),
    ifelse(heterosced_corrected_h, "Yes", "No"),
    round(shapiro_test_h$p.value, 3),
    round(bp_pval_final_h, 3),
    status_acf,
    paste0(round(slope_ref_h, 3), " ± ", round(se_slope_ref_h, 3))
  ),
  
  Interpretation = c(
    ifelse(final_smooth_h["p-value"] < 0.05, "Significant", "NS"),
    ifelse(final_smooth_h["p-value"] < 0.01, "**", "*"),
    "Curve complexity",
    "Model fit",
    "Vertical growth",
    "Accumulated change",
    "Without temporal correction",
    "With temporal correction",
    ifelse(delta_aic > 2, "CAR1 Justified", "CAR1 Indifferent"),
    ifelse(used_cor, "Dependency strength", "Not applied"),
    "Variance weights",
    ifelse(shapiro_test_h$p.value > 0.05, "Normal", "Deviation"),
    ifelse(bp_pval_final_h > 0.05, "Homogeneous", "Heterogeneous"),
    ifelse(length(lags_sig_h) == 0, "Success", "Alert"),
    "Reference"
  ),
  check.names = FALSE
)

print(results_height_simple, row.names = FALSE)
write.csv(results_height_simple, "outputs/tables/MeanHeight_Results_summary.csv", row.names = FALSE)

# ==============================================================================
# 5. MAIN FIGURE
# ==============================================================================

cat("\nGENERATING MAIN FIGURE (MEAN HEIGHT)\n")
cat("------------------------------------\n")

# Prediction grid
anos_pred_h <- seq(primeiro_ano, ultimo_ano_obs, length.out = 200)

# Anthropogenic predictions
pred_antrop_h <- data.frame(Ano = anos_pred_h, Area_Antrop = "1")
p_ant <- predict(final_model_h$gam, newdata = pred_antrop_h, se.fit = TRUE)
pred_antrop_h$fit   <- p_ant$fit
pred_antrop_h$lower <- p_ant$fit - 1.96 * p_ant$se.fit
pred_antrop_h$upper <- p_ant$fit + 1.96 * p_ant$se.fit

# Reference predictions
pred_ref_h <- data.frame(Ano = anos_pred_h)
p_ref <- predict(gam_height_ref, newdata = pred_ref_h, se.fit = TRUE)
pred_ref_h$fit   <- p_ref$fit
pred_ref_h$lower <- p_ref$fit - 1.96 * p_ref$se.fit
pred_ref_h$upper <- p_ref$fit + 1.96 * p_ref$se.fit

# Annual means
annual_means_h <- antrop_height_data %>%
  group_by(Ano) %>%
  summarise(
    mean = mean(AlturaMedia, na.rm = TRUE),
    se   = sd(AlturaMedia,  na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Add group labels
pred_antrop_h$Grupo         <- "Abandoned pastures"
pred_ref_h$Grupo            <- "Natural grassland"
annual_means_h$Grupo        <- "Abandoned pastures"
reference_height_data$Grupo <- "Natural grassland"

# Build figure
fig_height <- ggplot() +
  # Ribbons
  geom_ribbon(data = pred_antrop_h,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.22) +
  geom_ribbon(data = pred_ref_h,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.18) +
  
  # Lines
  geom_line(data = pred_antrop_h,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 2.0, lineend = "round") +
  geom_line(data = pred_ref_h,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 1.4) +
  
  # Points and error bars
  geom_errorbar(data = annual_means_h,
                aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.35, alpha = 0.7, linewidth = 0.8, show.legend = FALSE) +
  geom_point(data = annual_means_h,
             aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 3.6, stroke = 1.2) +
  geom_point(data = reference_height_data,
             aes(x = Ano, y = AlturaMedia, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.8, alpha = 0.7) +
  
  # Scales
  scale_color_manual(values = c("Abandoned pastures" = "#D55E00",
                                "Natural grassland"  = "#009E73")) +
  scale_fill_manual(values  = c("Abandoned pastures" = "#D55E00",
                                "Natural grassland"  = "#009E73")) +
  scale_linetype_manual(values = c("Abandoned pastures" = "solid",
                                   "Natural grassland"  = "23")) +
  scale_shape_manual(values = c("Abandoned pastures" = 21,
                                "Natural grassland"  = 23)) +
  
  guides(
    color    = guide_legend(override.aes = list(size = 3)),
    fill     = guide_legend(override.aes = list(alpha = 0.4)),
    linetype = guide_legend(),
    shape    = guide_legend()
  ) +
  
  # Axes
  scale_x_continuous(
    breaks = seq(floor(primeiro_ano/2)*2, 2025, by = 2),
    limits = c(primeiro_ano, 2025),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  
  # Labels and theme
  labs(x = "Year", y = "Mean height (m)",
       color = NULL, fill = NULL, linetype = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    axis.line  = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "bottom",
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.height = unit(0.9, "lines"),
    legend.key.width  = unit(1.6, "lines")
  )

ggsave(paste0(output_path, "Figure_3_mean_height.tiff"), fig_height,
       width = 14, height = 10, units = "cm", dpi = 300, compression = "lzw")
cat("Main figure saved:", paste0(output_path, "Figure_3_mean_height.tiff"), "\n")

# ==============================================================================
# 6. DIAGNOSTIC FIGURE
# ==============================================================================

cat("\nGENERATING DIAGNOSTIC FIGURE (MEAN HEIGHT)\n")
cat("------------------------------------------\n")

residuos_h <- resid(final_model_h$lme, type = "normalized")
valores_ajustados_h <- fitted(final_model_h$lme)

tiff(paste0(output_path, "Figure_S2_height_diagnostics.tiff"),
     width = 2400, height = 1800, res = 300, compression = "lzw")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2), cex.lab = 1.1, cex.axis = 1)

# A. Residuals vs Fitted
plot(valores_ajustados_h, residuos_h,
     xlab = "Fitted values", ylab = "Standardized residuals",
     main = "A. Residuals vs Fitted", pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
abline(h = 0, col = "red", lty = 2, lwd = 2)
lines(lowess(valores_ajustados_h, residuos_h), col = "blue", lwd = 2)

# B. Q-Q Plot
qqnorm(residuos_h, main = "B. Normal Q-Q Plot",
       pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
qqline(residuos_h, col = "red", lwd = 2)

# C. Histogram
hist(residuos_h, breaks = 15, probability = TRUE,
     main = "C. Distribution of Residuals",
     xlab = "Standardized residuals", col = "lightgray", border = "darkgray")
curve(dnorm(x, mean = mean(residuos_h), sd = sd(residuos_h)),
      col = "red", lwd = 2, add = TRUE)

# D. ACF
acf(residuos_h, main = "D. Autocorrelation Function", lwd = 2)

dev.off()
cat("Diagnostic figure saved:", paste0(output_path, "Figure_S2_height_diagnostics.tiff"), "\n")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("================================================================================\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
