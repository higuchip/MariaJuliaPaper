# ================================================================================
# TEMPORAL DYNAMICS OF SPECIES RICHNESS IN ANTHROPOGENIC GRASSLANDS
# Complete analysis with figures and diagnostics
# ================================================================================

set.seed(123)

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(nlme)

# Create output directories
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)

# ================================================================================
# 1. DATA PREPARATION
# ================================================================================

cat("================================================================================\n")
cat("RICHNESS ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

# Identify height columns
height_cols <- grep("^h\\d{4}$", names(df_processed), value = TRUE)

# Convert to long format
richness_long <- df_processed %>%
  select(area, parc, ni, especies, all_of(height_cols)) %>%
  pivot_longer(
    cols = all_of(height_cols),
    names_to = "Ano",
    values_to = "Altura",
    names_prefix = "h"
  ) %>%
  mutate(Ano = as.numeric(Ano)) %>%
  filter(!is.na(Altura))

# Calculate richness (unique species per area and year)
richness_data <- richness_long %>%
  group_by(Area = as.factor(area), Ano) %>%
  summarise(Riqueza = n_distinct(especies), .groups = 'drop')

# Separate datasets
antrop_richness_data <- richness_data %>%
  filter(Area %in% c("1", "2", "4")) %>%
  mutate(Area_Antrop = factor(Area, levels = c("1", "2", "4")))

reference_richness_data <- richness_data %>%
  filter(Area == "3")

cat("DATA SUMMARY\n")
cat("------------\n")
cat("Anthropogenic sites: n =", nrow(antrop_richness_data), "\n")
cat("Reference site: n =", nrow(reference_richness_data), "\n")
cat("Period:", min(richness_data$Ano), "-", max(richness_data$Ano), "\n\n")

# ================================================================================
# 2. STATISTICAL MODELS
# ================================================================================

cat("MODEL FITTING\n")
cat("-------------\n")

k_antrop <- min(length(unique(antrop_richness_data$Ano)) - 1, 4)

# Initial model
antrop_gamm_initial <- gamm(Riqueza ~ s(Ano, k = k_antrop),
                            random = list(Area_Antrop = ~1),
                            data = antrop_richness_data,
                            method = "REML")

# Test heteroscedasticity
residuos_init <- resid(antrop_gamm_initial$lme, type = "normalized")
valores_ajustados_init <- fitted(antrop_gamm_initial$lme)
bp_model_init <- lm(residuos_init^2 ~ valores_ajustados_init)
bp_stat_init <- length(residuos_init) * summary(bp_model_init)$r.squared
bp_pval_init <- pchisq(bp_stat_init, df = 1, lower.tail = FALSE)

cat("Initial model - Breusch-Pagan test: p =", round(bp_pval_init, 4), "\n")

# Fit corrected model if needed
heteroscedasticity_corrected <- FALSE
if(bp_pval_init < 0.05) {
  cat("Heteroscedasticity detected - fitting corrected model...\n")
  antrop_gamm_corrected <- gamm(Riqueza ~ s(Ano, k = k_antrop),
                                random = list(Area_Antrop = ~1),
                                weights = varIdent(form = ~ 1 | Area_Antrop),
                                data = antrop_richness_data,
                                method = "REML")
  
  # Test corrected model
  residuos_corr <- resid(antrop_gamm_corrected$lme, type = "normalized")
  valores_ajustados_corr <- fitted(antrop_gamm_corrected$lme)
  bp_model_corr <- lm(residuos_corr^2 ~ valores_ajustados_corr)
  bp_stat_corr <- length(residuos_corr) * summary(bp_model_corr)$r.squared
  bp_pval_corr <- pchisq(bp_stat_corr, df = 1, lower.tail = FALSE)
  
  cat("Corrected model - Breusch-Pagan test: p =", round(bp_pval_corr, 4), "\n")
  
  # Compare models
  aic_init <- AIC(antrop_gamm_initial$lme)
  aic_corr <- AIC(antrop_gamm_corrected$lme)
  
  if(aic_corr < aic_init && bp_pval_corr > bp_pval_init) {
    antrop_gamm_final <- antrop_gamm_corrected
    heteroscedasticity_corrected <- TRUE
    final_bp_pval <- bp_pval_corr
    het_pval <- bp_pval_corr
    cat("Using corrected model (improved AIC and residuals)\n")
  } else {
    antrop_gamm_final <- antrop_gamm_initial
    final_bp_pval <- bp_pval_init
    het_pval <- bp_pval_init
    cat("Keeping initial model\n")
  }
} else {
  antrop_gamm_final <- antrop_gamm_initial
  final_bp_pval <- bp_pval_init
  het_pval <- bp_pval_init
  cat("No heteroscedasticity detected\n")
}

# Reference model
ref_gam <- gam(Riqueza ~ s(Ano, k = 3),
               data = reference_richness_data,
               method = "REML")

# Linear models for rates
lm_antrop <- lm(Riqueza ~ Ano, data = antrop_richness_data)
lm_ref <- lm(Riqueza ~ Ano, data = reference_richness_data)

# ================================================================================
# 3. EXTRACT RESULTS
# ================================================================================

cat("\nEXTRACTING RESULTS\n")
cat("------------------\n")

# Final model statistics
final_summary <- summary(antrop_gamm_final$gam)
final_smooth <- final_summary$s.table["s(Ano)",]
residuos_final <- resid(antrop_gamm_final$lme, type = "normalized")

# Rates of change
slope_antrop <- coef(lm_antrop)[2]
se_slope_antrop <- summary(lm_antrop)$coefficients[2, 2]
slope_ref <- coef(lm_ref)[2]

# Calculate increment
primeiro_ano <- min(antrop_richness_data$Ano)
ultimo_ano <- max(antrop_richness_data$Ano)
pred_inicial <- predict(antrop_gamm_final$gam, 
                        newdata = data.frame(Ano = primeiro_ano, Area_Antrop = "1"))
pred_final <- predict(antrop_gamm_final$gam, 
                      newdata = data.frame(Ano = ultimo_ano, Area_Antrop = "1"))
incremento_pct <- ((pred_final - pred_inicial) / pred_inicial) * 100

# Diagnostic tests
shapiro_test <- shapiro.test(residuos_final)
acf_result <- acf(residuos_final, plot = FALSE)
lags_sig <- which(abs(acf_result$acf[-1]) > qnorm(0.975)/sqrt(length(residuos_final)))

# ================================================================================
# 4. RESULTS TABLE
# ================================================================================

cat("\n=== MAIN RESULTS ===\n")

results_richness <- data.frame(
  Parameter = c(
    "Temporal trend (F)",
    "p-value",
    "EDF",
    "Trend type",
    "R² adjusted",
    "Linear rate (species/year)",
    "Total increment (%)",
    "Period analyzed",
    "N observations",
    "Heteroscedasticity corrected",
    "Shapiro-Wilk (W)",
    "Normality (p-value)",
    "Temporal autocorrelation",
    "Homoscedasticity (p-value)",
    "Reference rate (species/year)"
  ),
  
  Value = c(
    round(final_smooth["F"], 2),
    formatC(final_smooth["p-value"], format = "e", digits = 2),
    round(final_smooth["edf"], 2),
    ifelse(final_smooth["edf"] > 1.5, "Non-linear", "Linear"),
    round(final_summary$r.sq, 3),
    paste0(round(slope_antrop, 2), " ± ", round(se_slope_antrop, 2)),
    round(incremento_pct, 1),
    paste(primeiro_ano, "-", ultimo_ano),
    nrow(antrop_richness_data),
    ifelse(heteroscedasticity_corrected, "Yes", "No"),
    round(shapiro_test$statistic, 3),
    round(shapiro_test$p.value, 3),
    ifelse(length(lags_sig) == 0, "Not detected", "Present"),
    round(het_pval, 3),
    round(slope_ref, 2)
  ),
  
  Interpretation = c(
    ifelse(final_smooth["p-value"] < 0.05, "Significant", "Not significant"),
    ifelse(final_smooth["p-value"] < 0.01, "Highly significant", "Significant"),
    ifelse(final_smooth["edf"] > 2, "Strong non-linearity", "Slight non-linearity"),
    "Pattern of change",
    ifelse(final_summary$r.sq > 0.7, "Excellent fit", "Good fit"),
    ifelse(slope_antrop > 0, "Increasing diversity", "Decreasing diversity"),
    ifelse(incremento_pct > 0, "Net gain", "Net loss"),
    paste0(ultimo_ano - primeiro_ano + 1, " years"),
    "3 anthropogenic areas",
    ifelse(heteroscedasticity_corrected, "Corrected model", "Standard model"),
    ifelse(shapiro_test$p.value > 0.05, "Normal", "Non-normal"),
    ifelse(shapiro_test$p.value > 0.05, "Assumption met", "Check"),
    ifelse(length(lags_sig) == 0, "Assumption met", "Check"),
    ifelse(het_pval > 0.05, "Assumption met", "Residual issue"),
    "Natural grassland pattern"
  )
)

print(results_richness, row.names = FALSE)
write.csv(results_richness, "outputs/tables/Richness_results.csv", row.names = FALSE)

# ================================================================================
# 5. MAIN FIGURE
# ================================================================================

cat("\nGENERATING MAIN FIGURE\n")
cat("----------------------\n")

# Prepare predictions
anos_pred <- seq(min(antrop_richness_data$Ano), 
                 max(antrop_richness_data$Ano), length.out = 200)

# Anthropogenic predictions
pred_antrop <- data.frame(Ano = anos_pred, Area_Antrop = "1")
pred_vals <- predict(antrop_gamm_final$gam, newdata = pred_antrop, se.fit = TRUE)
pred_antrop$fit <- pred_vals$fit
pred_antrop$lower <- pred_antrop$fit - 1.96 * pred_vals$se.fit
pred_antrop$upper <- pred_antrop$fit + 1.96 * pred_vals$se.fit
pred_antrop$Grupo <- "Anthropogenic grasslands"

# Reference predictions
pred_ref <- data.frame(Ano = anos_pred)
pred_vals_ref <- predict(ref_gam, newdata = pred_ref, se.fit = TRUE)
pred_ref$fit <- pred_vals_ref$fit
pred_ref$lower <- pred_ref$fit - 1.96 * pred_vals_ref$se.fit
pred_ref$upper <- pred_ref$fit + 1.96 * pred_vals_ref$se.fit
pred_ref$Grupo <- "Natural grassland"

# Annual means for anthropogenic
annual_means <- antrop_richness_data %>%
  group_by(Ano) %>%
  summarise(mean = mean(Riqueza),
            se = sd(Riqueza)/sqrt(n()),
            .groups = 'drop') %>%
  mutate(Grupo = "Anthropogenic grasslands")

reference_richness_data$Grupo <- "Natural grassland"

# Create figure
fig_richness <- ggplot() +
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
            linewidth = 2.0) +
  geom_line(data = pred_ref,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 1.4) +
  
  # Points
  geom_errorbar(data = annual_means,
                aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.35, alpha = 0.7, linewidth = 0.8, show.legend = FALSE) +
  geom_point(data = annual_means,
             aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 3.6, stroke = 1.2) +
  geom_point(data = reference_richness_data,
             aes(x = Ano, y = Riqueza, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.8, alpha = 0.7) +
  
  # Scales
  scale_color_manual(values = c("Anthropogenic grasslands" = "#D55E00",
                                "Natural grassland" = "#009E73")) +
  scale_fill_manual(values = c("Anthropogenic grasslands" = "#D55E00",
                               "Natural grassland" = "#009E73")) +
  scale_linetype_manual(values = c("Anthropogenic grasslands" = "solid",
                                   "Natural grassland" = "dashed")) +
  scale_shape_manual(values = c("Anthropogenic grasslands" = 21,
                                "Natural grassland" = 23)) +
  
  scale_x_continuous(breaks = seq(floor(primeiro_ano/2)*2, ceiling(ultimo_ano/2)*2, by = 2)) +
  labs(x = "Year", y = "Species richness", 
       color = NULL, fill = NULL, linetype = NULL, shape = NULL) +
  
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    legend.text = element_text(size = 10, face = "bold")
  )

# Save
ggsave("outputs/figures/Figure_2_richness.tiff", fig_richness,
       width = 14, height = 10, units = "cm", dpi = 300, compression = "lzw")
cat("Main figure saved: Figure_2_richness.tiff\n")

# ================================================================================
# 6. DIAGNOSTIC FIGURE
# ================================================================================

cat("\nGENERATING DIAGNOSTIC FIGURE\n")

residuos <- resid(antrop_gamm_final$lme, type = "normalized")
valores_ajustados <- fitted(antrop_gamm_final$lme)

tiff("outputs/figures/Figure_S2_richness_diagnostics.tiff",
     width = 2400, height = 1800, res = 300, compression = "lzw")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2), cex.lab = 1.1, cex.axis = 1)

plot(valores_ajustados, residuos,
     xlab = "Fitted values", ylab = "Standardized residuals",
     main = "A. Residuals vs Fitted",
     pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
abline(h = 0, col = "red", lty = 2, lwd = 2)
lines(lowess(valores_ajustados, residuos), col = "blue", lwd = 2)

qqnorm(residuos, main = "B. Normal Q-Q Plot",
       pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
qqline(residuos, col = "red", lwd = 2)

hist(residuos, breaks = 15, probability = TRUE,
     main = "C. Distribution of Residuals",
     xlab = "Standardized residuals",
     col = "lightgray", border = "darkgray")
curve(dnorm(x, mean = mean(residuos), sd = sd(residuos)),
      col = "red", lwd = 2, add = TRUE)

acf(residuos, main = "D. Autocorrelation Function", lwd = 2)

dev.off()
cat("Diagnostic figure saved: Figure_S2_richness_diagnostics.tiff\n")

# ================================================================================
# 7. FINAL SUMMARY
# ================================================================================

cat("\n================================================================================\n")
cat("RICHNESS ANALYSIS COMPLETED\n")
cat("================================================================================\n")
cat("• Temporal trend:", ifelse(final_smooth["p-value"] < 0.05, "SIGNIFICANT", "not significant"), "\n")
cat("• Pattern:", ifelse(final_smooth["edf"] > 1.5, "NON-LINEAR", "LINEAR"), "\n")
cat("• Rate of change:", round(slope_antrop, 2), "species/year\n")
cat("• Total change:", round(incremento_pct, 1), "%\n")
cat("• Model diagnostics:", ifelse(shapiro_test$p.value > 0.05 & length(lags_sig) == 0 & het_pval > 0.05,
                                   "All assumptions met", "Check diagnostic plots"), "\n")
cat("================================================================================\n")
