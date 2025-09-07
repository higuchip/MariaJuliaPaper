# ================================================================================
# TEMPORAL DYNAMICS OF MEAN HEIGHT IN ANTHROPOGENIC GRASSLANDS
# Complete analysis with comprehensive results table
# ================================================================================

set.seed(123)

# Packages -----------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(nlme)
})

# Output dirs -------------------------------------------------------------------
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)

# ================================================================================
# 1. DATA PREPARATION
# ================================================================================

cat("================================================================================\n")
cat("MEAN HEIGHT ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

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

# Mean height by Area x Ano (this is what we model/plot)
mean_height_data <- height_long %>%
  group_by(Area = as.factor(area), Ano) %>%
  summarise(
    AlturaMedia   = mean(Altura, na.rm = TRUE),
    N_Individuos  = n(),
    .groups = "drop"
  )

# Split anthropogenic (áreas 1,2,4) vs reference (área 3)
antrop_height_data <- mean_height_data %>%
  filter(Area %in% c("1","2","4")) %>%
  mutate(Area_Antrop = factor(Area, levels = c("1","2","4")))

reference_height_data <- mean_height_data %>%
  filter(Area == "3")

primeiro_ano <- min(mean_height_data$Ano, na.rm = TRUE)
ultimo_ano_obs <- max(mean_height_data$Ano, na.rm = TRUE)
ultimo_ano_plot <- 2025  # força o plot até 2025

cat("DATA SUMMARY\n")
cat("------------\n")
cat("Anthropogenic sites (mean height): n =", nrow(antrop_height_data), "\n")
cat("Reference site (mean height): n =", nrow(reference_height_data), "\n")
cat("Period (observed):", primeiro_ano, "-", ultimo_ano_obs, "\n\n")

# ================================================================================
# 2. STATISTICAL MODELS 
#    GAMM (global smooth in time) with random intercept for Area; test heteroscedasticity
# ================================================================================

cat("MODEL FITTING (MEAN HEIGHT)\n")
cat("---------------------------\n")

k_antrop_h <- max(3, min(4, length(unique(antrop_height_data$Ano)) - 1))

# Initial model (homogeneous variance)
gamm_height_initial <- gamm(
  AlturaMedia ~ s(Ano, k = k_antrop_h),
  random = list(Area_Antrop = ~ 1),
  data   = antrop_height_data,
  method = "REML"
)

# Breusch-Pagan style check (proxy) on squared normalized residuals ~ fitted
res_init_h   <- resid(gamm_height_initial$lme, type = "normalized")
fit_init_h   <- fitted(gamm_height_initial$lme)
bp_model_h   <- lm(res_init_h^2 ~ fit_init_h)
bp_stat_h    <- length(res_init_h) * summary(bp_model_h)$r.squared
bp_pval_h    <- pchisq(bp_stat_h, df = 1, lower.tail = FALSE)

cat("Initial model - Breusch-Pagan test (mean height): p =", round(bp_pval_h, 4), "\n")

# Fit corrected model if needed (variance by area)
heteroscedasticity_corrected_h <- FALSE
if (bp_pval_h < 0.05) {
  cat("Heteroscedasticity detected - fitting corrected model (varIdent by Area)...\n")
  gamm_height_corrected <- gamm(
    AlturaMedia ~ s(Ano, k = k_antrop_h),
    random  = list(Area_Antrop = ~ 1),
    weights = varIdent(form = ~ 1 | Area_Antrop),
    data    = antrop_height_data,
    method  = "REML"
  )
  
  # Re-test
  res_corr_h <- resid(gamm_height_corrected$lme, type = "normalized")
  fit_corr_h <- fitted(gamm_height_corrected$lme)
  bp_model_corr_h <- lm(res_corr_h^2 ~ fit_corr_h)
  bp_stat_corr_h  <- length(res_corr_h) * summary(bp_model_corr_h)$r.squared
  bp_pval_corr_h  <- pchisq(bp_stat_corr_h, df = 1, lower.tail = FALSE)
  
  cat("Corrected model - Breusch-Pagan test: p =", round(bp_pval_corr_h, 4), "\n")
  
  # Compare AIC and residuals
  aic_init_h <- AIC(gamm_height_initial$lme)
  aic_corr_h <- AIC(gamm_height_corrected$lme)
  
  if (aic_corr_h < aic_init_h && bp_pval_corr_h > bp_pval_h) {
    gamm_height_final <- gamm_height_corrected
    heteroscedasticity_corrected_h <- TRUE
    final_bp_pval_h <- bp_pval_corr_h
    cat("Using corrected model (improved AIC and residuals)\n")
  } else {
    gamm_height_final <- gamm_height_initial
    final_bp_pval_h <- bp_pval_h
    cat("Keeping initial model\n")
  }
} else {
  gamm_height_final <- gamm_height_initial
  final_bp_pval_h <- bp_pval_h
  cat("No heteroscedasticity detected\n")
}

# Reference model (single series)
k_ref_h <- max(3, min(4, length(unique(reference_height_data$Ano)) - 1))
gam_height_ref <- gam(
  AlturaMedia ~ s(Ano, k = k_ref_h),
  data = reference_height_data,
  method = "REML"
)

# Linear models for average rate (for reporting)
lm_height_antrop <- lm(AlturaMedia ~ Ano, data = antrop_height_data)
lm_height_ref    <- lm(AlturaMedia ~ Ano, data = reference_height_data)

# ================================================================================
# 3. EXTRACT RESULTS
# ================================================================================

cat("\nEXTRACTING RESULTS (MEAN HEIGHT)\n")
cat("--------------------------------\n")

final_summary_h <- summary(gamm_height_final$gam)
final_smooth_h  <- final_summary_h$s.table["s(Ano)",]
residuos_final_h <- resid(gamm_height_final$lme, type = "normalized")

# Rates (m/ano)
slope_antrop_h    <- coef(lm_height_antrop)[2]
se_slope_antrop_h <- summary(lm_height_antrop)$coefficients[2,2]
slope_ref_h       <- coef(lm_height_ref)[2]
se_slope_ref_h    <- summary(lm_height_ref)$coefficients[2,2]

# Increment (%) using smooth predictions at first/last observed year
pred_ini_h <- predict(gamm_height_final$gam,
                      newdata = data.frame(Ano = primeiro_ano, Area_Antrop = "1"))
pred_fin_h <- predict(gamm_height_final$gam,
                      newdata = data.frame(Ano = ultimo_ano_obs, Area_Antrop = "1"))
incremento_pct_h <- ((pred_fin_h - pred_ini_h) / pred_ini_h) * 100

# Diagnostics
shapiro_test_h <- shapiro.test(residuos_final_h)
acf_result_h   <- acf(residuos_final_h, plot = FALSE)
lags_sig_h     <- which(abs(acf_result_h$acf[-1]) > qnorm(0.975)/sqrt(length(residuos_final_h)))

# ================================================================================
# 4. RESULTS TABLE (simple)
# ================================================================================

cat("\n=== RESULTADOS PRINCIPAIS (ALTURA MÉDIA) ===\n")

results_height_simple <- data.frame(
  Parâmetro = c(
    "Tendência temporal (F)",
    "Valor p",
    "EDF (graus de liberdade efetivos)",
    "Tipo de tendência",
    "R² ajustado",
    "Taxa de mudança linear (m/ano)",
    "Incremento total (%)",
    "Período analisado",
    "N observações",
    "Heterocedasticidade corrigida",
    "Teste Shapiro-Wilk (W)",
    "Normalidade (p-valor)",
    "Autocorrelação temporal",
    "Homocedasticidade final (p-valor)",
    "Taxa referência (m/ano)"
  ),
  Valor = c(
    round(final_smooth_h["F"], 2),
    formatC(final_smooth_h["p-value"], format = "e", digits = 2),
    round(final_smooth_h["edf"], 2),
    ifelse(final_smooth_h["edf"] > 1.5, "Não-linear", "Linear"),
    round(final_summary_h$r.sq, 3),
    paste0(round(slope_antrop_h, 3), " ± ", round(se_slope_antrop_h, 3)),
    round(incremento_pct_h, 1),
    paste(primeiro_ano, "-", ultimo_ano_obs),
    nrow(antrop_height_data),
    ifelse(heteroscedasticity_corrected_h, "Sim", "Não"),
    round(shapiro_test_h$statistic, 3),
    round(shapiro_test_h$p.value, 3),
    ifelse(length(lags_sig_h) == 0, "Não detectada", "Presente"),
    round(final_bp_pval_h, 3),
    paste0(round(slope_ref_h, 3), " ± ", round(se_slope_ref_h, 3))
  ),
  Interpretação = c(
    ifelse(final_smooth_h["p-value"] < 0.05, "Significativo", "Não significativo"),
    ifelse(final_smooth_h["p-value"] < 0.01, "Altamente significativo", "Significativo"),
    ifelse(final_smooth_h["edf"] > 2, "Forte não-linearidade", "Leve não-linearidade"),
    "Padrão temporal ajustado por Área",
    ifelse(final_summary_h$r.sq > 0.7, "Excelente ajuste", "Bom ajuste"),
    "Variação média anual",
    "Mudança acumulada no período",
    paste0(ultimo_ano_obs - primeiro_ano + 1, " anos"),
    "3 áreas antropogênicas",
    ifelse(heteroscedasticity_corrected_h, "Modelo corrigido", "Modelo padrão"),
    ifelse(shapiro_test_h$p.value > 0.05, "Normal", "Não normal"),
    ifelse(shapiro_test_h$p.value > 0.05, "Pressuposto atendido", "Verificar"),
    ifelse(length(lags_sig_h) == 0, "Pressuposto atendido", "Verificar"),
    ifelse(final_bp_pval_h > 0.05, "Pressuposto atendido", "Problema residual"),
    "Referência próxima de estável"
  ),
  check.names = FALSE
)

print(results_height_simple, row.names = FALSE)
write.csv(results_height_simple, "outputs/tables/MeanHeight_Results_summary.csv", row.names = FALSE)

# ================================================================================
# 5. MAIN FIGURE (with LEGEND) - axis up to 2025
# ================================================================================

cat("\nGENERATING MAIN FIGURE (MEAN HEIGHT)\n")
cat("------------------------------------\n")

# Prediction grid (for plotting we extend to 2025 just for axis; predictions stay within observed range)
anos_pred_h <- seq(primeiro_ano, ultimo_ano_obs, length.out = 200)

# Anthropogenic predictions (need a dummy level for Area_Antrop)
pred_antrop_h <- data.frame(Ano = anos_pred_h, Area_Antrop = "1")
p_ant <- predict(gamm_height_final$gam, newdata = pred_antrop_h, se.fit = TRUE)
pred_antrop_h$fit   <- p_ant$fit
pred_antrop_h$lower <- p_ant$fit - 1.96 * p_ant$se.fit
pred_antrop_h$upper <- p_ant$fit + 1.96 * p_ant$se.fit

# Reference predictions
pred_ref_h <- data.frame(Ano = anos_pred_h)
p_ref <- predict(gam_height_ref, newdata = pred_ref_h, se.fit = TRUE)
pred_ref_h$fit   <- p_ref$fit
pred_ref_h$lower <- p_ref$fit - 1.96 * p_ref$se.fit
pred_ref_h$upper <- p_ref$fit + 1.96 * p_ref$se.fit

# Annual means across anthropogenic areas (to show points ± SE)
annual_means_h <- antrop_height_data %>%
  group_by(Ano) %>%
  summarise(
    mean = mean(AlturaMedia, na.rm = TRUE),
    se   = sd(AlturaMedia,  na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Mark groups for legend
pred_antrop_h$Grupo              <- "Anthropogenic grasslands"
pred_ref_h$Grupo                 <- "Natural grassland"
annual_means_h$Grupo             <- "Anthropogenic grasslands"
reference_height_data$Grupo      <- "Natural grassland"

fig_height <- ggplot() +
  # ribbons
  geom_ribbon(data = pred_antrop_h,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.22, show.legend = TRUE) +
  geom_ribbon(data = pred_ref_h,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.18, show.legend = TRUE) +
  # lines
  geom_line(data = pred_antrop_h,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 2.0, lineend = "round") +
  geom_line(data = pred_ref_h,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 1.4) +
  # points/bars anthropogenic (annual means)
  geom_errorbar(data = annual_means_h,
                aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.35, alpha = 0.7, linewidth = 0.8, show.legend = FALSE) +
  geom_point(data = annual_means_h,
             aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 3.6, stroke = 1.2) +
  # points reference (observed mean by area=3)
  geom_point(data = reference_height_data,
             aes(x = Ano, y = AlturaMedia, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.8, alpha = 0.7) +
  # scales & legend
  scale_color_manual(values = c("Anthropogenic grasslands" = "#D55E00",
                                "Natural grassland"        = "#009E73")) +
  scale_fill_manual(values  = c("Anthropogenic grasslands" = "#D55E00",
                                "Natural grassland"        = "#009E73")) +
  scale_linetype_manual(values = c("Anthropogenic grasslands" = "solid",
                                   "Natural grassland"        = "23")) +
  scale_shape_manual(values = c("Anthropogenic grasslands" = 21,
                                "Natural grassland"        = 23)) +
  guides(
    color    = guide_legend(override.aes = list(size = 3)),
    fill     = guide_legend(override.aes = list(alpha = 0.4)),
    linetype = guide_legend(),
    shape    = guide_legend()
  ) +
  scale_x_continuous(
    breaks = seq(floor(primeiro_ano/2)*2, 2025, by = 2),
    limits = c(primeiro_ano, 2025),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(x = "Year", y = "Mean height (m)",
       color = NULL, fill = NULL, linetype = NULL, shape = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text  = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    axis.line  = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.height = unit(0.9, "lines"),
    legend.key.width  = unit(1.6, "lines")
  )

ggsave("outputs/figures/Figure_2_mean_height.tiff", fig_height,
       width = 14, height = 10, units = "cm", dpi = 300, compression = "lzw")

cat("Main figure saved: outputs/figures/Figure_2_mean_height.tiff\n")

# ================================================================================
# 6. RESIDUAL DIAGNOSTICS FIGURE
# ================================================================================

cat("\nGENERATING DIAGNOSTIC FIGURE (MEAN HEIGHT)\n")
cat("-------------------------------------------\n")

residuos_h         <- resid(gamm_height_final$lme, type = "normalized")
valores_ajustados_h <- fitted(gamm_height_final$lme)

tiff("outputs/figures/Figure_S2_diagnostics_mean_height.tiff",
     width = 2400, height = 1800, res = 300, compression = "lzw")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2), cex.lab = 1.1, cex.axis = 1)

# A) Residuals vs Fitted
plot(valores_ajustados_h, residuos_h,
     xlab = "Fitted values", ylab = "Standardized residuals",
     main = "A. Residuals vs Fitted",
     pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
abline(h = 0, col = "red", lty = 2, lwd = 2)
lines(lowess(valores_ajustados_h, residuos_h), col = "blue", lwd = 2)

# B) Normal Q-Q
qqnorm(residuos_h, main = "B. Normal Q-Q Plot",
       pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
qqline(residuos_h, col = "red", lwd = 2)

# C) Histogram
hist(residuos_h, breaks = 15, probability = TRUE,
     main = "C. Distribution of Residuals",
     xlab = "Standardized residuals",
     col = "lightgray", border = "darkgray")
curve(dnorm(x, mean = mean(residuos_h), sd = sd(residuos_h)),
      col = "red", lwd = 2, add = TRUE)

# D) ACF
acf(residuos_h, main = "D. Autocorrelation Function", lwd = 2)

dev.off()
cat("Diagnostic figure saved: outputs/figures/Figure_S2_diagnostics_mean_height.tiff\n")

# ================================================================================
# 7. ASSUMPTIONS SUMMARY (console)
# ================================================================================

cat("\n================================================================================\n")
cat("ASSUMPTIONS SUMMARY (MEAN HEIGHT)\n")
cat("================================================================================\n")
cat("Normality:", ifelse(shapiro_test_h$p.value > 0.05, "✓ OK", "✗ Violated"), "\n")
cat("Homoscedasticity:", ifelse(final_bp_pval_h > 0.05, "✓ OK", "✗ Violated"), "\n")
cat("Independence (ACF):", ifelse(length(lags_sig_h) == 0, "✓ OK", "✗ Violated"), "\n")
cat("================================================================================\n")
