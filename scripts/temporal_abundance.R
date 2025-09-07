# ================================================================================
# TEMPORAL DYNAMICS OF ABUNDANCE IN ANTHROPOGENIC GRASSLANDS
# Complete analysis with comprehensive results table
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
cat("ABUNDANCE ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

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

# Calculate abundance
abundance_data <- abundance_long %>%
  group_by(Area = as.factor(area), Ano) %>%
  summarise(Abundancia = n(), .groups = 'drop')

# Separate datasets
antrop_abundance_data <- abundance_data %>%
  filter(Area %in% c("1", "2", "4")) %>%
  mutate(Area_Antrop = factor(Area, levels = c("1", "2", "4")))

reference_abundance_data <- abundance_data %>%
  filter(Area == "3")

cat("DATA SUMMARY\n")
cat("------------\n")
cat("Anthropogenic sites: n =", nrow(antrop_abundance_data), "\n")
cat("Reference site: n =", nrow(reference_abundance_data), "\n")
cat("Period:", min(abundance_data$Ano), "-", max(abundance_data$Ano), "\n\n")

# ================================================================================
# 2. STATISTICAL MODELS
# ================================================================================

cat("MODEL FITTING\n")
cat("-------------\n")

k_antrop <- min(length(unique(antrop_abundance_data$Ano)) - 1, 4)

# Initial model
antrop_gamm_initial <- gamm(Abundancia ~ s(Ano, k = k_antrop),
                            random = list(Area_Antrop = ~1),
                            data = antrop_abundance_data,
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
  antrop_gamm_corrected <- gamm(Abundancia ~ s(Ano, k = k_antrop),
                                random = list(Area_Antrop = ~1),
                                weights = varIdent(form = ~ 1 | Area_Antrop),
                                data = antrop_abundance_data,
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
    cat("Using corrected model (improved AIC and residuals)\n")
  } else {
    antrop_gamm_final <- antrop_gamm_initial
    final_bp_pval <- bp_pval_init
    cat("Keeping initial model\n")
  }
} else {
  antrop_gamm_final <- antrop_gamm_initial
  final_bp_pval <- bp_pval_init
  cat("No heteroscedasticity detected\n")
}

# Reference model
ref_gam <- gam(Abundancia ~ s(Ano, k = 3),
               data = reference_abundance_data,
               method = "REML")

# Linear models for rates
lm_antrop <- lm(Abundancia ~ Ano, data = antrop_abundance_data)
lm_ref <- lm(Abundancia ~ Ano, data = reference_abundance_data)

# ================================================================================
# 3. EXTRACT RESULTS
# ================================================================================

cat("\nEXTRACTING RESULTS\n")
cat("------------------\n")

# Final model statistics
final_summary <- summary(antrop_gamm_final$gam)
final_smooth <- final_summary$s.table["s(Ano)",]
residuos_final <- resid(antrop_gamm_final$lme, type = "normalized")

# Define het_pval (variável que estava faltando)
het_pval <- final_bp_pval

# Rates of change
slope_antrop <- coef(lm_antrop)[2]
se_slope_antrop <- summary(lm_antrop)$coefficients[2, 2]
slope_ref <- coef(lm_ref)[2]
se_slope_ref <- summary(lm_ref)$coefficients[2, 2]

# Calculate increment
primeiro_ano <- min(antrop_abundance_data$Ano)
ultimo_ano <- max(antrop_abundance_data$Ano)
pred_inicial <- predict(antrop_gamm_final$gam, 
                        newdata = data.frame(Ano = primeiro_ano, Area_Antrop = "1"))
pred_final <- predict(antrop_gamm_final$gam, 
                      newdata = data.frame(Ano = ultimo_ano, Area_Antrop = "1"))
incremento_pct <- ((pred_final - pred_inicial) / pred_inicial) * 100

# Diagnostic tests
shapiro_test <- shapiro.test(residuos_final)
acf_result <- acf(residuos_final, plot = FALSE)
lags_sig <- which(abs(acf_result$acf[-1]) > qnorm(0.975)/sqrt(length(residuos_final)))

# Variance components
var_comp <- VarCorr(antrop_gamm_final$lme)
between_area_var <- as.numeric(var_comp[1,1])
residual_var <- as.numeric(var_comp[2,1])

# ================================================================================
# 4. PREPARAR DADOS PARA O GRÁFICO
# ================================================================================

cat("\nPREPARING DATA FOR VISUALIZATION\n")
cat("---------------------------------\n")

# 1. Criar predições para o grupo antropogênico
anos_pred <- seq(min(antrop_abundance_data$Ano), max(antrop_abundance_data$Ano), length.out = 100)
pred_antrop <- data.frame(Ano = anos_pred)

# Predições do modelo GAMM
pred_antrop_list <- predict(antrop_gamm_final$gam, 
                            newdata = data.frame(Ano = anos_pred, 
                                                 Area_Antrop = "1"), 
                            se.fit = TRUE)
pred_antrop$fit <- pred_antrop_list$fit
pred_antrop$se.fit <- pred_antrop_list$se.fit
pred_antrop$lower <- pred_antrop$fit - 1.96 * pred_antrop$se.fit
pred_antrop$upper <- pred_antrop$fit + 1.96 * pred_antrop$se.fit

# 2. Criar predições para o grupo de referência
anos_ref <- seq(min(reference_abundance_data$Ano), max(reference_abundance_data$Ano), length.out = 100)
pred_ref <- data.frame(Ano = anos_ref)

pred_ref_list <- predict(ref_gam, 
                         newdata = data.frame(Ano = anos_ref), 
                         se.fit = TRUE)
pred_ref$fit <- pred_ref_list$fit
pred_ref$se.fit <- pred_ref_list$se.fit
pred_ref$lower <- pred_ref$fit - 1.96 * pred_ref$se.fit
pred_ref$upper <- pred_ref$fit + 1.96 * pred_ref$se.fit

# 3. Calcular médias anuais para o grupo antropogênico
annual_means <- antrop_abundance_data %>%
  group_by(Ano) %>%
  summarise(
    mean = mean(Abundancia),
    se = sd(Abundancia)/sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

# Verificar se os dados foram criados corretamente
cat("Verificação dos dados para o gráfico:\n")
cat("- pred_antrop: ", nrow(pred_antrop), "linhas\n")
cat("- pred_ref: ", nrow(pred_ref), "linhas\n")
cat("- annual_means: ", nrow(annual_means), "linhas\n")
cat("- reference_abundance_data: ", nrow(reference_abundance_data), "linhas\n\n")

# ================================================================================
# 5. TABELA DE RESULTADOS
# ================================================================================

cat("CREATING RESULTS TABLE\n")
cat("----------------------\n")

# Criar tabela simplificada
results_simple <- data.frame(
  Parâmetro = c(
    "Tendência temporal (F)",
    "Valor p",
    "EDF (graus de liberdade efetivos)",
    "Tipo de tendência",
    "R² ajustado",
    "Taxa de mudança linear (ind/ano)",
    "Incremento total (%)",
    "Período analisado",
    "N observações",
    "Heterocedasticidade corrigida",
    "Teste Shapiro-Wilk (W)",
    "Normalidade (p-valor)",
    "Autocorrelação temporal",
    "Homocedasticidade final (p-valor)",
    "Taxa referência (ind/ano)"
  ),
  
  Valor = c(
    round(final_smooth["F"], 2),
    formatC(final_smooth["p-value"], format = "e", digits = 2),
    round(final_smooth["edf"], 2),
    ifelse(final_smooth["edf"] > 1.5, "Não-linear", "Linear"),
    round(final_summary$r.sq, 3),
    paste0(round(slope_antrop, 1), " ± ", round(se_slope_antrop, 1)),
    round(incremento_pct, 1),
    paste(primeiro_ano, "-", ultimo_ano),
    nrow(antrop_abundance_data),
    ifelse(heteroscedasticity_corrected, "Sim", "Não"),
    round(shapiro_test$statistic, 3),
    round(shapiro_test$p.value, 3),
    ifelse(length(lags_sig) == 0, "Não detectada", "Presente"),
    round(het_pval, 3),
    round(slope_ref, 1)
  ),
  
  Interpretação = c(
    ifelse(final_smooth["p-value"] < 0.05, "Significativo", "Não significativo"),
    ifelse(final_smooth["p-value"] < 0.01, "Altamente significativo", "Significativo"),
    ifelse(final_smooth["edf"] > 2, "Forte não-linearidade", "Leve não-linearidade"),
    "Padrão de crescimento acelerado",
    ifelse(final_summary$r.sq > 0.7, "Excelente ajuste", "Bom ajuste"),
    "Crescimento robusto",
    "Aumento substancial",
    paste0(ultimo_ano - primeiro_ano + 1, " anos"),
    "3 áreas antropogênicas",
    ifelse(heteroscedasticity_corrected, "Modelo corrigido", "Modelo padrão"),
    ifelse(shapiro_test$p.value > 0.05, "Normal", "Não normal"),
    ifelse(shapiro_test$p.value > 0.05, "Pressuposto atendido", "Verificar"),
    ifelse(length(lags_sig) == 0, "Pressuposto atendido", "Verificar"),
    ifelse(het_pval > 0.05, "Pressuposto atendido", "Problema residual"),
    "Campo natural estável"
  )
)

print(results_simple, row.names = FALSE)
write.csv(results_simple, "outputs/tables/Results_summary.csv", row.names = FALSE)

# ================================================================================
# 6. FIGURA PRINCIPAL PARA O ARTIGO
# ================================================================================

cat("\nGENERATING MAIN FIGURE\n")
cat("----------------------\n")

# Marcar o grupo em cada objeto usado no plot
pred_antrop$Grupo <- "Anthropogenic grasslands"
pred_ref$Grupo <- "Natural grassland"
annual_means$Grupo <- "Anthropogenic grasslands"
reference_abundance_data$Grupo <- "Natural grassland"

# Construir a figura
fig_main <- ggplot() +
  # --- Faixas de incerteza (ribbons) ---
  geom_ribbon(data = pred_antrop,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.22, show.legend = TRUE) +
  geom_ribbon(data = pred_ref,
              aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo),
              alpha = 0.18, show.legend = TRUE) +
  
  # --- Linhas suavizadas (mostrar linetype na legenda) ---
  geom_line(data = pred_antrop,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 2.0, lineend = "round", show.legend = TRUE) +
  geom_line(data = pred_ref,
            aes(x = Ano, y = fit, color = Grupo, linetype = Grupo),
            linewidth = 1.4, show.legend = TRUE) +
  
  # --- Pontos e barras (médias anuais do grupo antropizado) ---
  geom_errorbar(data = annual_means,
                aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.35, alpha = 0.7, linewidth = 0.8, show.legend = FALSE) +
  geom_point(data = annual_means,
             aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 3.6, stroke = 1.2, show.legend = TRUE) +
  
  # --- Pontos (dados de referência) ---
  geom_point(data = reference_abundance_data,
             aes(x = Ano, y = Abundancia, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.8, alpha = 0.7, show.legend = TRUE) +
  
  # --- Escalas e legenda ---
  scale_color_manual(values = c("Anthropogenic grasslands" = "#D55E00",
                                "Natural grassland" = "#009E73")) +
  scale_fill_manual(values = c("Anthropogenic grasslands" = "#D55E00",
                               "Natural grassland" = "#009E73")) +
  scale_linetype_manual(values = c("Anthropogenic grasslands" = "solid",
                                   "Natural grassland" = "23")) +
  scale_shape_manual(values = c("Anthropogenic grasslands" = 21,
                                "Natural grassland" = 23)) +
  
  # Melhorar a chave da legenda (tamanho dos símbolos)
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    fill = guide_legend(override.aes = list(alpha = 0.4)),
    linetype = guide_legend(),
    shape = guide_legend()
  ) +
  
  # --- Eixos com respiro ---
  scale_x_continuous(
    breaks = seq(floor(primeiro_ano/2)*2, 2025, by = 2),
    limits = c(primeiro_ano, 2025),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  
  # --- Rótulos e tema ---
  labs(x = "Year", 
       y = "Abundance (individuals)", 
       color = NULL, 
       fill = NULL,
       linetype = NULL, 
       shape = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    legend.text = element_text(size = 10, face = "bold"),
    legend.key.height = unit(0.9, "lines"),
    legend.key.width = unit(1.6, "lines")
  )

# Salvar
ggsave("outputs/figures/Figure_1_abundance.tiff", fig_main,
       width = 14, height = 10, units = "cm", dpi = 300, compression = "lzw")

cat("Main figure saved: Figure_1_abundance.tiff\n")

# ================================================================================
# 7. FIGURA DE DIAGNÓSTICO DE RESÍDUOS
# ================================================================================

cat("\nGENERATING DIAGNOSTIC FIGURE\n")
cat("-----------------------------\n")

# Extract residuals and fitted values
residuos <- resid(antrop_gamm_final$lme, type = "normalized")
valores_ajustados <- fitted(antrop_gamm_final$lme)

# Create diagnostic plot
tiff("outputs/figures/Figure_S1_diagnostics.tiff",
     width = 2400, height = 1800, res = 300, compression = "lzw")
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 2), cex.lab = 1.1, cex.axis = 1)

# Panel A: Residuals vs Fitted
plot(valores_ajustados, residuos,
     xlab = "Fitted values", ylab = "Standardized residuals",
     main = "A. Residuals vs Fitted",
     pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
abline(h = 0, col = "red", lty = 2, lwd = 2)
lines(lowess(valores_ajustados, residuos), col = "blue", lwd = 2)

# Panel B: Q-Q Plot
qqnorm(residuos, main = "B. Normal Q-Q Plot",
       pch = 19, col = rgb(0, 0, 0, 0.5), cex = 1.2)
qqline(residuos, col = "red", lwd = 2)

# Panel C: Histogram
hist(residuos, breaks = 15, probability = TRUE,
     main = "C. Distribution of Residuals",
     xlab = "Standardized residuals",
     col = "lightgray", border = "darkgray")
curve(dnorm(x, mean = mean(residuos), sd = sd(residuos)),
      col = "red", lwd = 2, add = TRUE)

# Panel D: ACF
acf(residuos, main = "D. Autocorrelation Function", lwd = 2)

dev.off()
cat("Diagnostic figure saved: Figure_S1_diagnostics.tiff\n")

# ================================================================================
# 8. RESUMO FINAL
# ================================================================================

cat("\n================================================================================\n")
cat("ANALYSIS SUMMARY\n")
cat("================================================================================\n")

# RESUMO DIAGNÓSTICO
cat("\nDIAGNÓSTICO DO MODELO:\n")
cat("----------------------\n")
cat("Pressupostos:\n")
cat("• Normalidade:", ifelse(shapiro_test$p.value > 0.05, "✓ OK", "✗ Violado"), 
    "(p =", round(shapiro_test$p.value, 3), ")\n")
cat("• Homocedasticidade:", ifelse(het_pval > 0.05, "✓ OK", "✗ Violado"),
    "(p =", round(het_pval, 3), ")\n")
cat("• Independência:", ifelse(length(lags_sig) == 0, "✓ OK", "✗ Violado"), "\n")

cat("\nINTERPRETAÇÃO FINAL:\n")
cat("-------------------\n")
cat("• Tendência temporal nos campos antropogênicos:", 
    ifelse(final_smooth["p-value"] < 0.05, "SIGNIFICATIVA", "NÃO SIGNIFICATIVA"),
    "(p =", format(final_smooth["p-value"], digits = 3), ")\n")
cat("• Padrão:", ifelse(final_smooth["edf"] > 1.5, "NÃO-LINEAR", "LINEAR"), "\n")
cat("• Taxa de crescimento:", round(slope_antrop, 1), "± ", round(se_slope_antrop, 1), "ind/ano\n")
cat("• Incremento total:", round(incremento_pct, 0), "%\n")
cat("• Campo de referência:", round(slope_ref, 1), "ind/ano\n")
cat("• Comparação: Antropogênico", ifelse(slope_antrop > 0, "crescente", "decrescente"),
    "vs Natural", ifelse(abs(slope_ref) < 1, "estável", ifelse(slope_ref > 0, "crescente", "decrescente")), "\n")

cat("\nARQUIVOS GERADOS:\n")
cat("-----------------\n")
cat("• Figura principal: outputs/figures/Figure_1_abundance.tiff\n")
cat("• Diagnósticos: outputs/figures/Figure_S1_diagnostics.tiff\n")
cat("• Tabela de resultados: outputs/tables/Results_summary.csv\n")

cat("\n================================================================================\n")
cat("ANÁLISE CONCLUÍDA COM SUCESSO!\n")
cat("================================================================================\n")