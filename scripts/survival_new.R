# ================================================================================
# ANÁLISE DE SOBREVIVÊNCIA: CAMPOS ANTROPOGÊNICOS VS CAMPO NATURAL
# Comparação de curvas Kaplan-Meier entre tipos de campo
# ================================================================================

library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(ggpubr)

# ================================================================================
# 1. PREPARAÇÃO DOS DADOS
# ================================================================================

cat("================================================================================\n")
cat("SURVIVAL ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

# Renomear colunas
dados_brutos <- dados_orig
colnames(dados_brutos) <- gsub("^h", "", colnames(dados_brutos))
colnames(dados_brutos)[colnames(dados_brutos) == "especies"] <- "Especie"
colnames(dados_brutos)[colnames(dados_brutos) == "ni"] <- "ID_individuo"
colnames(dados_brutos)[colnames(dados_brutos) == "area"] <- "Area"

# Identificar colunas de anos
year_cols <- grep("^[0-9]{4}$", colnames(dados_brutos), value = TRUE)
essential_cols <- c("ID_individuo", "Especie", "Area")

# Converter para formato longo
dados_long <- dados_brutos %>%
  select(all_of(essential_cols), all_of(year_cols)) %>%
  mutate(across(all_of(year_cols), as.character)) %>%
  pivot_longer(
    cols = all_of(year_cols),
    names_to = "Ano",
    values_to = "Status_Altura",
    names_transform = list(Ano = as.integer)
  ) %>%
  filter(!is.na(Especie), Especie != "")

# ================================================================================
# 2. PROCESSAR DADOS DE SOBREVIVÊNCIA
# ================================================================================

survival_data <- dados_long %>%
  filter(!is.na(Status_Altura), Status_Altura != "") %>%
  group_by(ID_individuo, Especie, Area) %>%
  arrange(Ano) %>%
  summarise(
    Ano_Entrada = min(Ano),
    Ano_Saida = max(Ano),
    Morreu = any(toupper(Status_Altura) == "M"),
    Virou_Adulto = any(toupper(Status_Altura) == "AD"),
    Ano_Morte = if(Morreu) min(Ano[toupper(Status_Altura) == "M"]) else NA_integer_,
    Ano_Adulto = if(Virou_Adulto) min(Ano[toupper(Status_Altura) == "AD"]) else NA_integer_,
    .groups = 'drop'
  ) %>%
  mutate(
    # Evento = morte; Censura = adulto ou fim do estudo
    Status_Evento = if_else(Morreu, 1, 0),
    Tempo_Final = case_when(
      Morreu & !is.na(Ano_Morte) ~ Ano_Morte,
      !Morreu & Virou_Adulto & !is.na(Ano_Adulto) ~ Ano_Adulto,
      TRUE ~ Ano_Saida
    ),
    Tempo = Tempo_Final - Ano_Entrada,
    # CRIAR VARIÁVEL DE TIPO DE CAMPO
    Tipo_Campo = factor(
      ifelse(Area == 3, "Natural grassland", "Abandoned pastures"),
      levels = c("Abandoned pastures", "Natural grassland")
    )
  ) %>%
  filter(Tempo >= 0)

cat("DATA SUMMARY\n")
cat("------------\n")
cat("Total individuals:", nrow(survival_data), "\n")
cat("Abandoned pastures (Areas 1,2,4):", sum(survival_data$Tipo_Campo == "Abandoned pastures"), "\n")
cat("Natural grassland (Area 3):", sum(survival_data$Tipo_Campo == "Natural grassland"), "\n")
cat("Total deaths:", sum(survival_data$Status_Evento), "\n")
cat("Censored observations:", sum(survival_data$Status_Evento == 0), "\n\n")

# ================================================================================
# 3. ANÁLISE COMPARATIVA POR TIPO DE CAMPO
# ================================================================================

cat("COMPARATIVE ANALYSIS\n")
cat("-------------------\n")

# Criar objeto de sobrevivência
surv_obj <- Surv(time = survival_data$Tempo, event = survival_data$Status_Evento)

# Ajustar modelo Kaplan-Meier por tipo de campo
fit_km_tipo <- survfit(surv_obj ~ Tipo_Campo, data = survival_data)

# Teste Log-Rank
log_rank_test <- survdiff(surv_obj ~ Tipo_Campo, data = survival_data)
cat("Log-Rank Test:\n")
print(log_rank_test)

# ================================================================================
# 4. FIGURA PRINCIPAL - COMPARAÇÃO ENTRE TIPOS
# ================================================================================

cat("\nGENERATING MAIN COMPARISON FIGURE\n")
cat("----------------------------------\n")

# Figura principal comparando tipos de campo
fig_comparison <- ggsurvplot(
  fit_km_tipo,
  data = survival_data,
  title = "Survival curves: Abandoned pastures vs Natural grassland",
  conf.int = TRUE,
  conf.int.alpha = 0.2,
  pval = TRUE,
  pval.method = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  size = 1.2,
  palette = c("#D55E00", "#009E73"),
  legend.title = "Field type",
  legend.labs = c("Abandoned pastures", "Natural grassland"),
  xlab = "Time (years)",
  ylab = "Survival probability",
  break.time.by = 2,
  risk.table.height = 0.25,
  ggtheme = theme_classic(base_size = 11)
)

print(fig_comparison)
fig_comparison$data.survplot
fig_comparison$data.survtable
fig_comparison$plot
fig_comparison$table
# ================================================================================
# 5. ANÁLISE POR ESPÉCIE DENTRO DE CADA TIPO
# ================================================================================

cat("\nSPECIES-LEVEL ANALYSIS\n")
cat("----------------------\n")

# Top 5 espécies por tipo de campo
top_species_by_type <- survival_data %>%
  group_by(Tipo_Campo, Especie) %>%
  summarise(n_individuos = n(), .groups = 'drop') %>%
  group_by(Tipo_Campo) %>%
  slice_max(order_by = n_individuos, n = 5) %>%
  ungroup()

# Lista para armazenar gráficos
species_plots <- list()

for(tipo in unique(survival_data$Tipo_Campo)) {
  
  cat(paste("\nProcessing:", tipo, "\n"))
  
  # Filtrar top espécies
  top_species <- top_species_by_type %>%
    filter(Tipo_Campo == tipo) %>%
    pull(Especie)
  
  if(length(top_species) < 2) next
  
  # Dados filtrados
  data_tipo <- survival_data %>%
    filter(Tipo_Campo == tipo, Especie %in% top_species) %>%
    mutate(Especie = factor(Especie))
  
  if(nrow(data_tipo) < 5) next
  
  # Ajustar KM por espécie
  surv_obj_species <- Surv(time = data_tipo$Tempo, event = data_tipo$Status_Evento)
  fit_km_species <- survfit(surv_obj_species ~ Especie, data = data_tipo)
  
  # Teste log-rank entre espécies
  if(sum(data_tipo$Status_Evento) > 0) {
    lr_test_species <- survdiff(surv_obj_species ~ Especie, data = data_tipo)
    cat(sprintf("  Species differences (p-value): %.4f\n", 
                1 - pchisq(lr_test_species$chisq, length(lr_test_species$n) - 1)))
  }
  
  # Criar gráfico
  plot_species <- ggsurvplot(
    fit_km_species,
    data = data_tipo,
    title = paste("Top species -", tipo),
    conf.int = FALSE,
    pval = TRUE,
    risk.table = FALSE,
    size = 0.8,
    palette = "viridis",
    legend.title = "Species",
    xlab = "Time (years)",
    ylab = "Survival probability",
    ggtheme = theme_classic(base_size = 10),
    legend = "right"
  )
  
  species_plots[[tipo]] <- plot_species$plot
}

# ================================================================================
# 6. PAINEL COMBINADO
# ================================================================================

if(length(species_plots) > 0) {
  cat("\nGENERATING COMBINED PANEL\n")
  cat("-------------------------\n")
  
  # Combinar figura principal com análises por espécie
  combined_panel <- ggarrange(
    fig_comparison$plot,
    ggarrange(plotlist = species_plots, ncol = 2, labels = c("b", "c")),
    nrow = 2,
    labels = c("a", ""),
    heights = c(1.2, 1)
  )
  
  combined_panel <- annotate_figure(
    combined_panel,
    top = text_grob("Survival Analysis: Field Type Comparison", 
                    face = "bold", size = 14)
  )
  
  print(combined_panel)
  
  # Salvar
  ggsave("outputs/figures/Figure_5_survival_comparison.tiff", 
         plot = combined_panel, 
         width = 12, height = 10, dpi = 300, compression = "lzw")
}

# ================================================================================
# 7. TABELA DE RESULTADOS
# ================================================================================

cat("\n=== SURVIVAL STATISTICS ===\n")

# Estatísticas por tipo de campo
survival_stats <- survival_data %>%
  group_by(Tipo_Campo) %>%
  summarise(
    N_individuals = n(),
    N_deaths = sum(Status_Evento),
    Mortality_rate = mean(Status_Evento) * 100,
    Median_time = median(Tempo),
    Mean_time = mean(Tempo),
    SD_time = sd(Tempo),
    .groups = 'drop'
  )

print(survival_stats)

# Median survival time from KM
median_survival <- summary(fit_km_tipo)$table[,"median"]
cat("\nMedian survival time (Kaplan-Meier):\n")
print(median_survival)

# Salvar tabela
write.csv(survival_stats, "outputs/tables/Survival_statistics.csv", row.names = FALSE)

# ================================================================================
# 8. MODELO DE COX (ANÁLISE ADICIONAL)
# ================================================================================

cat("\n=== COX PROPORTIONAL HAZARDS MODEL ===\n")

# Modelo de Cox comparando tipos
cox_model <- coxph(surv_obj ~ Tipo_Campo, data = survival_data)
print(summary(cox_model))

# Hazard ratio
hr <- exp(coef(cox_model))
ci <- exp(confint(cox_model))
cat(sprintf("\nHazard Ratio (Natural vs Abandoned): %.2f (95%% CI: %.2f-%.2f)\n",
            hr, ci[1], ci[2]))

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETED\n")
cat("================================================================================\n")


# ================================================================================
# FIGURA DE SOBREVIVÊNCIA PARA PUBLICAÇÃO - CORREÇÃO FINAL
# ================================================================================

library(survival)
library(survminer)
library(ggplot2)

# Criar figura base
fig_publication <- ggsurvplot(
  fit_km_tipo,
  data = survival_data,
  
  conf.int = TRUE,
  conf.int.alpha = 0.15,
  palette = c("#D55E00", "#009E73"),
  
  pval = TRUE,
  pval.coord = c(7, 0.15),
  pval.size = 4,
  
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  risk.table.fontsize = 3.5,
  
  xlab = "Time (years)",
  ylab = "Survival probability",
  xlim = c(0, 11),
  break.time.by = 1,
  
  legend = c(0.75, 0.85),
  legend.title = "",
  legend.labs = c("Abandoned pastures (n=4,967)", "Natural grassland (n=68)"),
  
  size = 1.5,
  censor = FALSE,
  
  ggtheme = theme_classic(base_size = 11)
)

# Customizar o plot principal
fig_publication$plot <- fig_publication$plot +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank()
  ) +
  # Adicionar anotação
  annotate("text", 
           x = 2.5, y = 0.5,
           label = "86% mortality\nin 2 years",
           color = "#009E73", 
           size = 3.5, 
           fontface = "italic") +
  annotate("segment", 
           x = 2.3, xend = 2, 
           y = 0.45, yend = 0.15,
           arrow = arrow(length = unit(0.15, "cm")),
           color = "#009E73", 
           alpha = 0.6)

# FORMA CORRETA DE SALVAR ggsurvplot COM TABELA DE RISCO
# Opção 1: Usar a função própria do survminer
pdf("outputs/figures/Figure_Survival_Publication.pdf", width = 8, height = 7)
print(fig_publication)
dev.off()

# Opção 2: Extrair e combinar os componentes
library(gridExtra)
combined_plot <- arrange_ggsurvplots(
  list(fig_publication),
  print = FALSE,
  ncol = 1, nrow = 1
)

ggsave(
  filename = "outputs/figures/Figure_Survival_Publication.tiff",
  plot = combined_plot,
  width = 8,
  height = 7,
  dpi = 300,
  compression = "lzw"
)

# ================================================================================
# ALTERNATIVA MAIS SIMPLES: APENAS O PLOT PRINCIPAL (SEM TABELA)
# ================================================================================

# Extrair apenas o gráfico principal
plot_only <- fig_publication$plot

# Salvar apenas o gráfico (sem tabela de risco)
ggsave(
  filename = "outputs/figures/Figure_Survival_NoTable.tiff",
  plot = plot_only,
  width = 7,
  height = 5,
  dpi = 300,
  compression = "lzw"
)

# ================================================================================
# VERSÃO MANUAL: COMBINAR PLOT E TABELA
# ================================================================================

# Se quiser controle total sobre o layout
library(cowplot)

# Extrair componentes
survival_plot <- fig_publication$plot
risk_table <- fig_publication$table

# Combinar manualmente
final_combined <- plot_grid(
  survival_plot,
  risk_table,
  ncol = 1,
  align = "v",
  rel_heights = c(3, 1),
  labels = c("", ""),
  label_size = 12
)

# Salvar versão combinada manual
ggsave(
  filename = "outputs/figures/Figure_Survival_Combined.tiff",
  plot = final_combined,
  width = 8,
  height = 7,
  dpi = 300,
  compression = "lzw"
)

print(final_combined)
