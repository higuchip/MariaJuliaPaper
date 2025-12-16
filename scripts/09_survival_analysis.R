# ==============================================================================
# SURVIVAL ANALYSIS: ABANDONED PASTURES VS NATURAL GRASSLAND
# ==============================================================================
#
# Description: Compares survival curves between abandoned pastures and natural
#              grassland using Kaplan-Meier estimators and Cox proportional
#              hazards models. 
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
# Input:       dados_orig from 01_data_processing.R
#
# Output:      ../artigo/biotropica/review/Figure_6_survival.tiff
#              outputs/tables/Survival_statistics.csv
#
# Note:        Requires dados_orig from 01_data_processing.R
#              R project located in MariaJulia/dissertacao
#
# ==============================================================================

# --- Setup --------------------------------------------------------------------

set.seed(123)

# Required packages
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
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
cat("SURVIVAL ANALYSIS - ", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

cat("1. DATA PREPARATION\n")
cat("-------------------\n\n")

# Check dependency
if(!exists("dados_orig")) {
  stop("Object 'dados_orig' not found. Please run 01_data_processing.R first.")
}

# Rename columns
dados_brutos <- dados_orig
colnames(dados_brutos) <- gsub("^h", "", colnames(dados_brutos))
colnames(dados_brutos)[colnames(dados_brutos) == "especies"] <- "Especie"
colnames(dados_brutos)[colnames(dados_brutos) == "ni"] <- "ID_individuo"
colnames(dados_brutos)[colnames(dados_brutos) == "area"] <- "Area"

# Identify year columns
year_cols <- grep("^[0-9]{4}$", colnames(dados_brutos), value = TRUE)
essential_cols <- c("ID_individuo", "Especie", "Area")

# Convert to long format
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

# ==============================================================================
# 2. PROCESS SURVIVAL DATA
# ==============================================================================

cat("Processing survival data...\n")

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
    # Event = death; Censored = adult or end of study
    Status_Evento = if_else(Morreu, 1, 0),
    Tempo_Final = case_when(
      Morreu & !is.na(Ano_Morte) ~ Ano_Morte,
      !Morreu & Virou_Adulto & !is.na(Ano_Adulto) ~ Ano_Adulto,
      TRUE ~ Ano_Saida
    ),
    Tempo = Tempo_Final - Ano_Entrada,
    # Create field type variable
    Tipo_Campo = factor(
      ifelse(Area == 3, "Natural grassland", "Abandoned pastures"),
      levels = c("Abandoned pastures", "Natural grassland")
    )
  ) %>%
  filter(Tempo >= 0)

cat("\nDATA SUMMARY\n")
cat("------------\n")
cat("Total individuals:", nrow(survival_data), "\n")
cat("Abandoned pastures (Areas 1,2,4):", sum(survival_data$Tipo_Campo == "Abandoned pastures"), "\n")
cat("Natural grassland (Area 3):", sum(survival_data$Tipo_Campo == "Natural grassland"), "\n")
cat("Total deaths:", sum(survival_data$Status_Evento), "\n")
cat("Censored observations:", sum(survival_data$Status_Evento == 0), "\n\n")

# ==============================================================================
# 3. KAPLAN-MEIER ANALYSIS
# ==============================================================================

cat("================================================================================\n")
cat("3. KAPLAN-MEIER ANALYSIS\n")
cat("================================================================================\n\n")

# Create survival object
surv_obj <- Surv(time = survival_data$Tempo, event = survival_data$Status_Evento)

# Fit Kaplan-Meier by field type
fit_km_tipo <- survfit(surv_obj ~ Tipo_Campo, data = survival_data)

# Log-Rank test
log_rank_test <- survdiff(surv_obj ~ Tipo_Campo, data = survival_data)
cat("Log-Rank Test:\n")
print(log_rank_test)

# Extract p-value
log_rank_p <- 1 - pchisq(log_rank_test$chisq, df = 1)
cat(sprintf("\nLog-Rank p-value: %.2e\n", log_rank_p))

# ==============================================================================
# 4. COX PROPORTIONAL HAZARDS MODEL
# ==============================================================================

cat("\n================================================================================\n")
cat("4. COX PROPORTIONAL HAZARDS MODEL\n")
cat("================================================================================\n\n")

cox_model <- coxph(surv_obj ~ Tipo_Campo, data = survival_data)
print(summary(cox_model))

# Hazard ratio
hr <- exp(coef(cox_model))
ci <- exp(confint(cox_model))
cat(sprintf("\nHazard Ratio (Natural vs Abandoned): %.2f (95%% CI: %.2f-%.2f)\n",
            hr, ci[1], ci[2]))

# ==============================================================================
# 5. SURVIVAL STATISTICS TABLE
# ==============================================================================

cat("\n================================================================================\n")
cat("5. SURVIVAL STATISTICS\n")
cat("================================================================================\n\n")

# Statistics by field type
survival_stats <- survival_data %>%
  group_by(Tipo_Campo) %>%
  summarise(
    N_individuals = n(),
    N_deaths = sum(Status_Evento),
    Mortality_rate_pct = round(mean(Status_Evento) * 100, 1),
    Median_time = median(Tempo),
    Mean_time = round(mean(Tempo), 2),
    SD_time = round(sd(Tempo), 2),
    .groups = 'drop'
  )

print(survival_stats)

# Median survival time from KM
cat("\nMedian survival time (Kaplan-Meier):\n")
print(summary(fit_km_tipo)$table[,"median"])

# Save table
write.csv(survival_stats, "outputs/tables/Survival_statistics.csv", row.names = FALSE)

# Calculate 2-year mortality for both ecosystem types
mortality_2yr <- survival_data %>%
  group_by(Tipo_Campo) %>%
  summarise(
    total = n(),
    deaths_2yr = sum(Status_Evento == 1 & Tempo <= 2),
    mortality_2yr_pct = round(deaths_2yr / total * 100, 1),
    .groups = 'drop'
  )

cat("\n2-YEAR MORTALITY:\n")
print(mortality_2yr)

mort_abandoned <- mortality_2yr %>% 
  filter(Tipo_Campo == "Abandoned pastures") %>% 
  pull(mortality_2yr_pct)

mort_natural <- mortality_2yr %>% 
  filter(Tipo_Campo == "Natural grassland") %>% 
  pull(mortality_2yr_pct)

cat(sprintf("\nAbandoned pastures: %.1f%% mortality in 2 years\n", mort_abandoned))
cat(sprintf("Natural grassland: %.1f%% mortality in 2 years\n", mort_natural))

# ==============================================================================
# 6. PUBLICATION FIGURE
# ==============================================================================

cat("\n================================================================================\n")
cat("6. GENERATING PUBLICATION FIGURE\n")
cat("================================================================================\n\n")

# Get sample sizes for legend
n_abandoned <- sum(survival_data$Tipo_Campo == "Abandoned pastures")
n_natural <- sum(survival_data$Tipo_Campo == "Natural grassland")

# Create figure (without risk table)
fig_publication <- ggsurvplot(
  fit_km_tipo,
  data = survival_data,
  
  # Confidence intervals
  conf.int = TRUE,
  conf.int.alpha = 0.15,
  palette = c("#D55E00", "#009E73"),
  
  # P-value
  pval = TRUE,
  pval.coord = c(7, 0.15),
  pval.size = 4,
  
  # No risk table
  risk.table = FALSE,
  
  # Axes
  xlab = "Time (years)",
  ylab = "Survival probability",
  xlim = c(0, 11),
  break.time.by = 1,
  
  # Legend
  legend = c(0.75, 0.85),
  legend.title = "",
  legend.labs = c(paste0("Abandoned pastures (n=", n_abandoned, ")"), 
                  paste0("Natural grassland (n=", n_natural, ")")),
  
  # Line aesthetics
  size = 1.5,
  censor = FALSE,
  
  # Theme
  ggtheme = theme_classic(base_size = 11)
)

# Customize plot with calculated mortality values
fig_survival <- fig_publication$plot +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank()
  ) +
  # Annotation for natural grassland (high mortality)
  annotate("text", 
           x = 2.5, y = 0.45,
           label = paste0(86.4, "% mortality\nin 2 years"),
           color = "#009E73", 
           size = 3.5, 
           fontface = "italic") +
  annotate("segment", 
           x = 2.3, xend = 2, 
           y = 0.40, yend = 0.15,
           arrow = arrow(length = unit(0.15, "cm")),
           color = "#009E73", 
           alpha = 0.6) 

# Save figure
ggsave(
  filename = paste0(output_path, "Fig3.tiff"),
  plot = fig_survival,
  width = 14,
  height = 10,
  units = "cm",
  dpi = 300,
  compression = "lzw"
)

cat("Figure saved:", paste0(output_path, "Fig3.tiff"), "\n")

# Calculate survival at specific time points
surv_summary <- summary(fit_km_tipo, times = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

cat("\n================================================================================\n")
cat("7. DETAILED RESULTS\n")
cat("================================================================================\n\n")

# --- Sample sizes ---
cat("SAMPLE SIZES:\n")
cat(sprintf("  Abandoned pastures: n = %d individuals\n", n_abandoned))
cat(sprintf("  Natural grassland: n = %d individuals\n", n_natural))

# --- Log-rank test ---
cat("\nLOG-RANK TEST:\n")
cat(sprintf("  Chi-squared = %.1f\n", log_rank_test$chisq))
cat(sprintf("  p-value = %.4e\n", log_rank_p))

# --- Median survival times ---
median_surv_table <- summary(fit_km_tipo)$table
cat("\nMEDIAN SURVIVAL TIME (from Kaplan-Meier):\n")
cat(sprintf("  Abandoned pastures: %.0f years\n", median_surv_table["Tipo_Campo=Abandoned pastures", "median"]))
cat(sprintf("  Natural grassland: %.0f years\n", median_surv_table["Tipo_Campo=Natural grassland", "median"]))

# --- Survival probabilities at key time points (Kaplan-Meier estimates) ---
cat("\nSURVIVAL PROBABILITIES (Kaplan-Meier estimates):\n")
cat("------------------------------------------------\n")

# Create data frame with survival estimates
surv_df <- data.frame(
  Time = surv_summary$time,
  Strata = gsub("Tipo_Campo=", "", surv_summary$strata),
  n_risk = surv_summary$n.risk,
  n_event = surv_summary$n.event,
  Survival_prob = surv_summary$surv,
  Survival_pct = round(surv_summary$surv * 100, 1),
  SE = surv_summary$std.err,
  Lower_95CI = surv_summary$lower,
  Upper_95CI = surv_summary$upper
)

print(surv_df)

# --- Extract specific values for manuscript ---
cat("\n--- VALUES FOR MANUSCRIPT TEXT ---\n\n")

# Abandoned pastures
surv_abandoned <- surv_df %>% filter(Strata == "Abandoned pastures")

surv_ab_2yr <- surv_abandoned %>% filter(Time == 2) %>% pull(Survival_pct)
surv_ab_10yr <- surv_abandoned %>% filter(Time == 10) %>% pull(Survival_pct)

cat("ABANDONED PASTURES:\n")
cat(sprintf("  Survival at year 2: %.1f%%\n", surv_ab_2yr))
cat(sprintf("  Survival at year 10: %.1f%%\n", surv_ab_10yr))
cat(sprintf("  Median survival: %.0f years\n", median_surv_table["Tipo_Campo=Abandoned pastures", "median"]))

# Natural grassland
surv_natural_df <- surv_df %>% filter(Strata == "Natural grassland")

surv_nat_2yr <- surv_natural_df %>% filter(Time == 2) %>% pull(Survival_pct)

# Find year of complete mortality (survival = 0 or first NA)
year_complete_mortality <- surv_natural_df %>% 
  filter(Survival_prob == 0 | n_risk == 0) %>% 
  slice(1) %>% 
  pull(Time)

if(length(year_complete_mortality) == 0) {
  # Check for very low survival
  year_complete_mortality <- surv_natural_df %>% 
    filter(Survival_pct < 1) %>% 
    slice(1) %>% 
    pull(Time)
}

cat("\nNATURAL GRASSLAND:\n")
cat(sprintf("  Survival at year 2: %.1f%%\n", surv_nat_2yr))
cat(sprintf("  Year of complete mortality: %s\n", 
            ifelse(length(year_complete_mortality) > 0, 
                   paste0("Year ", year_complete_mortality), 
                   "Not reached")))
cat(sprintf("  Median survival: %.0f years\n", median_surv_table["Tipo_Campo=Natural grassland", "median"]))

# --- Cox model results ---
cat("\nCOX PROPORTIONAL HAZARDS:\n")
cat(sprintf("  Hazard Ratio (Natural vs Abandoned): %.2f\n", hr))
cat(sprintf("  95%% CI: %.2f - %.2f\n", ci[1], ci[2]))


# Save detailed table
write.csv(surv_df, "outputs/tables/Survival_KaplanMeier_estimates.csv", row.names = FALSE)
cat("\nDetailed estimates saved: outputs/tables/Survival_KaplanMeier_estimates.csv\n")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("================================================================================\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================