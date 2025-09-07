# ================================================================================
# COMPLETE PANEL FIGURE - TEMPORAL DYNAMICS
# 
# Four panels: Abundance, Height, Richness, Similarity
# ================================================================================

set.seed(123)

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(nlme)
library(lme4)
library(vegan)
library(patchwork)
library(cowplot)

# Create output directories
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", showWarnings = FALSE)

cat("================================================================================\n")
cat("GENERATING COMBINED PANEL FIGURE - COMPLETE PROCESSING\n")
cat("================================================================================\n\n")

# ================================================================================
# 1. PREPARE ALL BASE DATA
# ================================================================================

# Identify height columns
height_cols <- grep("^h\\d{4}$", names(df_processed), value = TRUE)

# Long format for all analyses
data_long <- df_processed %>%
  select(area, parc, ni, especies, all_of(height_cols)) %>%
  pivot_longer(
    cols = all_of(height_cols),
    names_to = "Ano",
    values_to = "Altura",
    names_prefix = "h"
  ) %>%
  mutate(Ano = as.numeric(Ano)) %>%
  filter(!is.na(Altura))

primeiro_ano <- min(data_long$Ano)
ultimo_ano <- max(data_long$Ano)

# ================================================================================
# 2. PROCESS ABUNDANCE DATA
# ================================================================================

cat("Processing ABUNDANCE data...\n")

# Calculate abundance
abundance_data <- data_long %>%
  group_by(Area = as.factor(area), Ano) %>%
  summarise(Abundancia = n(), .groups = 'drop')

# Separate datasets
antrop_abundance_data <- abundance_data %>%
  filter(Area %in% c("1", "2", "4")) %>%
  mutate(Area_Antrop = factor(Area))

reference_abundance_data <- abundance_data %>%
  filter(Area == "3")

# Fit models
k_abund <- min(length(unique(antrop_abundance_data$Ano)) - 1, 4)
gamm_abundance <- gamm(Abundancia ~ s(Ano, k = k_abund),
                       random = list(Area_Antrop = ~1),
                       data = antrop_abundance_data,
                       method = "REML")

gam_abundance_ref <- gam(Abundancia ~ s(Ano, k = 3),
                         data = reference_abundance_data,
                         method = "REML")

# Create predictions
anos_pred <- seq(primeiro_ano, ultimo_ano, length.out = 100)

# Abandoned pastures
pred_antrop_abund <- data.frame(Ano = anos_pred, Area_Antrop = "1")
pred_vals <- predict(gamm_abundance$gam, newdata = pred_antrop_abund, se.fit = TRUE)
pred_antrop_abund$fit <- pred_vals$fit
pred_antrop_abund$lower <- pred_antrop_abund$fit - 1.96 * pred_vals$se.fit
pred_antrop_abund$upper <- pred_antrop_abund$fit + 1.96 * pred_vals$se.fit
pred_antrop_abund$Grupo <- "Abandoned pastures"

# Natural grassland
pred_ref_abund <- data.frame(Ano = anos_pred)
pred_vals <- predict(gam_abundance_ref, newdata = pred_ref_abund, se.fit = TRUE)
pred_ref_abund$fit <- pred_vals$fit
pred_ref_abund$lower <- pred_ref_abund$fit - 1.96 * pred_vals$se.fit
pred_ref_abund$upper <- pred_ref_abund$fit + 1.96 * pred_vals$se.fit
pred_ref_abund$Grupo <- "Natural grassland"

# Annual means
annual_means_abund <- antrop_abundance_data %>%
  group_by(Ano) %>%
  summarise(mean = mean(Abundancia), se = sd(Abundancia)/sqrt(n()), .groups = 'drop') %>%
  mutate(Grupo = "Abandoned pastures")

reference_abundance_data$Grupo <- "Natural grassland"

# ================================================================================
# 3. PROCESS HEIGHT DATA
# ================================================================================

cat("Processing HEIGHT data...\n")

# Mean height by area and year
height_data <- data_long %>%
  group_by(Area = as.factor(area), Ano) %>%
  summarise(AlturaMedia = mean(Altura, na.rm = TRUE), .groups = 'drop')

antrop_height_data <- height_data %>%
  filter(Area %in% c("1", "2", "4")) %>%
  mutate(Area_Antrop = factor(Area))

reference_height_data <- height_data %>%
  filter(Area == "3")

# Fit models
k_height <- min(length(unique(antrop_height_data$Ano)) - 1, 4)
gamm_height <- gamm(AlturaMedia ~ s(Ano, k = k_height),
                    random = list(Area_Antrop = ~1),
                    data = antrop_height_data,
                    method = "REML")

gam_height_ref <- gam(AlturaMedia ~ s(Ano, k = 3),
                      data = reference_height_data,
                      method = "REML")

# Predictions
pred_antrop_height <- data.frame(Ano = anos_pred, Area_Antrop = "1")
pred_vals <- predict(gamm_height$gam, newdata = pred_antrop_height, se.fit = TRUE)
pred_antrop_height$fit <- pred_vals$fit
pred_antrop_height$lower <- pred_antrop_height$fit - 1.96 * pred_vals$se.fit
pred_antrop_height$upper <- pred_antrop_height$fit + 1.96 * pred_vals$se.fit
pred_antrop_height$Grupo <- "Abandoned pastures"

pred_ref_height <- data.frame(Ano = anos_pred)
pred_vals <- predict(gam_height_ref, newdata = pred_ref_height, se.fit = TRUE)
pred_ref_height$fit <- pred_vals$fit
pred_ref_height$lower <- pred_ref_height$fit - 1.96 * pred_vals$se.fit
pred_ref_height$upper <- pred_ref_height$fit + 1.96 * pred_vals$se.fit
pred_ref_height$Grupo <- "Natural grassland"

annual_means_height <- antrop_height_data %>%
  group_by(Ano) %>%
  summarise(mean = mean(AlturaMedia), se = sd(AlturaMedia)/sqrt(n()), .groups = 'drop') %>%
  mutate(Grupo = "Abandoned pastures")

reference_height_data$Grupo <- "Natural grassland"

# ================================================================================
# 4. PROCESS RICHNESS DATA
# ================================================================================

cat("Processing RICHNESS data...\n")

# Species richness
richness_data <- data_long %>%
  group_by(Area = as.factor(area), Ano) %>%
  summarise(Riqueza = n_distinct(especies), .groups = 'drop')

antrop_richness_data <- richness_data %>%
  filter(Area %in% c("1", "2", "4")) %>%
  mutate(Area_Antrop = factor(Area))

reference_richness_data <- richness_data %>%
  filter(Area == "3")

# Fit models
k_rich <- min(length(unique(antrop_richness_data$Ano)) - 1, 4)
gamm_richness <- gamm(Riqueza ~ s(Ano, k = k_rich),
                      random = list(Area_Antrop = ~1),
                      data = antrop_richness_data,
                      method = "REML")

gam_richness_ref <- gam(Riqueza ~ s(Ano, k = 3),
                        data = reference_richness_data,
                        method = "REML")

# Predictions
pred_antrop_richness <- data.frame(Ano = anos_pred, Area_Antrop = "1")
pred_vals <- predict(gamm_richness$gam, newdata = pred_antrop_richness, se.fit = TRUE)
pred_antrop_richness$fit <- pred_vals$fit
pred_antrop_richness$lower <- pred_antrop_richness$fit - 1.96 * pred_vals$se.fit
pred_antrop_richness$upper <- pred_antrop_richness$fit + 1.96 * pred_vals$se.fit
pred_antrop_richness$Grupo <- "Abandoned pastures"

pred_ref_richness <- data.frame(Ano = anos_pred)
pred_vals <- predict(gam_richness_ref, newdata = pred_ref_richness, se.fit = TRUE)
pred_ref_richness$fit <- pred_vals$fit
pred_ref_richness$lower <- pred_ref_richness$fit - 1.96 * pred_vals$se.fit
pred_ref_richness$upper <- pred_ref_richness$fit + 1.96 * pred_vals$se.fit
pred_ref_richness$Grupo <- "Natural grassland"

annual_means_richness <- antrop_richness_data %>%
  group_by(Ano) %>%
  summarise(mean = mean(Riqueza), se = sd(Riqueza)/sqrt(n()), .groups = 'drop') %>%
  mutate(Grupo = "Abandoned pastures")

reference_richness_data$Grupo <- "Natural grassland"

# ================================================================================
# 5. PROCESS SIMILARITY DATA
# ================================================================================

cat("Processing SIMILARITY data...\n")

# Load reference data
referencia <- read.table("dados/dados_talissa.csv", header=TRUE, dec=".", sep=",")
dados_ref <- referencia %>% 
  filter(DAP2 > 0, Area %in% c(1, 2, 3, 4)) %>%
  select(Area, Especie)

# Function to calculate Jaccard
calcular_jaccard <- function(area_alvo, ano_alvo) {
  if(area_alvo == 3) return(0)
  
  especies_ref <- dados_ref %>%
    filter(Area == area_alvo) %>%
    pull(Especie) %>%
    unique()
  
  nome_coluna <- paste0("h", ano_alvo)
  if (!(nome_coluna %in% names(df_processed))) return(NA)
  
  especies_reg <- df_processed %>%
    filter(area == area_alvo) %>%
    filter(!is.na(!!sym(nome_coluna))) %>%
    filter(!!sym(nome_coluna) > 0) %>%
    pull(especies) %>%
    unique()
  
  if (length(especies_reg) == 0) return(0)
  
  intersecao <- length(intersect(especies_ref, especies_reg))
  uniao <- length(union(especies_ref, especies_reg))
  
  if (uniao == 0) return(0)
  return(intersecao / uniao)
}

# Calculate similarity
anos_sim <- unique(data_long$Ano)
similarity_data <- expand.grid(Area = c(1, 2, 3, 4), Ano = anos_sim) %>%
  rowwise() %>%
  mutate(Similaridade = calcular_jaccard(Area, Ano)) %>%
  filter(!is.na(Similaridade))

antrop_similarity_data <- similarity_data %>%
  filter(Area %in% c(1, 2, 4))

reference_similarity_data <- similarity_data %>%
  filter(Area == 3)

# Since similarity is stable, use mean for predictions
mean_sim <- mean(antrop_similarity_data$Similaridade)
se_sim <- sd(antrop_similarity_data$Similaridade)/sqrt(nrow(antrop_similarity_data))

pred_antrop_sim <- data.frame(
  Ano = anos_pred,
  fit = rep(mean_sim, length(anos_pred)),
  lower = rep(mean_sim - 1.96*se_sim, length(anos_pred)),
  upper = rep(mean_sim + 1.96*se_sim, length(anos_pred)),
  Grupo = "Abandoned pastures"
)

pred_ref_sim <- data.frame(
  Ano = anos_pred,
  fit = 0,
  lower = 0,
  upper = 0,
  Grupo = "Natural grassland"
)

annual_means_sim <- antrop_similarity_data %>%
  group_by(Ano) %>%
  summarise(mean = mean(Similaridade), se = sd(Similaridade)/sqrt(n()), .groups = 'drop') %>%
  mutate(Grupo = "Abandoned pastures")

reference_similarity_data$Grupo <- "Natural grassland"

# ================================================================================
# 6. CREATE PANELS
# ================================================================================

cat("\nCreating panels...\n")

# Define colors for consistency
color_values <- c("Abandoned pastures" = "#D55E00",
                  "Natural grassland" = "#009E73")
fill_values <- color_values
linetype_values <- c("Abandoned pastures" = "solid",
                     "Natural grassland" = "23")
shape_values <- c("Abandoned pastures" = 21,
                  "Natural grassland" = 23)

# Common theme elements
theme_panel <- theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black", size = 9),
    axis.title.y = element_text(color = "black", size = 10),
    axis.title.x = element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 11, hjust = -0.1),
    legend.position = "none"
  )

# PANEL A: ABUNDANCE
panel_a <- ggplot() +
  geom_ribbon(data = pred_antrop_abund, aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo), alpha = 0.22) +
  geom_ribbon(data = pred_ref_abund, aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo), alpha = 0.18) +
  geom_line(data = pred_antrop_abund, aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), linewidth = 1.5) +
  geom_line(data = pred_ref_abund, aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), linewidth = 1.2) +
  geom_errorbar(data = annual_means_abund, aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.3, alpha = 0.7, linewidth = 0.6, show.legend = FALSE) +
  geom_point(data = annual_means_abund, aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.5, stroke = 0.8) +
  geom_point(data = reference_abundance_data, aes(x = Ano, y = Abundancia, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = fill_values) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = seq(2014, 2025, by = 2), limits = c(primeiro_ano, 2025)) +
  labs(x = NULL, y = "Abundance (individuals)") +
  ggtitle("(a)") +
  theme_panel

# PANEL B: HEIGHT
panel_b <- ggplot() +
  geom_ribbon(data = pred_antrop_height, aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo), alpha = 0.22) +
  geom_ribbon(data = pred_ref_height, aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo), alpha = 0.18) +
  geom_line(data = pred_antrop_height, aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), linewidth = 1.5) +
  geom_line(data = pred_ref_height, aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), linewidth = 1.2) +
  geom_errorbar(data = annual_means_height, aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.3, alpha = 0.7, linewidth = 0.6, show.legend = FALSE) +
  geom_point(data = annual_means_height, aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.5, stroke = 0.8) +
  geom_point(data = reference_height_data, aes(x = Ano, y = AlturaMedia, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = fill_values) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = seq(2014, 2025, by = 2), limits = c(primeiro_ano, 2025)) +
  labs(x = NULL, y = "Mean height (m)") +
  ggtitle("(b)") +
  theme_panel

# PANEL C: RICHNESS
panel_c <- ggplot() +
  geom_ribbon(data = pred_antrop_richness, aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo), alpha = 0.22) +
  geom_ribbon(data = pred_ref_richness, aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo), alpha = 0.18) +
  geom_line(data = pred_antrop_richness, aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), linewidth = 1.5) +
  geom_line(data = pred_ref_richness, aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), linewidth = 1.2) +
  geom_errorbar(data = annual_means_richness, aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.3, alpha = 0.7, linewidth = 0.6, show.legend = FALSE) +
  geom_point(data = annual_means_richness, aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.5, stroke = 0.8) +
  geom_point(data = reference_richness_data, aes(x = Ano, y = Riqueza, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = fill_values) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = seq(2014, 2025, by = 2), limits = c(primeiro_ano, 2025)) +
  labs(x = "Year", y = "Species richness") +
  ggtitle("(c)") +
  theme_panel

# PANEL D: SIMILARITY
panel_d <- ggplot() +
  geom_ribbon(data = pred_antrop_sim, aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo), alpha = 0.22) +
  geom_ribbon(data = pred_ref_sim, aes(x = Ano, ymin = lower, ymax = upper, fill = Grupo), alpha = 0.18) +
  geom_line(data = pred_antrop_sim, aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), linewidth = 1.5) +
  geom_line(data = pred_ref_sim, aes(x = Ano, y = fit, color = Grupo, linetype = Grupo), linewidth = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50", alpha = 0.5) +
  geom_errorbar(data = annual_means_sim, aes(x = Ano, ymin = mean - se, ymax = mean + se, color = Grupo),
                width = 0.3, alpha = 0.7, linewidth = 0.6, show.legend = FALSE) +
  geom_point(data = annual_means_sim, aes(x = Ano, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2.5, stroke = 0.8) +
  geom_point(data = reference_similarity_data, aes(x = Ano, y = Similaridade, color = Grupo, fill = Grupo, shape = Grupo),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = fill_values) +
  scale_linetype_manual(values = linetype_values) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = seq(2014, 2025, by = 2), limits = c(primeiro_ano, 2025)) +
  scale_y_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Year", y = "Jaccard similarity index") +
  ggtitle("(d)") +
  theme_panel

# ================================================================================
# 7. COMBINE PANELS WITH LEGEND
# ================================================================================

cat("Combining panels with legend...\n")

# Create a dummy plot to extract legend
legend_plot <- ggplot() +
  geom_point(data = data.frame(x = 1:2, y = 1:2,
                               Grupo = c("Abandoned pastures", "Natural grassland")),
             aes(x = x, y = y, color = Grupo, fill = Grupo, shape = Grupo), size = 3) +
  geom_line(data = data.frame(x = c(1, 1.5, 2, 2.5), y = c(1, 1.5, 2, 2.5),
                              Grupo = rep(c("Abandoned pastures", "Natural grassland"), each = 2)),
            aes(x = x, y = y, color = Grupo, linetype = Grupo), linewidth = 1.2) +
  scale_color_manual(values = color_values, name = NULL) +
  scale_fill_manual(values = fill_values, name = NULL) +
  scale_linetype_manual(values = linetype_values, name = NULL) +
  scale_shape_manual(values = shape_values, name = NULL) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.5, "cm"),
    legend.direction = "horizontal"
  )

# Extract legend
legend <- get_legend(legend_plot)

# Combine using patchwork
combined_figure <- (panel_a | panel_b) / (panel_c | panel_d) / legend +
  plot_layout(heights = c(1, 1, 0.1))

# Save
ggsave("outputs/figures/Figure_Temporal_Dynamics.tiff",
       combined_figure,
       width = 18, height = 16, units = "cm", dpi = 300, compression = "lzw")

cat("\n================================================================================\n")
cat("PANEL FIGURE COMPLETED SUCCESSFULLY!\n")
cat("Saved as: outputs/figures/Figure_Temporal_Dynamics.tiff\n")
cat("\nLegend nomenclature:\n")
cat("- Abandoned pastures: Former cattle pastures under passive restoration since 2007-2008\n")
cat("- Natural grassland: Continuously grassland, cattle exclusion since 2007-2008\n")
cat("================================================================================\n")

