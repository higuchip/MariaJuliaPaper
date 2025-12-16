# ==============================================================================
# FOREST INVENTORY DATA PROCESSING
# ==============================================================================
#
# Description: Processes long-term forest inventory data from permanent plots,
#              identifying mortality (M) and addition (AD) events, and preparing
#              data for population dynamics analyses.
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
# ==============================================================================

# --- Setup --------------------------------------------------------------------

# Required packages
library(dplyr)
library(tidyr)

# --- 1. Data Loading ----------------------------------------------------------

#' Load raw forest inventory data
#' 
#' Data format: CSV with semicolon separator and comma as decimal mark
#' Expected columns: area, parc, ni, especies, hYYYY (height measurements)
#' Special values in height columns:
#'   - "M"  = Mortality (individual died)
#'   - "AD" = Adult (individual reached adult size, i.e., DBH ≥ 5cm)

dados_orig <- read.csv2(
  "dados/dados_2025_atualizado.csv", 
  fileEncoding = "latin1", 
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "nan", "")
)

# --- 2. Height Column Processing ----------------------------------------------

#' Identify height columns
#' Format: hYYYY (e.g., h2014, h2015, h2017, h2019, h2021, h2025)
#' Note: Irregular sampling intervals due to field logistics

height_cols <- grep("^h\\d{4}$", names(dados_orig), value = TRUE)

#' Create indicator columns for demographic events
#' 
#' For each height column, create:
#'   - m_YYYY:  TRUE if individual died (recorded as "M")
#'   - ad_YYYY: TRUE if individual reached adult size (recorded as "AD")

for (col in height_cols) {
  # Define new column names (h2014 -> m_2014, ad_2014)
  m_col_name <- paste0("m_", sub("h", "", col))
  ad_col_name <- paste0("ad_", sub("h", "", col))
  
  # Clean and standardize values (trim spaces, convert to uppercase)
  original_col_cleaned <- trimws(toupper(dados_orig[[col]]))
  
  # Create boolean indicator columns
  dados_orig <- dados_orig %>%
    mutate(
      !!m_col_name := !is.na(original_col_cleaned) & original_col_cleaned == "M",
      !!ad_col_name := !is.na(original_col_cleaned) & original_col_cleaned == "AD"
    )
}

# --- 3. Numeric Data Conversion -----------------------------------------------

#' Convert height columns to numeric
#' Values "M" and "AD" become NA (already captured in indicator columns)

for (col in height_cols) {
  # Replace decimal comma with period
  dados_orig[[col]] <- gsub(",", ".", dados_orig[[col]])
  
  # Convert to numeric (suppressWarnings to ignore NA warnings)
  dados_orig[[col]] <- suppressWarnings(as.numeric(dados_orig[[col]]))
}

#' Convert spatial coordinates to numeric (if present)

if ("x" %in% names(dados_orig)) {
  dados_orig$x <- gsub(",", ".", dados_orig$x)
  dados_orig$x <- suppressWarnings(as.numeric(dados_orig$x))
}

if ("y" %in% names(dados_orig)) {
  dados_orig$y <- gsub(",", ".", dados_orig$y)
  dados_orig$y <- suppressWarnings(as.numeric(dados_orig$y))
}

# --- 4. Data Validation -------------------------------------------------------

#' Basic validation checks
#' Uncomment to run diagnostics

# cat("Number of individuals:", nrow(dados_orig), "\n")
# cat("Number of sites:", length(unique(dados_orig$area)), "\n")
# cat("Number of plots:", length(unique(paste(dados_orig$area, dados_orig$parc))), "\n")
# cat("Height columns processed:", paste(height_cols, collapse = ", "), "\n")

# --- 5. View Processed Data ---------------------------------------------------

View(dados_orig)

# ==============================================================================
# OUTPUT STRUCTURE
# ==============================================================================
#
# The processed dataframe (dados_orig) contains:
#
# Identification columns:
#   - area:     Site identifier (1 = abandoned pasture site A, etc.)
#   - parc:     Plot number within site (1-20)
#   - ni:       Individual identifier within plot
#   - especies: Species name
#
# Height measurements (numeric):
#   - h2014, h2015, h2017, h2019, h2021, h2025
#   - Values in meters; NA indicates no measurement (dead or not yet recruited)
#
# Demographic indicators (logical):
#   - m_YYYY:  TRUE if individual died in/before year YYYY
#   - ad_YYYY: TRUE if individual reached adult size in year YYYY
#
# Spatial coordinates (if present):
#   - x, y: Position within plot (meters)
#
# ==============================================================================
# END OF SCRIPT
# ==============================================================================
