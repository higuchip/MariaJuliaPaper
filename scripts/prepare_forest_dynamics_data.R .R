# ===============================================================
# FOREST INVENTORY DATA PROCESSING
# ===============================================================
# Description: Script to process forest inventory data,
# identifying mortality (M), addition (AD) events, and preparing 
# data for population dynamics analyses.
#
# Author: Pedro Higuchi
# Date: 2025
# ===============================================================

# --- Load Required Packages ---
library(dplyr)
library(tidyr)

# ===============================================================
# 1. DATA LOADING
# ===============================================================

# Load forest inventory data
# Format: CSV with ';' separator and ',' as decimal mark
dados_orig <- read.csv2(
  "dados/dados_2025_atualizado.csv", 
  fileEncoding = "latin1", 
  stringsAsFactors = FALSE, 
  na.strings = c("NA", "nan", "")
)

# Create copy for processing (preserve original data)
df_processed <- dados_orig

# ===============================================================
# 2. HEIGHT COLUMN PROCESSING
# ===============================================================

# Identify height columns (format: hYYYY, e.g., h2014, h2015)
height_cols <- grep("^h\\d{4}$", names(df_processed), value = TRUE)

# Create indicator columns for mortality (M) and addition (AD)
# NOTE: In the original file, cells with "M" indicate mortality
#       and cells with "AD" indicate adult individual additions

for (col in height_cols) {
  # Define new column names
  m_col_name <- paste0("m_", sub("h", "", col))   # e.g., h2014 -> m_2014
  ad_col_name <- paste0("ad_", sub("h", "", col)) # e.g., h2014 -> ad_2014
  
  # Clean and standardize values (trim spaces, convert to uppercase)
  original_col_cleaned <- trimws(toupper(df_processed[[col]]))
  
  # Create boolean columns for mortality and addition
  df_processed <- df_processed %>%
    mutate(
      # Mortality column: TRUE if cell contains "M"
      !!m_col_name := ifelse(
        !is.na(original_col_cleaned) & original_col_cleaned == "M", 
        TRUE, 
        FALSE
      ),
      # Addition column: TRUE if cell contains "AD"
      !!ad_col_name := ifelse(
        !is.na(original_col_cleaned) & original_col_cleaned == "AD", 
        TRUE, 
        FALSE
      )
    )
}

# ===============================================================
# 3. NUMERIC DATA CONVERSION
# ===============================================================

# Convert height columns to numeric
# Values "M" and "AD" will be converted to NA
for (col in height_cols) {
  # Replace decimal comma with period
  df_processed[[col]] <- gsub(",", ".", df_processed[[col]])
  
  # Convert to numeric (suppressWarnings to ignore NA warnings)
  df_processed[[col]] <- suppressWarnings(as.numeric(df_processed[[col]]))
}

# Convert spatial coordinates (x, y) to numeric, if they exist
if ("x" %in% names(df_processed)) {
  df_processed$x <- gsub(",", ".", df_processed$x)
  df_processed$x <- suppressWarnings(as.numeric(df_processed$x))
}

if ("y" %in% names(df_processed)) {
  df_processed$y <- gsub(",", ".", df_processed$y)
  df_processed$y <- suppressWarnings(as.numeric(df_processed$y))
}

# ===============================================================
# 4. VIEW PROCESSED DATA
# ===============================================================

# View processed data
View(df_processed)

# ===============================================================
# FINAL DATAFRAME STRUCTURE
# ===============================================================
# The processed dataframe contains:
# - Original identification columns (area, parc, ni, especies)
# - Height columns converted to numeric (hYYYY)
# - Mortality columns (m_YYYY): TRUE/FALSE
# - Addition columns (ad_YYYY): TRUE/FALSE
# - Spatial coordinates x, y in numeric format
#
# This format facilitates analyses of:
# - Population dynamics
# - Mortality and recruitment rates
# - Spatial patterns
# ===============================================================
