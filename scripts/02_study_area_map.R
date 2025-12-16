# ==============================================================================
# STUDY AREA LOCATION MAP
# ==============================================================================
#
# Description: Creates a multi-panel location map showing the study area at
#              three spatial scales: South America/Brazil, Santa Catarina state,
#              and Urubici municipality with sampling sites.
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
# Output:      ../artigo/biotropica/review/Fig1.tiff (180 x 190 mm, 300 dpi)
#              ../artigo/biotropica/review/study_areas_coordinates.csv
#
# Note:        R project located in MariaJulia/dissertacao
#
# ==============================================================================

# --- Setup --------------------------------------------------------------------

# Required packages
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(geobr)
library(patchwork)
library(ggspatial)
library(dplyr)
library(cowplot)

#' Custom theme for scientific publication
#' 
#' Based on theme_bw with modifications for clean, publication-ready maps

theme_publication <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      # Panels
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", size = 0.5),
      panel.grid.major = element_line(color = "gray80", size = 0.2),
      panel.grid.minor = element_blank(),
      
      # Axes
      axis.text = element_text(color = "black", size = rel(0.8)),
      axis.title = element_text(color = "black", size = rel(0.9)),
      axis.ticks = element_line(color = "black", size = 0.3),
      
      # Legend
      legend.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size = rel(0.8)),
      legend.title = element_text(size = rel(0.9), face = "bold"),
      
      # Margins
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    )
}

# --- 1. Prepare Coordinate Data -----------------------------------------------

#' Convert degrees-minutes-seconds to decimal degrees
dms_to_decimal <- function(degrees, minutes, seconds, direction) {
  decimal <- degrees + minutes/60 + seconds/3600
  if (direction %in% c("S", "W")) decimal <- -decimal
  return(decimal)
}

#' Study area coordinates
#' Areas 1-3: Abandoned pastures (former Araucaria Forest sites)
#' Area 4: Natural grassland (reference site, burned in 2018)

study_areas <- data.frame(
  area = c("Area 1", "Area 2", "Area 3", "Area 4"),
  lat = c(
    dms_to_decimal(28, 05, 41.5, "S"),
    dms_to_decimal(28, 04, 46.87, "S"),
    dms_to_decimal(28, 09, 49.19, "S"),
    dms_to_decimal(28, 08, 35, "S")
  ),
  lon = c(
    dms_to_decimal(49, 30, 14.71, "W"),
    dms_to_decimal(49, 30, 51.29, "W"),
    dms_to_decimal(49, 36, 47.56, "W"),
    dms_to_decimal(49, 38, 09, "W")
  ),
  altitude = c(1628, 1356, 1660, 1377)
)

# Convert to sf object
study_areas_sf <- st_as_sf(study_areas, coords = c("lon", "lat"), crs = 4326)

# --- 2. Load Geographic Data --------------------------------------------------

print("Loading geographic data...")

# South America and Brazil
world <- ne_countries(scale = "medium", returnclass = "sf")
south_america <- world[world$continent == "South America", ]
brazil <- world[world$name == "Brazil", ]

# Brazilian states
estados <- read_state(year = 2020)
sc <- estados[estados$abbrev_state == "SC", ]

# SC municipalities (to get Urubici)
municipios_sc <- read_municipality(code_muni = "SC", year = 2020)
urubici <- municipios_sc[municipios_sc$name_muni == "Urubici", ]

# Fallback if Urubici not found
if(nrow(urubici) == 0) {
  print("Creating approximate Urubici area...")
  urubici_center <- st_sfc(st_point(c(mean(study_areas$lon), mean(study_areas$lat))), crs = 4326)
  urubici <- st_buffer(urubici_center, dist = 0.15)
  urubici <- st_sf(geometry = urubici)
}

# --- 3. Panel A: South America/Brazil -----------------------------------------

print("Creating Panel A - Brazil...")

panel_a <- ggplot() +
  # South America (background)
  geom_sf(data = south_america, fill = "white", color = "black", size = 0.3) +
  
  # Brazil in gray
  geom_sf(data = brazil, fill = "gray80", color = "black", size = 0.5) +
  
  # Santa Catarina highlighted
  geom_sf(data = sc, fill = "black", color = "black", size = 0.5) +
  
  # Limits and coordinates
  coord_sf(xlim = c(-74, -34), ylim = c(-33, 5), expand = FALSE) +
  scale_x_continuous(breaks = seq(-70, -40, by = 10)) +
  scale_y_continuous(breaks = seq(-30, 0, by = 10)) +
  labs(x = NULL, y = NULL) +
  theme_publication() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

# --- 4. Panel B: Santa Catarina State -----------------------------------------

print("Creating Panel B - Santa Catarina...")

panel_b <- ggplot() +
  # Santa Catarina
  geom_sf(data = sc, fill = "gray90", color = "black", size = 0.5) +
  
  # Urubici highlighted
  geom_sf(data = urubici, fill = "gray50", color = "black", size = 0.5) +
  
  # Labels
  annotate("text", x = -49.5, y = -28.0, label = "Urubici", 
           size = 3.5, fontface = "bold") +
  annotate("text", x = -50.5, y = -26.5, label = "Santa Catarina", 
           size = 4, fontface = "bold") +
  
  # Limits
  coord_sf(xlim = c(-54, -48), ylim = c(-29.5, -26), expand = FALSE) +
  scale_x_continuous(breaks = seq(-54, -48, by = 2)) +
  scale_y_continuous(breaks = seq(-29, -26, by = 1)) +
  labs(x = NULL, y = NULL) +
  theme_publication() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

# --- 5. Panel C: Urubici Detail with Study Sites ------------------------------

print("Creating Panel C - Urubici Detail...")

# Calculate bounding box
if(class(urubici)[1] == "sf") {
  bbox_urubici <- st_bbox(urubici)
  bbox_points <- st_bbox(study_areas_sf)
  
  xlim_c <- c(min(bbox_urubici[1], bbox_points[1]) - 0.05, 
              max(bbox_urubici[3], bbox_points[3]) + 0.05)
  ylim_c <- c(min(bbox_urubici[2], bbox_points[2]) - 0.05, 
              max(bbox_urubici[4], bbox_points[4]) + 0.05)
} else {
  xlim_c <- c(min(study_areas$lon) - 0.15, max(study_areas$lon) + 0.15)
  ylim_c <- c(min(study_areas$lat) - 0.1, max(study_areas$lat) + 0.1)
}

panel_c <- ggplot() +
  theme(panel.background = element_rect(fill = "gray95")) +
  
  # Urubici municipality
  geom_sf(data = urubici, fill = "gray70", color = NA, alpha = 0.8) +
  geom_sf(data = urubici, fill = NA, color = "black", size = 0.8, linetype = "solid") +
  
  # Study area points
  geom_sf(data = study_areas_sf, shape = 21, 
          fill = "black", color = "white", size = 3, stroke = 1) +
  
  # Area labels
  geom_text(data = study_areas,
            aes(x = lon, y = lat, label = area),
            hjust = -0.3, vjust = 0.5, size = 3, fontface = "bold") +
  
  # Municipality label
  annotate("text", 
           x = mean(study_areas$lon), 
           y = max(study_areas$lat) + 0.03,
           label = "URUBICI", 
           size = 4, fontface = "bold", color = "gray30") +
  
  # Map limits
  coord_sf(xlim = xlim_c, ylim = ylim_c, expand = FALSE) +
  
  # Scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.7,
    line_width = 0.8,
    height = unit(0.15, "cm"),
    pad_x = unit(0.3, "cm"),
    pad_y = unit(0.3, "cm"),
    style = "bar",
    bar_cols = c("black", "white")
  ) +
  
  # North arrow
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(0.8, "cm"),
    width = unit(0.6, "cm"),
    pad_x = unit(0.2, "cm"),
    pad_y = unit(0.2, "cm"),
    style = north_arrow_orienteering(
      fill = c("black", "white"),
      line_col = "black",
      text_size = 8
    )
  ) +
  
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = seq(-50, -49, by = 0.25)) +
  scale_y_continuous(breaks = seq(-28.5, -27.5, by = 0.25)) +
  theme_publication() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray90", size = 0.2)
  )

# --- 6. Combine Panels and Save -----------------------------------------------

print("Creating final layout...")

# Combine panels: A and B on top, C on bottom
mapa_final <- (panel_a | panel_b) / panel_c +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    caption = "Datum: SIRGAS 2000",
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5),
      plot.caption = element_text(size = 8, hjust = 1, margin = margin(t = 5))
    )
  )

# Save figure
print("Saving figure...")

# Output path (relative to R project in MariaJulia/dissertacao)
output_path <- "../artigo/biotropica/review/"

ggsave(
  filename = paste0(output_path, "Fig1.tiff"),
  plot = mapa_final,
  width = 180,
  height = 190,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)

print("Map created successfully!")
print(paste0("File saved: ", output_path, "Fig1.tiff"))

# --- 7. Export Coordinate Table -----------------------------------------------

coord_table <- study_areas %>%
  mutate(
    latitude = sprintf("%.2f°S", abs(lat)),
    longitude = sprintf("%.2f°W", abs(lon)),
    altitude_m = altitude
  ) %>%
  select(area, latitude, longitude, altitude_m)

print("\nCoordinate table for manuscript:")
print(coord_table)

write.csv(coord_table, paste0(output_path, "study_areas_coordinates.csv"), row.names = FALSE)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================