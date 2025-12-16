# Vegetation dynamics in abandoned Atlantic Forest highland pastures

Data and R scripts for:

> **Cruz MJC et al.** (in review) Forest recovery or shrubland assembly? Vegetation dynamics in abandoned Atlantic Forest highland pastures. *Biotropica*.

## Repository Structure

```
├── dados/
│   ├── dados_2025_atualizado.csv                   # Permanent plot data (2014-2025)
│   └── HIGUCHI-SILVA LABDENDRO DATABASE  2025.csv  # Reference forest data
│
├── scripts/
│   ├── 01_data_processing.R
│   ├── 02_study_area_map.R
│   ├── 03_floristic_similarity.R
│   ├── 04_abundance_dynamics.R
│   ├── 05_height_dynamics.R
│   ├── 06_richness_dynamics.R
│   ├── 07_panel_figure.R
│   ├── 08_compositional_analysis.R
│   └── 09_survival_analysis.R
│
└── README.md
```

## Usage

Run scripts sequentially (01 → 09). Each script depends on objects from previous scripts.

## Contact

Pedro Higuchi - higuchip@gmail.com
