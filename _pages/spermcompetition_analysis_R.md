---
layout: single
title: "Sperm Competition Analysis"
permalink: /code/sperm_competition/
author_profile: true
---

# Sperm Competition Analysis

This script processes sperm competition datasets to quantify the proportion of offspring sired by the second male (P2). It integrates multiple experimental blocks, computes summary statistics, and generates publication-quality visualizations comparing mating outcomes across genotypes.

The workflow supports two experimental designs:
- XO vs XY male competition
- XgY vs XdY male competition

---

## Script

```r
# =========================
# Libraries
# =========================
library(tidyverse)
library(ggplot2)

# =========================
# Helper: process XO vs XY datasets
# =========================
process_file <- function(path, dataset_label) {
  
  df <- read.csv(path, header = TRUE) %>% 
    select(1:7)
  
  names(df) <- c("female", "second_mate", "date",
                 "m_dark", "f_dark", "m_genome", "f_genome")
  
  df <- df %>%
    mutate(
      male_total      = m_dark + m_genome,
      female_total    = f_dark + f_genome,
      total_offspring = male_total + female_total
    )
  
  summary <- df %>%
    group_by(female, second_mate) %>%
    summarize(
      female_dark_total   = sum(f_dark, na.rm = TRUE),
      female_genome_total = sum(f_genome, na.rm = TRUE),
      female_total        = sum(f_dark + f_genome, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      prop_second_male = ifelse(second_mate == "XO",
                                female_dark_total / female_total,
                                female_genome_total / female_total),
      dataset = dataset_label
    ) %>%
    filter(
      !is.na(prop_second_male),
      is.finite(prop_second_male),
      female_total > 0,
      !(second_mate == "XY" & female_genome_total == 0)
    )
  
  return(summary)
}

# =========================
# Load and combine datasets
# =========================
files <- tibble::tribble(
  ~path, ~dataset,
  "data/input/block1.csv", "Block 1",
  "data/input/block2.csv", "Block 2",
  "data/input/block3.csv", "Block 3"
)

sperm_summary <- purrr::map2_dfr(files$path, files$dataset, process_file)

sperm_summary <- sperm_summary %>%
  mutate(second_mate = factor(second_mate, levels = c("XO", "XY")))

plot_df_main <- sperm_summary %>%
  filter(prop_second_male > 0)

# =========================
# Statistical test
# =========================
pval_main <- with(plot_df_main, t.test(prop_second_male ~ second_mate)$p.value)

p_star_main <- ifelse(pval_main < 0.001, "***",
                      ifelse(pval_main < 0.01, "**",
                             ifelse(pval_main < 0.05, "*", "ns")))

# =========================
# Plot: XO vs XY
# =========================
ggplot(plot_df_main,
       aes(x = second_mate,
           y = prop_second_male)) +
  
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4) +
  
  geom_boxplot(aes(color = second_mate),
               outlier.shape = NA,
               width = 0.55,
               fill = NA,
               size = 1) +
  
  geom_jitter(aes(color = dataset),
              width = 0.25,
              alpha = 0.8,
              size = 3) +
  
  annotate("text", x = 1.5, y = 1.05, label = p_star_main, size = 5) +
  
  labs(
    title = "Sperm Competition (XO vs XY)",
    x = "Second Mate Type",
    y = "Proportion Sired by Second Male",
    color = "Block"
  ) +
  
  coord_cartesian(ylim = c(0, 1.1)) +
  theme_minimal(base_size = 18)

# =========================
# XgY vs XdY analysis
# =========================
process_xgy <- function(path) {
  
  df <- read.csv(path, header = TRUE) %>%
    select(1:7)
  
  names(df) <- c("female", "second_mate", "date",
                 "m_dark", "f_dark", "m_genome", "f_genome")
  
  df <- df %>%
    mutate(female_total = f_dark + f_genome)
  
  summary <- df %>%
    group_by(female, second_mate) %>%
    summarize(
      female_dark_total   = sum(f_dark, na.rm = TRUE),
      female_genome_total = sum(f_genome, na.rm = TRUE),
      female_total        = sum(f_dark + f_genome, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      prop_second_male = ifelse(second_mate == "XdY",
                                female_dark_total / female_total,
                                female_genome_total / female_total)
    ) %>%
    filter(
      !is.na(prop_second_male),
      is.finite(prop_second_male),
      female_total > 5
    )
  
  return(summary)
}

xgy_df <- bind_rows(
  process_xgy("data/input/xgy_block1.csv"),
  process_xgy("data/input/xgy_block2.csv")
)

pval_xgy <- with(xgy_df, t.test(prop_second_male ~ second_mate)$p.value)

p_star_xgy <- ifelse(pval_xgy < 0.001, "***",
                     ifelse(pval_xgy < 0.01, "**",
                            ifelse(pval_xgy < 0.05, "*", "ns")))

# =========================
# Plot: XgY vs XdY
# =========================
ggplot(xgy_df,
       aes(x = second_mate,
           y = prop_second_male,
           color = second_mate)) +
  
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
  
  annotate("text", x = 1.5, y = 1.05, label = p_star_xgy, size = 5) +
  
  labs(
    title = "Sperm Competition (XgY vs XdY)",
    x = "Second Mate Type",
    y = "Proportion Sired by Second Male"
  ) +
  
  coord_cartesian(ylim = c(0, 1.1)) +
  theme_minimal(base_size = 18)
```

---

## Notes

- P2 represents the proportion of offspring sired by the second male  
- Data are aggregated per female to account for replicate matings  
- Statistical significance is assessed using a two-sample t-test  
- Visualization includes both distribution (boxplot) and individual observations (jitter)  
- Designed for comparing mating success across genetic backgrounds  
