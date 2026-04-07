---
title: "Sperm Competition Analysis Pipeline"
---

# Sperm Competition Analysis Pipeline

This script processes sperm competition datasets, calculates the proportion of offspring sired by the second male, and generates publication-quality plots with statistical comparisons for both XO vs XY and XgY vs XdY experiments.

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
      male_total        = m_dark + m_genome,
      female_total      = f_dark + f_genome,
      total_offspring   = male_total + female_total
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
# Load + combine datasets
# =========================
files <- tibble::tribble(
  ~path, ~dataset,
  "data/SpermCompetitionPilot.csv",  "Pilot",
  "data/SpermCompetitionRound1.csv", "Round 1",
  "data/SpermCompetitionRound2.csv", "Round 2"
)

sperm_summary_clean_all <- purrr::map2_dfr(files$path, files$dataset, process_file)

# =========================
# Relabel datasets → Blocks
# =========================
sperm_summary_clean_all <- sperm_summary_clean_all %>%
  mutate(
    dataset = recode(dataset,
                     "Pilot"   = "Block 1",
                     "Round 1" = "Block 2",
                     "Round 2" = "Block 3"),
    second_mate = factor(second_mate, levels = c("XO", "XY"))
  )

# =========================
# Filter
# =========================
plot_df_main <- sperm_summary_clean_all %>%
  filter(prop_second_male > 0)

# =========================
# Stats: XO vs XY
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
              size = 3,
              shape = 16) +
  
  annotate("segment", x = 1, xend = 2, y = 1.08, yend = 1.08) +
  annotate("segment", x = 1, xend = 1, y = 1.06, yend = 1.08) +
  annotate("segment", x = 2, xend = 2, y = 1.06, yend = 1.08) +
  annotate("text", x = 1.5, y = 1.095, label = p_star_main, size = 5) +
  
  scale_color_manual(values = c(
    "Block 1" = "#1b9e77",
    "Block 2" = "#d95f02",
    "Block 3" = "#7570b3"
  )) +
  
  labs(
    title = "P1–P2 Sperm Competition",
    x = "Second Mate Type",
    y = "Proportion Females Sired by Second Male",
    color = "Block"
  ) +
  
  coord_cartesian(ylim = c(0, 1.12)) +
  theme_minimal(base_size = 25)

# =========================
# XgY vs XdY processing
# =========================
process_xgy <- function(path) {
  
  df <- read.csv(path, header = TRUE) %>%
    select(1:7)
  
  names(df) <- c("female", "second_mate", "date",
                 "m_dark", "f_dark", "m_genome", "f_genome")
  
  df <- df %>%
    mutate(
      female_total = f_dark + f_genome
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
      prop_second_male = ifelse(second_mate == "XdY",
                                female_dark_total / female_total,
                                female_genome_total / female_total)
    ) %>%
    filter(
      !is.na(prop_second_male),
      is.finite(prop_second_male),
      female_total > 5,
      !(second_mate == "XgY" & female_genome_total == 0)
    )
  
  return(summary)
}

xgy1 <- process_xgy("data/XgYvsXdY.csv")
xgy2 <- process_xgy("data/XgYvsXdY_block2.csv")

plot_df_xgy <- bind_rows(xgy1, xgy2) %>%
  filter(prop_second_male > 0)

# =========================
# Stats: XgY vs XdY
# =========================
pval_xgy <- with(plot_df_xgy, t.test(prop_second_male ~ second_mate)$p.value)

p_star_xgy <- ifelse(pval_xgy < 0.001, "***",
                     ifelse(pval_xgy < 0.01, "**",
                            ifelse(pval_xgy < 0.05, "*", "ns")))

# =========================
# Plot: XgY vs XdY
# =========================
ggplot(plot_df_xgy,
       aes(x = second_mate,
           y = prop_second_male,
           color = second_mate)) +
  
  geom_hline(yintercept = 0.5, linetype = "dotted", linewidth = 0.5) +
  
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = NA, size = 1) +
  
  geom_jitter(width = 0.12, alpha = 0.6, size = 2.6) +
  
  annotate("segment", x = 1, xend = 2, y = 1.08, yend = 1.08) +
  annotate("segment", x = 1, xend = 1, y = 1.06, yend = 1.08) +
  annotate("segment", x = 2, xend = 2, y = 1.06, yend = 1.08) +
  annotate("text", x = 1.5, y = 1.095, label = p_star_xgy, size = 5) +
  
  scale_color_manual(values = c(
    "XdY" = "#E69F00",
    "XgY" = "#56B4E9"
  ), name = "Genotype") +
  
  scale_x_discrete(limits = c("XdY", "XgY")) +
  
  coord_cartesian(ylim = c(0, 1.12)) +
  
  labs(
    title = "P1–P2 Sperm Competition",
    x = "Second Mate Type",
    y = "Proportion Sired by Second Male"
  ) +
  
  theme_minimal(base_size = 25)
```
