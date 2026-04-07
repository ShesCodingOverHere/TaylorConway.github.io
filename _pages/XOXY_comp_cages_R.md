---
layout: single
title: "XO/XY Competition Cage Analysis"
permalink: /code/xoxy_cage_competition/
author_profile: true
---

# XO/XY Competition Cage Analysis

This script analyzes the frequency of XO males across experimental population cages over time. It summarizes genotype frequencies, applies filtering based on marker presence, and visualizes temporal trends across replicate cages.

The workflow calculates the proportion of XO males within each cage and collection timepoint, and tracks changes over time using both categorical collection points and continuous time (days).

---

## Script

```r
library(tidyverse)
library(lubridate)

# =========================
# Load data (generic path)
# =========================
df <- read.csv("data/input/cage_metadata.csv")

# =========================
# Initial summary by collection
# =========================
summary_df <- df %>%
  group_by(Cage, Collection) %>%
  summarise(
    freq_XO = mean(Y_primer == 0),
    .groups = "drop"
  )

# Plot: frequency over collections
ggplot(summary_df, aes(x = Collection, y = freq_XO, color = factor(Cage))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Collection",
    y = "Frequency of XO males",
    color = "Cage"
  ) +
  theme_classic(base_size = 14)

# Smoothed trend
ggplot(summary_df, aes(x = Collection, y = freq_XO, color = factor(Cage))) +
  geom_line(alpha = 0.6) +
  geom_smooth(se = FALSE, method = "loess") +
  theme_classic()

# =========================
# Filter: only X-positive individuals
# =========================
df_filtered <- df %>%
  filter(X_primer == 1)

summary_df <- df_filtered %>%
  group_by(Cage, Collection) %>%
  summarise(
    freq_XO = mean(Y_primer == 0),
    n = n(),
    .groups = "drop"
  )

# Plot with sample size
ggplot(summary_df, aes(x = Collection, y = freq_XO, color = factor(Cage))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Collection",
    y = "Frequency of XO males",
    color = "Cage"
  ) +
  theme_classic(base_size = 14)

ggplot(summary_df, aes(x = Collection, y = freq_XO,
                       color = factor(Cage), size = n)) +
  geom_line() +
  geom_point() +
  theme_classic()

# =========================
# Convert collection dates to time scale
# =========================
df$Collection_date <- mdy(df$Collection.date)

start_date <- mdy("01/01/2025")  # generic start date

df$days <- as.numeric(df$Collection_date - start_date)

df_filtered <- df %>%
  filter(X_primer == 1)

summary_df <- df_filtered %>%
  group_by(Cage, days) %>%
  summarise(
    freq_XO = mean(Y_primer == 0),
    n = n(),
    .groups = "drop"
  )

# Final plot: frequency over time
ggplot(summary_df, aes(x = days, y = freq_XO, color = factor(Cage))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    x = "Days since start",
    y = "Frequency of XO males",
    color = "Cage"
  ) +
  theme_classic(base_size = 25)
```

---

## Notes

- XO frequency is calculated as the proportion of individuals lacking the Y marker  
- Filtering for `X_primer == 1` ensures only valid individuals are included  
- Temporal trends are visualized both by discrete collection events and continuous time  
- Sample size (`n`) can be incorporated into visualization to assess data robustness  
- Designed for tracking genotype frequency changes across experimental populations  
