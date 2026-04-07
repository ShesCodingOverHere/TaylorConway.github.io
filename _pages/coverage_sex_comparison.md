---
layout: single
title: "Sex-Specific Coverage Analysis (Short and Long Reads)"
permalink: /code/coverage_sex_comparison/
author_profile: true
---

# Sex-Specific Coverage Analysis (Short and Long Reads)

This script compares sequencing coverage between male and female samples across genomic scaffolds using both short-read and long-read datasets. Coverage differences are summarized using log2 ratios and visualized to identify sex-linked regions (e.g., Y-linked scaffolds).

The workflow standardizes coverage calculations, removes low-quality regions, and highlights candidate sex-specific scaffolds.

---

## Script

```r
library(tidyverse)
library(ggplot2)
library(cowplot)

# =========================
# Helper function
# =========================
process_coverage <- function(male_path, female_path) {
  
  male <- read.table(male_path)
  female <- read.table(female_path)
  
  names(male)   <- paste("m.", c("chrom","start","end","reads","cov"), sep = "")
  names(female) <- paste("f.", c("chrom","start","end","reads","cov"), sep = "")
  
  df <- cbind(male, female) %>%
    filter(!m.chrom %in% c("Unknown", "Y_unplaced")) %>%
    filter(m.cov + f.cov > 0) %>%
    mutate(
      log2_ratio = log2((m.cov + 1) / (f.cov + 1)),
      highlight = ifelse(m.chrom == "Y_candidate", "Y_candidate", "other")
    )
  
  return(df)
}

# =========================
# Chromosome labels
# =========================
chrom_labels <- c(
  "ChrX" = "Chr X",
  "Chr2" = "Chr 2",
  "Chr3" = "Chr 3",
  "Chr4" = "Chr 4",
  "Chr5" = "Chr 5",
  "Y_candidate" = "Y-linked scaffold"
)

# =========================
# Short-read analysis
# =========================
short_df <- process_coverage(
  "data/input/male_shortreads.cov",
  "data/input/female_shortreads.cov"
)

ggplot(short_df, aes(x = m.chrom, y = log2_ratio, color = highlight)) +
  geom_jitter(width = 0.3, size = 1.5, alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
  geom_hline(yintercept = c(1, -1, 0), linetype = "dashed") +
  scale_color_manual(values = c("Y_candidate" = "#D55E00", "other" = "grey60")) +
  coord_cartesian(ylim = c(-10, 10)) +
  scale_x_discrete(labels = chrom_labels) +
  labs(
    x = "Chromosome",
    y = expression(log[2]~"(Male / Female)"),
    title = "Short-read coverage by chromosome"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# =========================
# Long-read analysis
# =========================
long_df <- process_coverage(
  "data/input/male_longreads.cov",
  "data/input/female_longreads.cov"
)

ggplot(long_df, aes(x = m.chrom, y = log2_ratio, color = highlight)) +
  geom_jitter(width = 0.3, size = 1.5, alpha = 0.3) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
  geom_hline(yintercept = c(1, -1, 0), linetype = "dashed") +
  scale_color_manual(values = c("Y_candidate" = "#0072B2", "other" = "grey60")) +
  coord_cartesian(ylim = c(-10, 10)) +
  scale_x_discrete(labels = chrom_labels) +
  labs(
    x = "Chromosome",
    y = expression(log[2]~"(Male / Female)"),
    title = "Long-read coverage by chromosome"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
```

---

## Notes

- Coverage is normalized using log2(Male / Female) ratios  
- A pseudocount (+1) is added to avoid division by zero  
- Low-quality scaffolds and zero-coverage regions are removed  
- Candidate sex-linked scaffolds can be highlighted for interpretation  
- Both short-read and long-read datasets are processed using the same workflow  

---

## Key Insight

Regions with strongly biased coverage between sexes (e.g., consistently high or low log2 ratios) are candidates for sex-linked genomic regions, such as Y-linked scaffolds.
