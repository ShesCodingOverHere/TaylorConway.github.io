---
title: "NUMTs Detection and Analysis Pipeline"
---

# NUMTs Detection and Analysis Pipeline

This pipeline describes a complete workflow for identifying nuclear mitochondrial DNA segments (NUMTs) from whole genome sequencing data. The approach integrates read preprocessing, genome alignment, sequence similarity searches, and filtering to generate high-confidence NUMT candidates. It is designed to be reproducible, portable, and free of machine-specific paths.

The workflow starts with raw sequencing reads and reference genomes. Inputs include paired-end FASTQ files, a nuclear genome FASTA, and a mitochondrial genome FASTA. Outputs include alignment files (BAM), candidate NUMT regions (BED), and optional summary statistics.

A recommended project structure:

```
project/
├── data/
│   ├── raw_reads/
│   ├── reference/
│   └── intermediate/
├── results/
│   ├── alignments/
│   ├── numt_candidates/
│   └── plots/
├── scripts/
└── logs/
```

Required software: fastqc, trimmomatic (or equivalent), bwa or bowtie2, samtools, bedtools, BLAST+ (blastn), and R (tidyverse, ggplot2).

Raw read quality control:

```
fastqc data/raw_reads/*.fastq.gz -o results/plots/
```

Adapter and quality trimming:

```
trimmomatic PE \
  input_R1.fastq.gz input_R2.fastq.gz \
  trimmed_R1.fastq.gz unpaired_R1.fastq.gz \
  trimmed_R2.fastq.gz unpaired_R2.fastq.gz \
  ILLUMINACLIP:adapters.fa:2:30:10 \
  SLIDINGWINDOW:4:20 \
  MINLEN:50
```

Index reference genomes:

```
bwa index data/reference/nuclear_genome.fa
bwa index data/reference/mitochondrial_genome.fa
```

Align reads to the nuclear genome:

```
bwa mem data/reference/nuclear_genome.fa \
  trimmed_R1.fastq.gz trimmed_R2.fastq.gz | \
  samtools sort -o results/alignments/nuclear.bam

samtools index results/alignments/nuclear.bam
```

Extract unmapped reads that may contain mitochondrial-origin sequence:

```
samtools view -b -f 4 results/alignments/nuclear.bam > unmapped.bam
samtools fastq unmapped.bam > unmapped.fastq
```

Align unmapped reads to the mitochondrial genome:

```
bwa mem data/reference/mitochondrial_genome.fa \
  unmapped.fastq | \
  samtools sort -o results/alignments/mt.bam

samtools index results/alignments/mt.bam
```

Identify candidate NUMTs by aligning the mitochondrial genome to the nuclear genome:

```
blastn \
  -query data/reference/mitochondrial_genome.fa \
  -subject data/reference/nuclear_genome.fa \
  -outfmt 6 \
  > results/numt_candidates/mt_vs_nuclear.tsv
```

Convert BLAST output into BED format:

```
awk '{print $2"\t"$9"\t"$10"\t"$1}' \
  results/numt_candidates/mt_vs_nuclear.tsv \
  > results/numt_candidates/numts.bed
```

Filter candidate NUMTs by minimum length and quality:

```
awk '($3 - $2) >= 100' \
  results/numt_candidates/numts.bed \
  > results/numt_candidates/filtered_numts.bed
```

Assess read support for candidate NUMTs:

```
bedtools coverage \
  -a results/numt_candidates/filtered_numts.bed \
  -b results/alignments/nuclear.bam \
  > results/numt_candidates/coverage.txt
```

Summarize NUMT characteristics in R:

```
library(tidyverse)

numts <- read.table("results/numt_candidates/coverage.txt", header = FALSE)
colnames(numts) <- c("chr", "start", "end", "id", "coverage")

summary_stats <- numts %>%
  mutate(length = end - start) %>%
  summarize(
    total_numts = n(),
    mean_length = mean(length),
    mean_coverage = mean(coverage)
  )

print(summary_stats)
```

Optional visualization of NUMT length distribution:

```
ggplot(numts, aes(x = end - start)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(
    title = "NUMT Length Distribution",
    x = "Length (bp)",
    y = "Count"
  )
```

This pipeline provides a general framework for NUMT detection. Parameters such as alignment thresholds, filtering criteria, and minimum coverage should be adjusted based on organism, sequencing depth, and research goals. Additional validation using long-read sequencing or PCR is recommended for confirming candidate NUMTs.
