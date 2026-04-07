---
layout: single
title: "Gene Coordinate Processing"
permalink: /files/genes/
---

# Gene Overlap Classification Pipeline

This script determines how structural mutations overlap with annotated gene regions. It compares mutation coordinates against gene coordinates and classifies each mutation as occurring within a gene, partially overlapping (left or right), or spanning an entire gene.

The workflow reads a gene annotation file and a mutation cluster file, then evaluates positional relationships between mutations and genes.

---

## Script

```python
# =========================
# Input files (generic paths)
# =========================
gene_file = open("data/reference/gene_coordinates.txt", 'r')
mutation_file = open("data/input/mutation_clusters.out", 'r')

# Store gene annotations
genes = []

for line in gene_file:
    genes.append(line)


# =========================
# Compare mutations to genes
# =========================
for line in mutation_file:

    A = line.split()

    if not A:
        continue

    chrom = A[1]
    mut_start = int(A[5])
    mut_stop = int(A[3])

    for gene in genes:

        B = gene.split()

        if not B:
            continue

        gene_start = int(B[4])
        gene_stop = int(B[5])

        # Extract chromosome from annotation
        chrom_info = B[2].split("=")
        gene_chrom = chrom_info[1]

        if gene_chrom == chrom:

            # Mutation completely within gene
            if gene_start < mut_start and gene_stop > mut_stop:
                print(B[0], chrom,
                      str(mut_start) + "." + str(mut_stop),
                      str(gene_start) + "." + str(gene_stop),
                      "WithinMutation")

            # Mutation overlaps right side of gene
            elif mut_start > gene_start and mut_start < gene_stop and mut_stop > gene_stop:
                print(B[0], chrom,
                      str(mut_start) + "." + str(mut_stop),
                      str(gene_start) + "." + str(gene_stop),
                      "PartialRightMutation")

            # Mutation overlaps left side of gene
            elif mut_start < gene_start and mut_stop > gene_start and mut_stop < gene_stop:
                print(B[0], chrom,
                      str(mut_start) + "." + str(mut_stop),
                      str(gene_start) + "." + str(gene_stop),
                      "PartialLeftMutation")

            # Mutation spans entire gene
            elif mut_start < gene_start and gene_stop < mut_stop:
                print(B[0], chrom,
                      str(mut_start) + "." + str(mut_stop),
                      str(gene_start) + "." + str(gene_stop),
                      "WholeMutation")

            else:
                continue
```

---

## Output

Each line contains:

```
gene_id  chrom  mutation_range  gene_range  classification
```

Where classification is one of:
- `WithinMutation`
- `PartialRightMutation`
- `PartialLeftMutation`
- `WholeMutation`

---

## Notes

- Assumes gene annotation file contains chromosome and coordinate information  
- Mutation coordinates are compared directly to gene boundaries  
- Classification is based on positional overlap logic  
- Designed for annotating structural variants relative to gene regions  
- Can be extended to include strand or functional annotation data  
