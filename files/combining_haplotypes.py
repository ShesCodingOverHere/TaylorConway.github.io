---
title: "Haplotype Integration and Annotation Pipeline"
---

# Haplotype Integration and Annotation Pipeline

This script integrates mutation data with haplotype blocks across multiple samples to determine whether mutations fall within heterozygous or inbred regions. It processes two mutation input files and compares them against haplotype annotations, assigning each mutation a classification based on its genomic context.

The workflow reads haplotype files, organizes them by chromosome for efficient lookup, and then annotates each mutation by determining whether it falls within a heterozygous haplotype block.

---

## Script

```python
import os

# =========================
# Input files (generic paths)
# =========================
InFile1 = open("data/input/compared_depths.out", 'r')
InFile2 = open("data/input/zero_coverage.out", 'r')

# Directory containing haplotype files
haplo_dir = "data/haplotypes/"

# List all haplotype files
haplo_files = os.listdir(haplo_dir)

# =========================
# Organize haplotypes by chromosome
# =========================
chromDict = {'2L':[], '2R':[], '3L':[], '3R':[], '4':[], 'X':[]}

for hap_file in haplo_files:
    InFile = open(haplo_dir + hap_file)

    for line in InFile:
        splitLine = line.split()
        chrom = splitLine[0]

        # Store entire line for later parsing
        chromDict[chrom].append(line.rstrip())


# =========================
# Optional: strains to exclude
# =========================
exclude_strains = [
    "sampleA", "sampleB", "sampleC"
]


# =========================
# Function to annotate mutations
# =========================
def annotate_mutations(InFile):

    for line in InFile:
        A = line.split()

        strain = A[0]

        # Skip excluded strains
        if strain in exclude_strains:
            continue

        # Parse mutation coordinates
        if ":" in A[1]:
            chrom_part = A[1].split(':')
            chrom = chrom_part[0]
            coords = chrom_part[1].split('-')
            start = int(coords[0])
            stop = int(coords[1])
            hap1 = A[5]
        else:
            chrom = A[1]
            start = int(A[2])
            stop = int(A[3])
            hap1 = A[4]

        # Default classification = inbred
        inbred = 1

        # Search matching haplotype blocks for this chromosome
        for entry in chromDict.get(chrom, []):
            X = entry.split()

            hapStrain = X[-1]

            if hapStrain == strain:
                start2 = int(X[1])
                stop2 = int(X[2])

                # Check if mutation lies within haplotype block
                if start >= start2 and stop <= stop2:

                    # If haplotype is heterozygous, update classification
                    if X[3] == "hetero":
                        inbred = 2

                    break

        # Output annotated mutation
        print(strain, chrom, start, stop, hap1, inbred)


# =========================
# Run annotation on both input files
# =========================
annotate_mutations(InFile1)
annotate_mutations(InFile2)
```

---

## Output

Each line of output contains:

```
strain  chrom  start  stop  haplotype  classification
```

Where:
- `classification = 1` → inbred region  
- `classification = 2` → heterozygous region  

---

## Notes

- Haplotype files are pre-organized by chromosome for efficient lookup  
- Mutations are matched to haplotype blocks based on coordinate overlap  
- Strain filtering can be customized depending on dataset  
- Assumes haplotype files contain chromosome, start, stop, genotype, and strain ID  
- Designed for batch processing of multiple mutation datasets  
