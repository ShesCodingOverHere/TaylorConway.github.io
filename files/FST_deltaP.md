---
title: "Allele Frequency Comparison and FST Calculation"
---

# Allele Frequency Comparison and FST Calculation

This script compares allele frequencies between two populations and calculates both ΔP (allele frequency difference) and FST for each genomic region. It takes precomputed allele frequency files as input, parses genomic coordinates and counts, and outputs summary statistics for population differentiation.

The workflow reads two input files representing separate populations, computes allele frequencies, and then calculates ΔP and FST for each shared genomic interval.

---

## Script

```python
from __future__ import division
import math

# =========================
# Input / Output (generic paths)
# =========================
InFile1 = open("data/input/population1_allele_freqs.out", 'r')
InFile2 = open("data/input/population2_allele_freqs.out", 'r')
OutFile = open("results/output/deltaP_fst_results.tsv", 'w')

p1Dict = {}
p2Dict = {}

# =========================
# Step 1: Read allele frequency files
# =========================
def readDict(Infile, mydict):

    for line in Infile:
        A = line.split()

        # Parse genomic coordinates
        span = A[0].split('.')
        chromosome = span[0]
        start = int(span[1])
        stop = int(span[2])

        # Parse counts
        B = A[1].split('/')
        num = int(B[0])
        denom = int(B[1])

        key = A[0]

        # Calculate allele frequency
        if denom == 0:
            p = 0
        else:
            p = float(num) / float(denom)

        mydict[key] = [chromosome, start, stop, num, denom, p]


readDict(InFile1, p1Dict)
readDict(InFile2, p2Dict)


# =========================
# Step 2: Compute ΔP and FST
# =========================
for key in p1Dict:

    chrom = p1Dict[key][0]
    start = p1Dict[key][1]
    stop = p1Dict[key][2]

    num1 = p1Dict[key][3]
    denom1 = p1Dict[key][4]
    freq1 = p1Dict[key][5]

    num2 = p2Dict[key][3]
    denom2 = p2Dict[key][4]
    freq2 = p2Dict[key][5]

    # Midpoint of region
    midpoint = math.trunc((start + stop) / 2)

    if (denom1 + denom2) != 0:

        # Combined allele frequency
        p_bar = (num1 + num2) / (denom1 + denom2)
        q_bar = 1 - p_bar

        # Population weights
        c1 = denom1 / (denom1 + denom2)
        c2 = denom2 / (denom1 + denom2)

        p1 = freq1
        p2 = freq2
        q1 = 1 - p1
        q2 = 1 - p2

        total_num = num1 + num2
        total_denom = denom1 + denom2
        total_freq = total_num / total_denom

        # FST calculation
        if (p_bar * q_bar) == 0:
            fst = 0
        else:
            fst = ((p_bar * q_bar) - ((p1 * q1 * c1) + (p2 * q2 * c2))) / (p_bar * q_bar)

        # Delta P
        deltaP = freq1 - freq2

        # Output
        OutFile.write(
            chrom + "\t" +
            str(start) + "\t" +
            str(stop) + "\t" +
            str(midpoint) + "\t" +
            str(num1) + "\t" +
            str(denom1) + "\t" +
            str(freq1) + "\t" +
            str(num2) + "\t" +
            str(denom2) + "\t" +
            str(freq2) + "\t" +
            str(deltaP) + "\t" +
            str(fst) + "\t" +
            str(total_freq) + "\n"
        )
```

---

## Output

The output file contains tab-separated values:

```
chrom  start  stop  midpoint  num1  denom1  freq1  num2  denom2  freq2  deltaP  fst  total_freq
```

---

## Notes

- ΔP represents the difference in allele frequency between two populations  
- FST quantifies population differentiation at each genomic region  
- Regions with zero total coverage are excluded from calculations  
- Input files must share the same genomic coordinate keys  
- Designed for downstream population genetics analyses  
