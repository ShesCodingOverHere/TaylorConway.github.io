---
title: "Mutation Depth Calculation Pipeline"
---

# Mutation Depth Calculation Pipeline

This script calculates sequencing depth across candidate mutation regions for multiple samples. It uses `samtools depth` to extract per-base coverage across specified genomic intervals and computes the average depth per mutation per sample.

The workflow reads a list of candidate mutation regions, iterates through associated samples, and calculates average sequencing depth for each region.

---

## Script

```python
from __future__ import division
import subprocess

# =========================
# Function: calculate depth across a region
# =========================
def calculate_depth(chrom, start, stop, bam_file):

    region = chrom + ":" + str(start) + "-" + str(stop)
    length = stop - start + 1

    command = ["samtools", "depth", "-r", region, bam_file]
    process = subprocess.Popen(command, stdout=subprocess.PIPE)

    total_depth = 0

    for line in process.stdout:
        fields = line.split()
        depth = int(fields[2])
        total_depth += depth

    return region, total_depth, length


# =========================
# Input file (generic path)
# =========================
InFile = open("data/input/mutation_clusters.out", 'r')


# =========================
# Process each mutation
# =========================
for line in InFile:

    A = line.split()

    mutation_id = A[0]
    chrom = A[4]
    start = int(A[5])
    stop = int(A[3])

    # List of samples containing this mutation
    samples = A[7].split(',')

    depths = []

    for sample in samples:

        bam_file = "data/bam/" + sample + ".bam"

        region, total_depth, length = calculate_depth(chrom, start, stop, bam_file)

        depths.append(total_depth)

        # Calculate average depth across region
        average_depth = sum(depths) / length

        print(sample, region, average_depth)
```

---

## Output

The script prints:

```
sample  chrom:start-stop  average_depth
```

---

## Notes

- Uses `samtools depth` to calculate per-base coverage  
- Average depth is computed as total depth divided by region length  
- Input mutation file must include genomic coordinates and associated samples  
- BAM files must be indexed and accessible  
- Designed for multi-sample depth comparison across candidate mutation regions  
