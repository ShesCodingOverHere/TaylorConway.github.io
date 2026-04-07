---
title: "Site Frequency Spectrum Projection"
---

# Site Frequency Spectrum Projection

This script computes the site frequency spectrum (SFS) from allele frequency data by projecting observed allele counts to a standardized sample size using the hypergeometric distribution. This allows comparison across datasets with different sequencing depths or sample sizes.

The workflow reads allele frequency data, applies hypergeometric projection to a fixed sample size, and outputs a normalized SFS.

---

## Script

```python
from __future__ import division
from scipy import special
import numpy as np
import sys

# =========================
# Helper function: binomial coefficient
# =========================
def choose(n, k):
    return special.binom(n, k)


# =========================
# Parameters
# =========================
max_sample_size = 19

# Initialize SFS vector
SFS = [0] * (max_sample_size + 1)

# Input file (generic path)
InFile = open("data/input/allele_frequencies.out", 'r')


# =========================
# Compute Site Frequency Spectrum
# =========================
for line in InFile:

    # Skip header lines
    if line.startswith("Chr"):
        continue

    A = line.split()

    # Parse numerator / denominator
    freq = A[1].split('/')
    numerator = int(freq[0])
    denominator = int(freq[1])

    # Skip empty sites
    if numerator == 0 and denominator == 0:
        continue

    # Only project if sufficient sample size
    if denominator >= max_sample_size:

        for i in range(0, max_sample_size + 1):

            SFS[i] += (
                choose(numerator, i) *
                choose(denominator - numerator, max_sample_size - i) /
                choose(denominator, max_sample_size)
            )

            # Check for numerical issues
            if np.isnan(SFS[i]):
                print(line)
                sys.exit()


# Remove zero-frequency class
SFS = SFS[1:]

# Normalize SFS
total = sum(SFS)

for value in SFS:
    print(value / total)
```

---

## Output

The script prints a normalized site frequency spectrum:

```
f1
f2
f3
...
fn
```

Where each value represents the proportion of sites at a given allele frequency.

---

## Notes

- Uses hypergeometric projection to standardize allele counts  
- Requires allele frequency input in `numerator/denominator` format  
- Skips sites with zero coverage  
- Output is normalized so values sum to 1  
- Designed for population genetic inference and comparison across datasets  
