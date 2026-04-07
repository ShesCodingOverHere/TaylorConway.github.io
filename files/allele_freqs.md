---
layout: single
title: "Allele Frequency Calculation"
permalink: /files/allele_freqs/
---

# Allele Frequency Calculation

This script calculates allele frequencies from haplotype output files. It can be used for both combined population datasets and subpopulation-specific datasets following haplotype calling.

The script parses an input file containing haplotype counts, aggregates counts by genomic position, and computes allele frequencies as the sum of numerators over denominators.

---

## Usage Notes

- Use a combined haplotype file (e.g., `combined_haplotypes.out`) to calculate overall allele frequencies.
- Use subpopulation-specific files (e.g., `subpopulation_*.out`) to calculate allele frequencies within subsets.
- When calculating subpopulation allele frequencies, remove or comment out the filtering step that excludes zero-numerator sites.

---

## Script

```python
# Open input haplotype file
InFile = open("data/input/combined_haplotypes.out", 'r')

# Dictionary to store:
# key = chrom.start.stop
# value = list of (numerator, denominator, sample_id)
MyDict = {}

# Parse input file and populate dictionary
for line in InFile:
    A = line.split()

    # Construct key from genomic coordinates
    mykey = A[1] + '.' + A[2] + '.' + A[3]

    # Append values to dictionary
    if mykey in MyDict:
        MyDict[mykey] = MyDict[mykey] + [(A[-2], A[-1], A[0])]
    else:
        MyDict[mykey] = [(A[-2], A[-1], A[0])]

# Iterate through each genomic position
for key in MyDict:

    # Lists to store numerators and denominators
    num = []
    denom = []

    # Process each entry for the position
    for i in MyDict[key]:

        # Coverage correction:
        # If numerator = 2 and denominator = 1 (invalid case),
        # correct numerator to 1
        if int(i[0]) == 2 and int(i[1]) == 1:
            num.append(1)
            denom.append(1)
        else:
            num.append(int(i[0]))
            denom.append(int(i[1]))

    # Sum counts across samples
    numerator = sum(num)
    denominator = sum(denom)

    # IMPORTANT:
    # When running on combined populations, exclude sites where numerator = 0
    # (e.g., to avoid reference contamination effects).
    # Comment out this block when analyzing subpopulations.
    if numerator == 0:
        continue
    else:
        print(key, str(numerator) + '/' + str(denominator))
```

---

## Output

The script prints allele frequencies in the format:

```
chrom.start.stop    numerator/denominator
```

---

## Notes

- This script assumes consistent formatting of haplotype output files.
- The coverage correction step prevents biologically impossible counts from propagating.
- Filtering behavior should be adjusted depending on whether data represent combined populations or subpopulations.
