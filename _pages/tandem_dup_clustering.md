---
layout: single
title: "Tandem Duplication Clustering Pipeline"
permalink: /code/tandem_dup_clustering/
author_profile: true
---

# Tandem Duplication Clustering Pipeline

This script clusters putative mutations across multiple samples based on genomic proximity. It parses per-sample mutation files, bins mutations by genomic position, and then merges nearby events into clusters shared across samples.

The workflow aggregates mutation calls, identifies clusters within a defined distance threshold, and outputs consolidated regions with the number of supporting samples.

NOTE: This script was originally written for Python 2.

---

## Script

```python
from __future__ import division

# Dictionaries for storing data at different stages
Dict1 = {}  # raw binned mutations
Dict2 = {}  # clustered mutations
Dict3 = {}  # final clustered regions

# Tracks positions already considered during clustering
PosList = []

# Bin size used for grouping nearby mutations
binsize = 325


# =========================
# Step 1: Load mutation files into bins
# =========================
def firstdict(filename):

    # Generic input path
    InFile = open("data/input/" + filename, 'r')

    for line in InFile:
        A = line.split()

        # Parse chromosome and coordinates
        parts = A[0].split(":")
        chrom = parts[0]

        # Expected format: chrom:start.stop
        span = parts[1].split(".")
        minimum = int(span[0])
        maximum = int(span[1])

        # Supporting read count
        count = int(A[1])

        # Sample identifier derived from filename
        sample = filename.split("_")[0]

        # Assign bins
        mybin1 = minimum - minimum % binsize
        mybin2 = maximum - maximum % binsize

        mykey = chrom + ":" + str(mybin1) + "." + str(mybin2)

        # Keep only mutations supported by ≥3 reads
        if count >= 3:
            if mykey in Dict1:
                Dict1[mykey].append([sample, minimum, maximum])
            else:
                Dict1[mykey] = [[sample, minimum, maximum]]


# =========================
# Step 2: Load all samples
# =========================
files = [
    "sample01_divreads.out",
    "sample02_divreads.out",
    "sample03_divreads.out",
    "sample04_divreads.out",
    "sample05_divreads.out"
]

for f in files:
    firstdict(f)


# =========================
# Step 3: Cluster mutations across neighboring bins
# =========================
def binchecking(key1, key2):

    if key1 in Dict1 and key2 in Dict1:

        run = 0
        for i in Dict1[key1]:

            count = 0
            sample1 = Dict1[key1][run][0]
            pos1a = Dict1[key1][run][1]
            pos1b = Dict1[key1][run][2]

            for j in Dict1[key2]:

                sample2 = Dict1[key2][count][0]
                pos2a = Dict1[key2][count][1]
                pos2b = Dict1[key2][count][2]

                if pos2a not in PosList and pos2b not in PosList:

                    # Cluster if positions are within bin distance
                    if abs(pos1a - pos2a) < binsize and abs(pos1b - pos2b) < binsize:

                        PosList.extend([pos1a, pos1b, pos2a, pos2b])

                        if key1 in Dict2:
                            Dict2[key1].append([pos2a, pos2b, sample2])
                        else:
                            Dict2[key1] = [
                                [pos1a, pos1b, sample1],
                                [pos2a, pos2b, sample2]
                            ]

                count += 1

            run += 1


# =========================
# Step 4: Compare neighboring bins
# =========================
for key in Dict1.keys():

    chrom, span = key.split(':')
    bin1, bin2 = map(int, span.split('.'))

    neighbor_keys = [
        chrom + ':' + str(bin1 - binsize) + '.' + str(bin2),
        chrom + ':' + str(bin1 + binsize) + '.' + str(bin2),
        chrom + ':' + str(bin1) + '.' + str(bin2 + binsize),
        chrom + ':' + str(bin1) + '.' + str(bin2 - binsize),
        chrom + ':' + str(bin1 - binsize) + '.' + str(bin2 + binsize),
        chrom + ':' + str(bin1 + binsize) + '.' + str(bin2 + binsize),
        chrom + ':' + str(bin1 - binsize) + '.' + str(bin2 - binsize),
        chrom + ':' + str(bin1 + binsize) + '.' + str(bin2 - binsize)
    ]

    for neighbor in neighbor_keys:
        binchecking(key, neighbor)


# =========================
# Step 5: Collapse clusters into final regions
# =========================
for key in Dict2.keys():

    chrom = key.split(":")[0]

    front = []
    back = []
    count = 0

    for entry in Dict2[key]:
        front.append(entry[0])
        back.append(entry[1])
        count += 1

    minimum = min(front)
    maximum = max(back)

    # Keep clusters smaller than threshold
    if maximum - minimum < 25000:
        newkey = chrom + ":" + str(minimum) + "." + str(maximum)
        Dict3[newkey] = count


# =========================
# Step 6: Output results
# =========================
for outkey in Dict3.keys():
    print(outkey, Dict3[outkey])
```

---

## Output

The script outputs clustered mutation regions in the format:

```
chrom:start.stop    count
```

Where:
- `chrom:start.stop` = merged mutation interval  
- `count` = number of supporting samples  

---

## Notes

- Binning is used as a heuristic to group nearby mutations before clustering  
- Clustering thresholds are currently simplistic and may be improved  
- Designed for multi-sample mutation comparison  
- Written for Python 2; minor modifications are needed for Python 3 compatibility  
