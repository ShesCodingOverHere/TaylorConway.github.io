---
layout: single
title: "Population Counts"
permalink: /code/population_counts/
---

# Population Mutation Count Analysis

This script quantifies mutation counts across multiple populations and mutation types. It processes both non-chimeric and chimeric mutation datasets, categorizes mutations by type, and summarizes counts at two levels: per mutation and per strain.

The workflow assigns each strain to a population group and tracks counts for different mutation classes, including within-gene, whole-gene, partial, and chimeric mutations.

---

## Script

```python
# =========================
# Input files (generic paths)
# =========================
NonChimeric = open("data/input/non_chimeric_mutations.out", 'r')
Chimeric = open("data/input/chimeric_mutations.out", 'r')


# =========================
# Population group definitions
# =========================
mainland_samples = {"sampleA1", "sampleA2", "sampleA3"}
island_samples = {"sampleB1", "sampleB2", "sampleB3"}
outgroup_samples = {"sampleC1", "sampleC2", "sampleC3"}


# =========================
# Initialize counters
# =========================
def initialize_counts():
    return {
        "total": 0,
        "within": 0,
        "left": 0,
        "right": 0,
        "whole": 0,
        "chimeric": 0
    }

mainland = initialize_counts()
island = initialize_counts()
outgroup = initialize_counts()


# =========================
# Function to update counts
# =========================
def update_counts(group, muttype):

    group["total"] += 1

    if muttype == "WithinMutation":
        group["within"] += 1
    elif muttype == "WholeMutation":
        group["whole"] += 1
    elif muttype == "PartialRightMutation":
        group["right"] += 1
    elif muttype == "PartialLeftMutation":
        group["left"] += 1
    elif muttype == "chimeric":
        group["chimeric"] += 1


# =========================
# Process non-chimeric mutations
# =========================
for line in NonChimeric:

    A = line.strip().split()
    muttype = A[5]
    strains = A[6].split(',')

    for sample in strains:

        if sample in mainland_samples:
            update_counts(mainland, muttype)

        elif sample in island_samples:
            update_counts(island, muttype)

        elif sample in outgroup_samples:
            update_counts(outgroup, muttype)

        break  # count once per mutation


# =========================
# Process chimeric mutations
# =========================
for line in Chimeric:

    A = line.strip().split()
    muttype = "chimeric"
    strains = A[12].split(',')

    for sample in strains:

        if sample in mainland_samples:
            update_counts(mainland, muttype)

        elif sample in island_samples:
            update_counts(island, muttype)

        elif sample in outgroup_samples:
            update_counts(outgroup, muttype)

        break  # count once per mutation


# =========================
# Output results
# =========================
def print_results(name, group):

    print(name + " Total Mutations: " + str(group["total"]))
    print(name + " Whole Mutations: " + str(group["whole"]))
    print(name + " Within Mutations: " + str(group["within"]))
    print(name + " Partial Right Mutations: " + str(group["right"]))
    print(name + " Partial Left Mutations: " + str(group["left"]))
    print(name + " Chimeric Mutations: " + str(group["chimeric"]))
    print("")


print("Mutation Counts per Population\n")

print_results("Mainland", mainland)
print_results("Island", island)
print_results("Outgroup", outgroup)
```

---

## Output

The script prints summary statistics for each population:

```
Population Total Mutations
Population Whole Mutations
Population Within Mutations
Population Partial Right Mutations
Population Partial Left Mutations
Population Chimeric Mutations
```

---

## Notes

- Populations are defined by sets of sample identifiers  
- Each mutation is counted once per population  
- Chimeric and non-chimeric mutations are processed separately  
- Mutation categories are based on positional relationships within genes  
- Designed for comparative population-level mutation analysis  
