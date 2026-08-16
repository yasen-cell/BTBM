#!/bin/bash
# COMPASS analysis for bulk transcriptomic data
# Example datasets: TNBC, MTC84, NCI60, ccRCC

# Input expression matrix (rows: genes, columns: samples)
DATA="expression.tsv"

# Number of parallel processes
NUM_PROCESSES=10

# Species (homo_sapiens for human)
SPECIES="homo_sapiens"

# Run COMPASS
compass --data "$DATA" --num-processes "$NUM_PROCESSES" --species "$SPECIES"