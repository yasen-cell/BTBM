#!/bin/bash
# COMPASS analysis for spatial transcriptomic data
# Example datasets: Human Striatum,Mouse SubNigra

# Input expression matrix (rows: genes, columns: samples)
DATA="expression.tsv"

# Number of parallel processes
NUM_PROCESSES=10

# Species (homo_sapiens for human/mus_musculus for mouse)
SPECIES="homo_sapiens"
#SPECIES="mus_musculus"
# Run COMPASS
compass --data "$DATA" --num-processes "$NUM_PROCESSES" --species "$SPECIES"