#!/bin/bash
# COMPASS analysis for single cell transcriptomic data
# Example datasets: Mouse Liver, Melanoma, HNSCC

# Input expression matrix (rows: genes, columns: samples)
DATA="expression.tsv"

# Number of parallel processes
NUM_PROCESSES=10

# Species (homo_sapiens for human/mus_musculus for mouse)
SPECIES="homo_sapiens"
#SPECIES="mus_musculus"
# Run COMPASS
compass --data "$DATA" --num-processes "$NUM_PROCESSES" --species "$SPECIES"