# Benchmarking transcriptome-based metabolic inference methods across bulk, single-cell, and spatial transcriptomic data
This GitHub repository contains the core inference scripts and downstream analysis code for our manuscript entitled "Benchmarking transcriptome-based metabolic inference methods across bulk, single-cell, and spatial transcriptomic data".

## Overview

This project systematically benchmarks **eight transcriptome-based metabolic inference methods** across **14 datasets** spanning bulk, single-cell, and spatial transcriptomics. The benchmark evaluates:

- Correlation with ground-truth metabolomic measurements (e.g., SUV, exchange fluxes, LC–MS/GC–MS, MALDI–MSI)
- Ability to resolve cell-type and subtype heterogeneity
- Recovery of pathway-level differential metabolism
- Spatial pattern detection (dopamine biosynthesis, histidine metabolism)
- Computational efficiency

## Repository Structure

BTBM/
├── README.md                         # This file
├── src/                              # Core scripts for running each method
│   ├── bulk/                         # Bulk RNA-seq inference scripts
│   │   ├── run_scMetabolism.R
│   │   ├── run_scFEA.ipynb
│   │   ├── run_METAFlux.R
│   │   ├── run_Compass.sh
│   │   └── run_scFBA.m
│   ├── single_cell/                  # Single-cell inference scripts
│   │   ├── run_scMetabolism.R
│   │   ├── run_scFEA.ipynb
│   │   ├── run_METAFlux.R
│   │   └── run_Compass.sh
│   └── spatial/                      # Spatial transcriptomics inference scripts
│       ├── run_scMetabolism.R
│       ├── run_scFEA.ipynb
│       ├── run_METAFlux.R
│       └── run_Compass.sh
└── analysis/                         # Downstream analysis scripts (to be added)
    └── (empty in this initial release)


