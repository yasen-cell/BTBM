# ============================================================
# METAFlux analysis for single-cell transcriptomic data
# (Mouse gene symbols have been converted to human gene symbols)
# ============================================================

# Load required package
library(METAFlux)

# Load reference data
load("data/human_gem.rda")      # human GEM database
load("data/human_blood.rda")    # medium for human samples (used for all cells)

# Load single-cell expression data
# Format: rows = genes (human gene symbols), columns = cells
load("sc_data.rda")   # object `sc_expr` is a gene x cell matrix

# Calculate reaction scores from single-cell expression
scores <- calculate_reaction_score(sc_data)

# Compute metabolic fluxes
# All samples/cells are analyzed using human_blood medium
flux <- compute_flux(mras = scores, medium = human_blood)

# Cube-root transformation to stabilize variance and reduce skewness
cbrt <- function(x) {
    sign(x) * abs(x)^(1/3)
}
flux <- cbrt(flux)

