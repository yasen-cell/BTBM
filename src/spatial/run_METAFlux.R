# ============================================================
# METAFlux analysis for spatial transcriptomic data
# (Mouse gene symbols have been converted to human gene symbols)
# ============================================================

# Load required package
library(METAFlux)

# Load reference data
load("data/human_gem.rda")      # human GEM database
load("data/human_blood.rda")    # medium for human samples (used for all spots)

# Load spatial expression data
# Format: rows = genes (human gene symbols), columns = spots
load("spatial_data.rda")   # object `spatial_data` is a gene x spots matrix

# Calculate reaction scores from spatial expression
scores <- calculate_reaction_score(spatial_data)

# Compute metabolic fluxes
# All spots are analyzed using human_blood medium
flux <- compute_flux(mras = scores, medium = human_blood)

# Cube-root transformation to stabilize variance and reduce skewness
cbrt <- function(x) {
    sign(x) * abs(x)^(1/3)
}
flux <- cbrt(flux)

