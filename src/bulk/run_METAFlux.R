# Load required package
library(METAFlux)

# Load reference data
load("data/human_gem.rda")      # human GEM database
load("data/cell_medium.rda")    # culture medium for cell line samples (e.g., NCI60)
load("data/human_blood.rda")    # medium for human tissue/primary samples (e.g., TNBC, MTC84, ccRCC)

# Load bulk transcriptomic data (e.g., TNBC, MTC84, NCI60, or ccRCC)
load("bulk_data.rda")

# Calculate reaction scores from bulk expression
scores <- calculate_reaction_score(bulk_data)

# Compute metabolic fluxes using the appropriate medium.
# For cell line data (e.g., NCI60), use `cell_medium`.
# For human tissue samples (e.g., TNBC, MTC84, ccRCC), use `human_blood`.
flux <- compute_flux(mras = scores, medium = human_blood)

# Cube-root transformation to stabilize variance and reduce skewness
cbrt <- function(x) {
    sign(x) * abs(x)^(1/3)
}
flux <- cbrt(flux)