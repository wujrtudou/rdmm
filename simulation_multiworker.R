setwd("C:/Users/uqjwu15/Desktop/rdmm_simulation")

Sys.setenv(
  PROJECT_ROOT = "C:/Users/uqjwu15/Desktop/rdmm_simulation",
  N_REP = 500,
  CORES = 12
)

source("scripts/01_sim1_parameter_recovery.R")
source("scripts/02_dimension_misspecification.R")
source("scripts/03_gaussian_contamination.R")
source("scripts/04_simulation_visuals_and_convergence.R")
source("scripts/07_covariance_sensitivity.R")

# Runtime benchmarking should be sequential for interpretable wall-clock times.
Sys.setenv(COST_N_REP = 20)
source("scripts/08_computational_cost.R")

source("scripts/06_validate_outputs.R")
source("scripts/05_make_tables_and_figures.R")