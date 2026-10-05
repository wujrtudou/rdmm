# Run only the two new reviewer-driven analyses, then rebuild tables/figures.
#
# Covariance sensitivity uses the standard Monte Carlo variables N_REP,
# REP_START, REP_END, COND and CORES (default N_REP=500).
# Computational cost uses separate COST_* variables and is always sequential
# (default COST_N_REP=20).

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), mustWork = TRUE)
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "scripts", "07_covariance_sensitivity.R"))
source(file.path(project_root, "scripts", "08_computational_cost.R"))
source(file.path(project_root, "scripts", "05_make_tables_and_figures.R"))
source(file.path(project_root, "scripts", "06_validate_outputs.R"))
