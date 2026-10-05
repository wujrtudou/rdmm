# Convenience runner. For publication-scale Monte Carlo experiments, running
# the main studies as separate/sharded jobs is recommended.  The computational
# cost benchmark is intentionally sequential and uses COST_N_REP (default 20).
project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), mustWork = TRUE)
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "scripts", "01_sim1_parameter_recovery.R"))
source(file.path(project_root, "scripts", "02_dimension_misspecification.R"))
source(file.path(project_root, "scripts", "03_gaussian_contamination.R"))
source(file.path(project_root, "scripts", "04_simulation_visuals_and_convergence.R"))
source(file.path(project_root, "scripts", "07_covariance_sensitivity.R"))
source(file.path(project_root, "scripts", "08_computational_cost.R"))
source(file.path(project_root, "scripts", "05_make_tables_and_figures.R"))
