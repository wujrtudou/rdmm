# Simulation 2 revised for Reviewer 1: clean observations are generated from
# a true Gaussian DGMM, then contaminated by nested multivariate t4 perturbations.

project_root <- Sys.getenv("PROJECT_ROOT", unset = ".")
source(file.path(project_root, "R", "load_models.R"))
load_rdmm_models(project_root)
source(file.path(project_root, "R", "simulation_parameters.R"))
source(file.path(project_root, "R", "simulation_helpers.R"))

n_rep <- read_env_int("N_REP", 500L)
rep_start <- read_env_int("REP_START", 1L)
rep_end <- read_env_int("REP_END", n_rep)
cores <- read_env_int("CORES", 1L)
max_iter <- read_env_int("MAX_ITER", 300L)
eps_grid <- c(0, .1, .2, .3, .4, .5)
pars <- make_manuscript_parameters()

worker <- function(replication) {
  data_seed <- 910000L + replication
  contamination_seed <- 920000L + replication
  fit_seed <- 930000L + replication

  clean <- simulate_gaussian_dgmm(n = 1000L, pars = pars, seed = data_seed)
  contam <- make_nested_t_contamination(
    clean$y, eps_grid = eps_grid, df = 4, scale_multiplier = 2,
    seed = contamination_seed
  )

  rows <- vector("list", length(eps_grid))
  for (j in seq_along(eps_grid)) {
    eps <- eps_grid[j]
    sim <- clean
    sim$y <- contam$levels[[j]]$y

    # A moderate/high starting df allows RDMM to approach the Gaussian boundary
    # under clean data while retaining the ability to adapt downward under contamination.
    fb <- fit_both_models(
      sim, r = c(5L, 2L), seed = fit_seed,
      max_iter = max_iter, eps = 1e-4,
      initial_nu = 30, nu_structure = "common",
      min_iter = 100L, moving_window = 20L
    )
    rr <- fb$results
    rr$Contamination <- eps
    rr$Replication <- replication
    rr$Data_seed <- data_seed
    rr$Contamination_seed <- contamination_seed
    rr$Fit_seed <- fit_seed
    rows[[j]] <- rr
  }
  do.call(rbind, rows)
}

cat(sprintf("\nGaussian contamination simulation, reps %d--%d\n", rep_start, rep_end))
ans <- run_parallel(seq.int(rep_start, rep_end), worker, cores)
out <- do.call(rbind, ans)

raw_dir <- file.path(project_root, "results", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- sprintf("rep%04d_%04d", rep_start, rep_end)
write.csv(out, file.path(raw_dir, paste0("gaussian_contamination_", suffix, ".csv")), row.names = FALSE)
cat("Gaussian-contamination shard written.\n")
