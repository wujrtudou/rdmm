# Sensitivity to the diagonal residual-covariance assumption.
#
# Data are generated from the same two-layer RDMM used in Simulation 1, except
# that Psi_l is allowed to have AR(1) off-diagonal correlation.  Both RDMM and
# DGMM are then fitted with the manuscript's current diagonal-Psi specification.
# This directly measures robustness to misspecification of the diagonal
# residual-covariance assumption without changing the fitted model.

project_root <- Sys.getenv("PROJECT_ROOT", unset = ".")
source(file.path(project_root, "R", "load_models.R"))
load_rdmm_models(project_root)
source(file.path(project_root, "R", "simulation_parameters.R"))
source(file.path(project_root, "R", "simulation_helpers.R"))

n_rep <- read_env_int("N_REP", 500L)
rep_start <- read_env_int("REP_START", 1L)
rep_end <- read_env_int("REP_END", n_rep)
cores <- read_env_int("CORES", 1L)
cond_index <- read_env_int("COND", 0L)
max_iter <- read_env_int("MAX_ITER", 300L)

# rho = 0 reproduces the manuscript assumption; larger values progressively
# violate it while preserving all marginal residual variances.
rho_grid <- c(0.0, 0.2, 0.4, 0.6)
if (cond_index > 0L) rho_grid <- rho_grid[cond_index]

nu_true <- 4
base_pars <- make_manuscript_parameters()
all_fit <- list()
all_rec <- list()

for (rho in rho_grid) {
  cat(sprintf("\nResidual-covariance sensitivity: rho = %.1f, reps %d--%d\n",
              rho, rep_start, rep_end))
  pars <- make_correlated_residual_parameters(rho1 = rho, rho2 = rho,
                                              pars = base_pars)

  worker <- function(replication) {
    # Use the same underlying RNG stream across rho values within a replication.
    # This makes the comparison across correlation strengths more nearly paired.
    data_seed <- 820000L + replication
    fit_seed <- 920000L + replication
    sim <- simulate_student_rdmm(n = 1000L, nu = nu_true,
                                 pars = pars, seed = data_seed)

    fb <- fit_both_models(
      sim, r = c(5L, 2L), seed = fit_seed,
      max_iter = max_iter, eps = 1e-4,
      initial_nu = 6, nu_structure = "common",
      min_iter = 100L, moving_window = 20L
    )

    rr <- fb$results
    rr$Residual_rho <- rho
    rr$Nu_true <- nu_true
    rr$Replication <- replication
    rr$Data_seed <- data_seed
    rr$Fit_seed <- fit_seed

    rec <- NULL
    if (!is.null(fb$fits$RDMM)) {
      rec <- tryCatch(parameter_recovery_metrics(fb$fits$RDMM, sim), error = identity)
      if (inherits(rec, "error")) {
        rec <- data.frame(
          Nu_hat = NA_real_, Nu_error = NA_real_, Pi1_RMSE = NA_real_, Pi2_RMSE = NA_real_,
          Eta1_RMSE = NA_real_, Psi1_RMSE = NA_real_, Loading1_Subspace_Error = NA_real_,
          PathMean_RMSE = NA_real_, PathScale_RelFrob = NA_real_,
          Recovery_Error = conditionMessage(rec), stringsAsFactors = FALSE
        )
      } else {
        rec$Recovery_Error <- ""
      }
      rec$Residual_rho <- rho
      rec$Nu_true <- nu_true
      rec$Replication <- replication
      rec$Data_seed <- data_seed
      rec$Fit_seed <- fit_seed
    }

    list(fit = rr, recovery = rec)
  }

  ans <- run_parallel(seq.int(rep_start, rep_end), worker, cores)
  all_fit[[as.character(rho)]] <- do.call(rbind, lapply(ans, function(x) x$fit))
  recs <- Filter(Negate(is.null), lapply(ans, function(x) x$recovery))
  all_rec[[as.character(rho)]] <- if (length(recs)) do.call(rbind, recs) else NULL
}

fit_out <- do.call(rbind, all_fit)
rec_out <- do.call(rbind, Filter(Negate(is.null), all_rec))

raw_dir <- file.path(project_root, "results", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- sprintf("rep%04d_%04d", rep_start, rep_end)
if (cond_index > 0L) suffix <- paste0("cond", cond_index, "_", suffix)
write.csv(fit_out, file.path(raw_dir, paste0("covariance_sensitivity_fit_", suffix, ".csv")),
          row.names = FALSE)
write.csv(rec_out, file.path(raw_dir, paste0("covariance_sensitivity_recovery_", suffix, ".csv")),
          row.names = FALSE)
cat("Residual-covariance sensitivity shard written.\n")
