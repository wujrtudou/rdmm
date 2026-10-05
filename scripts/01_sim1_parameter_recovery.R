# Simulation 1: original tail-weight experiment, expanded to recover parameters
# beyond degrees of freedom, as requested by Reviewer 1.

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

nu_grid <- c(3, 3.5, 4, 4.5, 5, 5.5)
if (cond_index > 0L) nu_grid <- nu_grid[cond_index]
pars <- make_manuscript_parameters()

all_fit <- list()
all_rec <- list()

for (nu in nu_grid) {
  cat(sprintf("\nSimulation 1: nu = %.1f, replications %d--%d\n", nu, rep_start, rep_end))

  worker <- function(replication) {
    data_seed <- 100000L + as.integer(round(100 * nu)) * 10000L + replication
    fit_seed <- 300000L + as.integer(round(100 * nu)) * 10000L + replication
    sim <- simulate_student_rdmm(n = 1000L, nu = nu, pars = pars, seed = data_seed)

    fb <- fit_both_models(
      sim, r = c(5L, 2L), seed = fit_seed,
      max_iter = max_iter, eps = 1e-4,
      initial_nu = 6, nu_structure = "common",
      min_iter = 100L, moving_window = 20L
    )

    rr <- fb$results
    rr$Nu_true <- nu
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
      rec$Nu_true <- nu
      rec$Replication <- replication
      rec$Data_seed <- data_seed
      rec$Fit_seed <- fit_seed
    }
    list(fit = rr, recovery = rec)
  }

  ans <- run_parallel(seq.int(rep_start, rep_end), worker, cores)
  all_fit[[as.character(nu)]] <- do.call(rbind, lapply(ans, function(x) x$fit))
  recs <- Filter(Negate(is.null), lapply(ans, function(x) x$recovery))
  all_rec[[as.character(nu)]] <- if (length(recs)) do.call(rbind, recs) else NULL
}

fit_out <- do.call(rbind, all_fit)
rec_out <- do.call(rbind, Filter(Negate(is.null), all_rec))

raw_dir <- file.path(project_root, "results", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- sprintf("rep%04d_%04d", rep_start, rep_end)
if (cond_index > 0L) suffix <- paste0("cond", cond_index, "_", suffix)
write.csv(fit_out, file.path(raw_dir, paste0("sim1_fit_", suffix, ".csv")), row.names = FALSE)
write.csv(rec_out, file.path(raw_dir, paste0("sim1_recovery_", suffix, ".csv")), row.names = FALSE)
cat("Simulation 1 shard written.\n")
