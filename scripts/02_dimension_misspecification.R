# Reviewer-requested sensitivity to slight misspecification of r1 and r2.
# Data are generated at the true (r1,r2)=(5,2), nu=4.

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

# One-at-a-time slight under/over specification plus the true architecture.
designs <- data.frame(
  Architecture = c("r1 under", "r2 under", "true", "r1 over", "r2 over"),
  r1 = c(4L, 5L, 5L, 6L, 5L),
  r2 = c(2L, 1L, 2L, 2L, 3L),
  stringsAsFactors = FALSE
)
if (cond_index > 0L) designs <- designs[cond_index, , drop = FALSE]

pars <- make_manuscript_parameters()
nu_true <- 4
out <- list()

for (ii in seq_len(nrow(designs))) {
  dd <- designs[ii, ]
  cat(sprintf("\nDimension sensitivity: %s (%d,%d), reps %d--%d\n",
              dd$Architecture, dd$r1, dd$r2, rep_start, rep_end))

  worker <- function(replication) {
    # Same generated data for every fitted architecture within a replication.
    data_seed <- 510000L + replication
    fit_seed <- 710000L + replication
    sim <- simulate_student_rdmm(n = 1000L, nu = nu_true, pars = pars, seed = data_seed)
    fb <- fit_both_models(
      sim, r = c(dd$r1, dd$r2), seed = fit_seed,
      max_iter = max_iter, eps = 1e-4,
      initial_nu = 6, nu_structure = "common",
      min_iter = 100L, moving_window = 20L
    )
    rr <- fb$results
    rr$Architecture <- dd$Architecture
    rr$r1_fit <- dd$r1
    rr$r2_fit <- dd$r2
    rr$r1_true <- 5L
    rr$r2_true <- 2L
    rr$Nu_true <- nu_true
    rr$Replication <- replication
    rr$Data_seed <- data_seed
    rr$Fit_seed <- fit_seed
    rr
  }

  ans <- run_parallel(seq.int(rep_start, rep_end), worker, cores)
  out[[ii]] <- do.call(rbind, ans)
}

out <- do.call(rbind, out)
raw_dir <- file.path(project_root, "results", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- sprintf("rep%04d_%04d", rep_start, rep_end)
if (cond_index > 0L) suffix <- paste0("cond", cond_index, "_", suffix)
write.csv(out, file.path(raw_dir, paste0("dimension_misspec_", suffix, ".csv")), row.names = FALSE)
cat("Dimension-misspecification shard written.\n")
