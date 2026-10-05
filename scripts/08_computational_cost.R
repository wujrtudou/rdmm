# Empirical computational-cost / scalability study for Reviewer 1, Comment 5.
#
# The fitted architecture is held fixed at K=(2,2), r=(5,2).  We vary n at
# p=20 and p at n=1000, timing DGMM and RDMM separately.  Timing benchmarks are
# intentionally run sequentially: parallel replication would make per-fit wall
# clock times depend on CPU contention and would be difficult to interpret.

project_root <- Sys.getenv("PROJECT_ROOT", unset = ".")
source(file.path(project_root, "R", "load_models.R"))
load_rdmm_models(project_root)
source(file.path(project_root, "R", "simulation_parameters.R"))
source(file.path(project_root, "R", "simulation_helpers.R"))

cost_n_rep <- read_env_int("COST_N_REP", 20L)
rep_start <- read_env_int("COST_REP_START", 1L)
rep_end <- read_env_int("COST_REP_END", cost_n_rep)
max_iter <- read_env_int("COST_MAX_ITER", 300L)
min_iter <- read_env_int("COST_MIN_ITER", 100L)
moving_window <- read_env_int("COST_MOVING_WINDOW", 20L)
cond_index <- read_env_int("COST_COND", 0L)

n_grid <- read_env_int_vector("COST_N_GRID", c(250L, 500L, 1000L, 2000L))
p_grid <- read_env_int_vector("COST_P_GRID", c(20L, 40L, 80L))
base_n <- read_env_int("COST_BASE_N", 1000L)
base_p <- read_env_int("COST_BASE_P", 20L)

if (any(n_grid < 20L)) stop("COST_N_GRID contains an impractically small sample size.")
if (any(p_grid < 20L)) stop("COST_P_GRID must use p >= 20 for this benchmark design.")
if (!(base_n %in% n_grid)) n_grid <- sort(unique(c(n_grid, base_n)))
if (!(base_p %in% p_grid)) p_grid <- sort(unique(c(p_grid, base_p)))

# Unique conditions: vary n at p=base_p; vary p at n=base_n.  The common
# baseline (base_n, base_p) is run only once and contributes to both summaries.
conditions <- unique(rbind(
  data.frame(n = n_grid, p = base_p),
  data.frame(n = base_n, p = p_grid)
))
conditions <- conditions[order(conditions$p, conditions$n), , drop = FALSE]
conditions$Condition <- seq_len(nrow(conditions))
if (cond_index > 0L) {
  if (cond_index > nrow(conditions)) stop("COST_COND exceeds the number of benchmark conditions.")
  conditions <- conditions[conditions$Condition == cond_index, , drop = FALSE]
}

nu_true <- 4
base_pars <- make_manuscript_parameters()
raw_dir <- file.path(project_root, "results", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# Save environment information once so reported runtimes can be interpreted.
sys_lines <- c(
  paste("Benchmark timestamp:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R version:", R.version.string),
  paste("Platform:", R.version$platform),
  paste("OS:", paste(Sys.info()[c("sysname", "release", "version", "machine")], collapse = " | ")),
  paste("Detected logical cores:", parallel::detectCores(logical = TRUE)),
  paste("Detected physical cores:", parallel::detectCores(logical = FALSE)),
  paste("Benchmark execution: sequential (one fit at a time)"),
  paste("DGMM/RDMM architecture: K=(2,2), r=(5,2)"),
  paste("nu for generated RDMM data:", nu_true),
  paste("max_iter:", max_iter, "min_iter:", min_iter, "moving_window:", moving_window)
)
writeLines(sys_lines, file.path(raw_dir, "computational_cost_system_info.txt"))

all_rows <- list(); z <- 1L
for (cc in seq_len(nrow(conditions))) {
  n_now <- as.integer(conditions$n[cc])
  p_now <- as.integer(conditions$p[cc])
  condition_id <- as.integer(conditions$Condition[cc])
  pars <- make_runtime_parameters(p = p_now, pars = base_pars)

  cat(sprintf("\nComputational cost: condition %d, n=%d, p=%d, reps %d--%d\n",
              condition_id, n_now, p_now, rep_start, rep_end))

  for (replication in seq.int(rep_start, rep_end)) {
    data_seed <- 1200000L + 10000L * condition_id + replication
    fit_seed <- 1300000L + 10000L * condition_id + replication

    tgen <- proc.time()[[3L]]
    sim <- simulate_student_rdmm(n = n_now, nu = nu_true,
                                 pars = pars, seed = data_seed)
    generation_seconds <- proc.time()[[3L]] - tgen

    # Reduce contamination from previous allocations before timing the fits.
    invisible(gc())
    fb <- fit_both_models(
      sim, r = c(5L, 2L), seed = fit_seed,
      max_iter = max_iter, eps = 1e-4,
      initial_nu = 6, nu_structure = "common",
      min_iter = min_iter, moving_window = moving_window
    )

    rr <- fb$results
    rr$n <- n_now
    rr$p <- p_now
    rr$Condition <- condition_id
    rr$Replication <- replication
    rr$Nu_true <- nu_true
    rr$Data_seed <- data_seed
    rr$Fit_seed <- fit_seed
    rr$Generation_seconds <- generation_seconds
    rr$Seconds_per_iteration <- ifelse(
      is.finite(rr$Iterations) & rr$Iterations > 0,
      rr$Elapsed_seconds / rr$Iterations,
      NA_real_
    )
    all_rows[[z]] <- rr
    z <- z + 1L
  }
}

out <- do.call(rbind, all_rows)
suffix <- sprintf("rep%04d_%04d", rep_start, rep_end)
if (cond_index > 0L) suffix <- paste0("cond", cond_index, "_", suffix)
write.csv(out,
          file.path(raw_dir, paste0("computational_cost_", suffix, ".csv")),
          row.names = FALSE)
cat("Computational-cost benchmark shard written.\n")
