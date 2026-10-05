# Sanity checks before manuscript use.
project_root <- Sys.getenv("PROJECT_ROOT", unset = ".")
source(file.path(project_root, "R", "output_helpers.R"))
raw_dir <- file.path(project_root, "results", "raw")

sim1 <- dedupe_rows(read_csv_shards(raw_dir, "^sim1_fit_.*\\.csv$"),
                    c("Nu_true", "Replication", "Model"))
mis <- dedupe_rows(read_csv_shards(raw_dir, "^dimension_misspec_.*\\.csv$"),
                   c("Architecture", "Replication", "Model"))
gc <- dedupe_rows(read_csv_shards(raw_dir, "^gaussian_contamination_.*\\.csv$"),
                  c("Contamination", "Replication", "Model"))

cat("Simulation 1 rows:", nrow(sim1), "\n")
print(with(sim1, table(Nu_true, Model)))
cat("\nDimension sensitivity rows:\n")
print(with(mis, table(Architecture, Model)))
cat("\nGaussian contamination rows:\n")
print(with(gc, table(Contamination, Model)))

cat("\nFailure rates:\n")
print(aggregate(!Success ~ Model, sim1, mean))
print(aggregate(!Success ~ Architecture + Model, mis, mean))
print(aggregate(!Success ~ Contamination + Model, gc, mean))

# Expected paper-scale counts if N_REP=500 and the full design was run.
expected_sim1 <- 6L * 500L * 2L
expected_mis <- 5L * 500L * 2L
expected_gc <- 6L * 500L * 2L
cat(sprintf("\nExpected full-run rows: Sim1=%d, dimension=%d, Gaussian contamination=%d\n",
            expected_sim1, expected_mis, expected_gc))

if (nrow(sim1) < expected_sim1 || nrow(mis) < expected_mis || nrow(gc) < expected_gc) {
  warning("At least one experiment has fewer rows than the 500-replication full design.")
}

# Optional revision experiments.
cov_files <- list.files(raw_dir, pattern = "^covariance_sensitivity_fit_.*\\.csv$", full.names = TRUE)
if (length(cov_files)) {
  covs <- dedupe_rows(read_csv_shards(raw_dir, "^covariance_sensitivity_fit_.*\\.csv$"),
                      c("Residual_rho", "Replication", "Model"))
  cat("\nCovariance-sensitivity rows:\n")
  print(with(covs, table(Residual_rho, Model)))
  cat("Covariance-sensitivity failure rates:\n")
  print(aggregate(!Success ~ Residual_rho + Model, covs, mean))
}

cost_files <- list.files(raw_dir, pattern = "^computational_cost_.*\\.csv$", full.names = TRUE)
if (length(cost_files)) {
  cost <- dedupe_rows(read_csv_shards(raw_dir, "^computational_cost_.*\\.csv$"),
                      c("n", "p", "Replication", "Model"))
  cat("\nComputational-cost rows:\n")
  print(with(cost, table(n, p, Model)))
  cat("Computational-cost failure rates:\n")
  print(aggregate(!Success ~ n + p + Model, cost, mean))
  cat("Mean fit times (seconds):\n")
  print(aggregate(Elapsed_seconds ~ n + p + Model,
                  cost[cost$Success %in% TRUE, ], mean))
}

