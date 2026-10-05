# Reviewer-requested visual diagnostics from simulated data:
# (i) second-layer latent coordinates and observed-space PCs;
# (ii) multiple-start MCEM/SEM traces, especially degrees of freedom.

project_root <- Sys.getenv("PROJECT_ROOT", unset = ".")
source(file.path(project_root, "R", "load_models.R"))
load_rdmm_models(project_root)
source(file.path(project_root, "R", "simulation_parameters.R"))
source(file.path(project_root, "R", "simulation_helpers.R"))

pars <- make_manuscript_parameters()
sim <- simulate_student_rdmm(n = 1000L, nu = 4, pars = pars, seed = 20261002L)
pc <- prcomp(sim$y, center = TRUE, scale. = TRUE)$x[, 1:2, drop = FALSE]

viz <- data.frame(
  Observation = seq_len(nrow(sim$y)),
  Z2_1 = sim$z2[, 1L], Z2_2 = sim$z2[, 2L],
  PC1 = pc[, 1L], PC2 = pc[, 2L],
  Layer1 = factor(sim$layer1), Layer2 = factor(sim$layer2),
  Pathway = factor(sim$path_id)
)
raw_dir <- file.path(project_root, "results", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(viz, file.path(raw_dir, "simulated_visualization_data.csv"), row.names = FALSE)

starts <- c(101L, 203L, 307L, 409L, 503L, 607L)
traces <- lapply(starts, function(s) {
  tr <- fit_rdmm_trace(
    sim$y, r = c(5L, 2L), seed = s, it = 300L, eps = 1e-4,
    initial_nu = 6, min_iter = 100L, moving_window = 20L
  )
  tr$Start <- factor(s)
  tr
})
traces <- do.call(rbind, traces)
write.csv(traces, file.path(raw_dir, "rdmm_multiple_start_traces.csv"), row.names = FALSE)
cat("Visualization and convergence raw data written.\n")
