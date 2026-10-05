project_root <- Sys.getenv("PROJECT_ROOT", unset = ".")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
source(file.path(project_root, "R", "output_helpers.R"))

raw_dir <- file.path(project_root, "results", "raw")
tab_dir <- file.path(project_root, "results", "tables")
fig_dir <- file.path(project_root, "results", "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =====================
# Simulation 1: clustering + nu recovery
# =====================
sim1 <- read_csv_shards(raw_dir, "^sim1_fit_.*\\.csv$")
sim1 <- dedupe_rows(sim1, c("Nu_true", "Replication", "Model"))
sim1$Model <- factor(sim1$Model, levels = c("DGMM", "RDMM"))

nu_vals <- sort(unique(sim1$Nu_true))
tab1 <- list(); z <- 1L
for (nu in nu_vals) {
  for (model in c("DGMM", "RDMM")) {
    d <- sim1[sim1$Nu_true == nu & sim1$Model == model & sim1$Success %in% TRUE, , drop = FALSE]
    tab1[[z]] <- data.frame(
      Nu = sprintf("%.1f", nu), Method = model,
      Estimated_nu = if (model == "RDMM") mean_sd_string(d$Estimated_nu) else "--",
      First_layer_ARI = mean_sd_string(d$FirstLayer_ARI),
      First_layer_MR = mean_sd_string(d$FirstLayer_MR),
      Pathway_ARI = mean_sd_string(d$Pathway_ARI),
      Pathway_MR = mean_sd_string(d$Pathway_MR),
      stringsAsFactors = FALSE
    )
    z <- z + 1L
  }
}
tab1 <- do.call(rbind, tab1)
write.csv(tab1, file.path(tab_dir, "Table_Sim1_clustering_nu.csv"), row.names = FALSE)
write_latex_table(tab1, file.path(tab_dir, "Table_Sim1_clustering_nu.tex"),
  "Recovery of the common degrees of freedom and clustering performance in Simulation 1. Entries are Monte Carlo means with standard deviations in parentheses.",
  "tab:study1-results", align = "c l c c c c c")

# Clustering boxplots (same role as current simulation1.png)
cl_long <- long_metrics(sim1[sim1$Success %in% TRUE, ],
                        c("Nu_true", "Model", "Replication"),
                        c("FirstLayer_ARI", "FirstLayer_MR", "Pathway_ARI", "Pathway_MR"))
cl_long$Metric <- factor(cl_long$Metric,
  levels = c("FirstLayer_ARI", "FirstLayer_MR", "Pathway_ARI", "Pathway_MR"),
  labels = c("First-layer ARI", "First-layer MR", "Pathway ARI", "Pathway MR"))
p <- ggplot2::ggplot(cl_long, ggplot2::aes(x = factor(Nu_true), y = Value, fill = Model)) +
  ggplot2::geom_boxplot(outlier.size = 0.3, position = ggplot2::position_dodge(width = 0.8)) +
  ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  ggplot2::labs(x = "True degrees of freedom", y = NULL, fill = "Method") +
  ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
save_plot_both(p, file.path(fig_dir, "Figure_Sim1_clustering"), 9, 6.5)

# =====================
# Simulation 1: other parameter recovery
# =====================
rec <- read_csv_shards(raw_dir, "^sim1_recovery_.*\\.csv$")
rec <- dedupe_rows(rec, c("Nu_true", "Replication"))
metrics_rec <- c("Pi1_RMSE", "Pi2_RMSE", "Eta1_RMSE", "Psi1_RMSE",
                 "Loading1_Subspace_Error", "PathMean_RMSE", "PathScale_RelFrob")
rec_sum <- summary_long(rec, "Nu_true", metrics_rec)
write.csv(rec_sum, file.path(tab_dir, "Table_Sim1_parameter_recovery_long.csv"), row.names = FALSE)

# Compact table, mean (SD).
tab2 <- data.frame(Nu = sprintf("%.1f", nu_vals), stringsAsFactors = FALSE)
for (m in metrics_rec) {
  tab2[[m]] <- vapply(nu_vals, function(nu) mean_sd_string(rec[[m]][rec$Nu_true == nu]), character(1))
}
write.csv(tab2, file.path(tab_dir, "Table_Sim1_parameter_recovery.csv"), row.names = FALSE)
write_latex_table(tab2, file.path(tab_dir, "Table_Sim1_parameter_recovery.tex"),
  "Recovery of additional RDMM parameters in Simulation 1. Entries are Monte Carlo means with standard deviations in parentheses. PathScale RelFrob is the relative Frobenius error of the implied observed-level pathway scale matrix; the loading error is rotation-invariant.",
  "tab:study1-parameter-recovery")

# Recovery metric figure.
rec_plot <- rec_sum
rec_plot$Metric <- factor(rec_plot$Metric, levels = metrics_rec,
  labels = c("Layer-1 mixing weight RMSE", "Layer-2 mixing weight RMSE",
             "Layer-1 location RMSE", "Layer-1 residual variance RMSE",
             "Layer-1 loading subspace error", "Pathway location RMSE",
             "Pathway scale relative Frobenius error"))
p <- ggplot2::ggplot(rec_plot,
  ggplot2::aes(x = Nu_true, y = Mean)) +
  ggplot2::geom_line() + ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - 1.96 * SE, ymax = Mean + 1.96 * SE), width = 0.08) +
  ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  ggplot2::labs(x = "True degrees of freedom", y = "Recovery error (mean ± 1.96 SE)") +
  ggplot2::theme_bw()
save_plot_both(p, file.path(fig_dir, "Figure_Sim1_parameter_recovery"), 10, 10)

# Degrees-of-freedom recovery as its own concise figure.
rd_only <- sim1[sim1$Model == "RDMM" & sim1$Success %in% TRUE, ]
p <- ggplot2::ggplot(rd_only, ggplot2::aes(x = factor(Nu_true), y = Estimated_nu)) +
  ggplot2::geom_boxplot(outlier.size = 0.4) +
  ggplot2::geom_point(data = data.frame(Nu_true = factor(nu_vals), Estimated_nu = nu_vals),
                      ggplot2::aes(x = Nu_true, y = Estimated_nu), shape = 4, size = 3,
                      inherit.aes = FALSE) +
  ggplot2::labs(x = "True degrees of freedom", y = "Estimated degrees of freedom") +
  ggplot2::theme_bw()
save_plot_both(p, file.path(fig_dir, "Figure_Sim1_nu_recovery"), 7.5, 5.5)

# =====================
# Latent dimension misspecification
# =====================
mis <- read_csv_shards(raw_dir, "^dimension_misspec_.*\\.csv$")
mis <- dedupe_rows(mis, c("Architecture", "Replication", "Model"))
arch_levels <- c("r1 under", "r2 under", "true", "r1 over", "r2 over")
mis$Architecture <- factor(mis$Architecture, levels = arch_levels)
mis$Model <- factor(mis$Model, levels = c("DGMM", "RDMM"))

metrics_mis <- c("Estimated_nu", "FirstLayer_ARI", "FirstLayer_MR", "Pathway_ARI", "Pathway_MR")
tab3 <- list(); z <- 1L
for (a in arch_levels) {
  for (model in c("DGMM", "RDMM")) {
    d <- mis[mis$Architecture == a & mis$Model == model & mis$Success %in% TRUE, , drop = FALSE]
    if (!nrow(d)) next
    tab3[[z]] <- data.frame(
      Architecture = a,
      r1 = unique(d$r1_fit)[1L], r2 = unique(d$r2_fit)[1L], Method = model,
      Estimated_nu = if (model == "RDMM") mean_sd_string(d$Estimated_nu) else "--",
      First_layer_ARI = mean_sd_string(d$FirstLayer_ARI),
      First_layer_MR = mean_sd_string(d$FirstLayer_MR),
      Pathway_ARI = mean_sd_string(d$Pathway_ARI),
      Pathway_MR = mean_sd_string(d$Pathway_MR),
      stringsAsFactors = FALSE
    )
    z <- z + 1L
  }
}
tab3 <- do.call(rbind, tab3)
write.csv(tab3, file.path(tab_dir, "Table_dimension_misspecification.csv"), row.names = FALSE)
write_latex_table(tab3, file.path(tab_dir, "Table_dimension_misspecification.tex"),
  "Sensitivity to slight misspecification of the latent dimensions. Data are generated with (r1,r2)=(5,2) and nu=4. Entries are Monte Carlo means with standard deviations in parentheses.",
  "tab:dimension-sensitivity")

mis_sum <- summary_long(mis, c("Architecture", "Model"),
                        c("FirstLayer_ARI", "FirstLayer_MR", "Pathway_ARI", "Pathway_MR"),
                        success_col = "Success")
mis_sum$Metric <- factor(mis_sum$Metric,
  levels = c("FirstLayer_ARI", "FirstLayer_MR", "Pathway_ARI", "Pathway_MR"),
  labels = c("First-layer ARI", "First-layer MR", "Pathway ARI", "Pathway MR"))
p <- ggplot2::ggplot(mis_sum,
  ggplot2::aes(x = Architecture, y = Mean, group = Model, shape = Model, linetype = Model)) +
  ggplot2::geom_line() + ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - 1.96 * SE, ymax = Mean + 1.96 * SE), width = 0.12) +
  ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  ggplot2::labs(x = "Fitted latent dimensions", y = "Mean ± 1.96 SE", shape = "Method", linetype = "Method") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), legend.position = "bottom")
save_plot_both(p, file.path(fig_dir, "Figure_dimension_misspecification"), 9, 7)

# =====================
# Gaussian clean base + contamination
# =====================
gc <- read_csv_shards(raw_dir, "^gaussian_contamination_.*\\.csv$")
gc <- dedupe_rows(gc, c("Contamination", "Replication", "Model"))
gc$Model <- factor(gc$Model, levels = c("DGMM", "RDMM"))
eps_vals <- sort(unique(gc$Contamination))

tab4 <- list(); z <- 1L
for (eps in eps_vals) {
  for (model in c("DGMM", "RDMM")) {
    d <- gc[gc$Contamination == eps & gc$Model == model & gc$Success %in% TRUE, , drop = FALSE]
    tab4[[z]] <- data.frame(
      Contamination = sprintf("%.2f", eps), Method = model,
      Estimated_nu = if (model == "RDMM") mean_sd_string(d$Estimated_nu) else "--",
      First_layer_ARI = mean_sd_string(d$FirstLayer_ARI),
      First_layer_MR = mean_sd_string(d$FirstLayer_MR),
      Pathway_ARI = mean_sd_string(d$Pathway_ARI),
      Pathway_MR = mean_sd_string(d$Pathway_MR), stringsAsFactors = FALSE
    )
    z <- z + 1L
  }
}
tab4 <- do.call(rbind, tab4)
write.csv(tab4, file.path(tab_dir, "Table_Gaussian_contamination.csv"), row.names = FALSE)
write_latex_table(tab4, file.path(tab_dir, "Table_Gaussian_contamination.tex"),
  "Clustering performance when the uncontaminated data are generated from a Gaussian DGMM and increasingly contaminated by multivariate t perturbations. Entries are Monte Carlo means with standard deviations in parentheses.",
  "tab:gaussian-contamination", align = "c l c c c c c")

gc_long <- long_metrics(gc[gc$Success %in% TRUE, ],
                         c("Contamination", "Model", "Replication"),
                         c("FirstLayer_ARI", "FirstLayer_MR", "Pathway_ARI", "Pathway_MR"))
gc_long$Metric <- factor(gc_long$Metric,
  levels = c("FirstLayer_ARI", "FirstLayer_MR", "Pathway_ARI", "Pathway_MR"),
  labels = c("First-layer ARI", "First-layer MR", "Pathway ARI", "Pathway MR"))
p <- ggplot2::ggplot(gc_long,
  ggplot2::aes(x = factor(Contamination), y = Value, fill = Model)) +
  ggplot2::geom_boxplot(outlier.size = 0.3, position = ggplot2::position_dodge(width = 0.8)) +
  ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  ggplot2::labs(x = "Contamination proportion", y = NULL, fill = "Method") +
  ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
save_plot_both(p, file.path(fig_dir, "Figure_Gaussian_contamination_clustering"), 9, 6.5)

rdgc <- gc[gc$Model == "RDMM" & gc$Success %in% TRUE, ]
p <- ggplot2::ggplot(rdgc,
  ggplot2::aes(x = factor(Contamination), y = Estimated_nu)) +
  ggplot2::geom_boxplot(outlier.size = 0.35) +
  ggplot2::labs(x = "Contamination proportion", y = "Estimated common degrees of freedom") +
  ggplot2::theme_bw()
save_plot_both(p, file.path(fig_dir, "Figure_Gaussian_contamination_nu"), 7.5, 5.5)

# =====================
# Simulated data visualization
# =====================
viz_file <- file.path(raw_dir, "simulated_visualization_data.csv")
if (file.exists(viz_file)) {
  vz <- read.csv(viz_file, stringsAsFactors = FALSE)
  zz <- rbind(
    data.frame(X = vz$Z2_1, Y = vz$Z2_2, Layer1 = factor(vz$Layer1),
               Pathway = factor(vz$Pathway), Space = "Deepest latent coordinates (z2)"),
    data.frame(X = vz$PC1, Y = vz$PC2, Layer1 = factor(vz$Layer1),
               Pathway = factor(vz$Pathway), Space = "First two PCs of observed data")
  )
  p <- ggplot2::ggplot(zz, ggplot2::aes(x = X, y = Y, shape = Layer1)) +
    ggplot2::geom_point(alpha = 0.55, size = 1.3) +
    ggplot2::facet_wrap(~Space, scales = "free", ncol = 2) +
    ggplot2::labs(x = NULL, y = NULL, shape = "First-layer component") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  save_plot_both(p, file.path(fig_dir, "Figure_simulated_data_visualization"), 10, 4.8)
}

# =====================
# Multiple-start convergence traces
# =====================
trace_file <- file.path(raw_dir, "rdmm_multiple_start_traces.csv")
if (file.exists(trace_file)) {
  tr <- read.csv(trace_file, stringsAsFactors = FALSE)
  tr_long <- rbind(
    data.frame(Iteration = tr$Iteration, Start = factor(tr$Start),
               Metric = "Observed-data log-likelihood", Value = tr$LogLik),
    data.frame(Iteration = tr$Iteration, Start = factor(tr$Start),
               Metric = "Estimated degrees of freedom", Value = tr$Nu)
  )
  p <- ggplot2::ggplot(tr_long,
    ggplot2::aes(x = Iteration, y = Value, group = Start, linetype = Start)) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "SEM iteration", y = NULL, linetype = "Random start") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  save_plot_both(p, file.path(fig_dir, "Figure_multiple_start_convergence"), 8.5, 7)

  last_rows <- do.call(rbind, lapply(split(tr, tr$Start), function(d) d[nrow(d), ]))
  write.csv(last_rows, file.path(tab_dir, "Table_multiple_start_convergence.csv"), row.names = FALSE)
}

# Failure-rate audit table: important for transparent simulation reporting.
audit <- rbind(
  data.frame(Experiment = "Simulation 1", Model = levels(sim1$Model),
             Failure_rate = vapply(levels(sim1$Model), function(m) mean(!(sim1$Success[sim1$Model == m] %in% TRUE)), numeric(1))),
  data.frame(Experiment = "Dimension sensitivity", Model = levels(mis$Model),
             Failure_rate = vapply(levels(mis$Model), function(m) mean(!(mis$Success[mis$Model == m] %in% TRUE)), numeric(1))),
  data.frame(Experiment = "Gaussian contamination", Model = levels(gc$Model),
             Failure_rate = vapply(levels(gc$Model), function(m) mean(!(gc$Success[gc$Model == m] %in% TRUE)), numeric(1)))
)
write.csv(audit, file.path(tab_dir, "simulation_failure_rates.csv"), row.names = FALSE)

# =====================
# Residual-covariance misspecification sensitivity
# =====================
cov_files <- list.files(raw_dir, pattern = "^covariance_sensitivity_fit_.*\\.csv$", full.names = TRUE)
if (length(cov_files)) {
  covs <- read_csv_shards(raw_dir, "^covariance_sensitivity_fit_.*\\.csv$")
  covs <- dedupe_rows(covs, c("Residual_rho", "Replication", "Model"))
  covs$Model <- factor(covs$Model, levels = c("DGMM", "RDMM"))
  rho_vals <- sort(unique(covs$Residual_rho))

  tab_cov <- list(); z <- 1L
  for (rho in rho_vals) {
    for (model in c("DGMM", "RDMM")) {
      d <- covs[covs$Residual_rho == rho & covs$Model == model & covs$Success %in% TRUE, , drop = FALSE]
      if (!nrow(d)) next
      tab_cov[[z]] <- data.frame(
        Residual_correlation = sprintf("%.1f", rho), Method = model,
        Estimated_nu = if (model == "RDMM") mean_sd_string(d$Estimated_nu) else "--",
        First_layer_ARI = mean_sd_string(d$FirstLayer_ARI),
        First_layer_MR = mean_sd_string(d$FirstLayer_MR),
        Pathway_ARI = mean_sd_string(d$Pathway_ARI),
        Pathway_MR = mean_sd_string(d$Pathway_MR),
        stringsAsFactors = FALSE
      )
      z <- z + 1L
    }
  }
  tab_cov <- do.call(rbind, tab_cov)
  write.csv(tab_cov, file.path(tab_dir, "Table_covariance_sensitivity.csv"), row.names = FALSE)
  write_latex_table(tab_cov, file.path(tab_dir, "Table_covariance_sensitivity.tex"),
    "Sensitivity to violation of the diagonal residual-covariance assumption. Data are generated with AR(1) residual correlation rho at both layers while preserving the manuscript marginal residual variances; both fitted models retain diagonal residual covariance matrices. Entries are Monte Carlo means with standard deviations in parentheses.",
    "tab:covariance-sensitivity", align = "c l c c c c c")

  cov_long <- long_metrics(covs[covs$Success %in% TRUE, ],
                           c("Residual_rho", "Model", "Replication"),
                           c("FirstLayer_ARI", "FirstLayer_MR", "Pathway_ARI", "Pathway_MR"))
  cov_long$Metric <- factor(cov_long$Metric,
    levels = c("FirstLayer_ARI", "FirstLayer_MR", "Pathway_ARI", "Pathway_MR"),
    labels = c("First-layer ARI", "First-layer MR", "Pathway ARI", "Pathway MR"))
  p <- ggplot2::ggplot(cov_long,
    ggplot2::aes(x = Residual_rho, y = Value, group = Model, shape = Model, linetype = Model)) +
    ggplot2::stat_summary(fun = mean, geom = "line") +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 2) +
    ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = "True AR(1) residual correlation (rho)", y = "Monte Carlo mean",
                  shape = "Method", linetype = "Method") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  save_plot_both(p, file.path(fig_dir, "Figure_covariance_sensitivity_clustering"), 9, 7)

  rec_files <- list.files(raw_dir, pattern = "^covariance_sensitivity_recovery_.*\\.csv$", full.names = TRUE)
  if (length(rec_files)) {
    covrec <- read_csv_shards(raw_dir, "^covariance_sensitivity_recovery_.*\\.csv$")
    covrec <- dedupe_rows(covrec, c("Residual_rho", "Replication"))
    rec_tab <- data.frame(Residual_correlation = sprintf("%.1f", rho_vals), stringsAsFactors = FALSE)
    rec_tab$Estimated_nu <- vapply(rho_vals, function(rho) mean_sd_string(covrec$Nu_hat[covrec$Residual_rho == rho]), character(1))
    rec_tab$PathMean_RMSE <- vapply(rho_vals, function(rho) mean_sd_string(covrec$PathMean_RMSE[covrec$Residual_rho == rho]), character(1))
    rec_tab$PathScale_RelFrob <- vapply(rho_vals, function(rho) mean_sd_string(covrec$PathScale_RelFrob[covrec$Residual_rho == rho]), character(1))
    write.csv(rec_tab, file.path(tab_dir, "Table_covariance_sensitivity_recovery.csv"), row.names = FALSE)
    write_latex_table(rec_tab, file.path(tab_dir, "Table_covariance_sensitivity_recovery.tex"),
      "RDMM recovery under increasing off-diagonal residual correlation. PathScale RelFrob measures the relative Frobenius error of the implied observed-level pathway scale matrix.",
      "tab:covariance-sensitivity-recovery")
  }
}

# =====================
# Empirical computational cost / scalability
# =====================
cost_files <- list.files(raw_dir, pattern = "^computational_cost_.*\\.csv$", full.names = TRUE)
if (length(cost_files)) {
  cost <- read_csv_shards(raw_dir, "^computational_cost_.*\\.csv$")
  cost <- dedupe_rows(cost, c("n", "p", "Replication", "Model"))
  cost$Model <- factor(cost$Model, levels = c("DGMM", "RDMM"))

  conds <- unique(cost[c("n", "p")])
  conds <- conds[order(conds$p, conds$n), , drop = FALSE]
  tab_cost <- list(); z <- 1L
  for (i in seq_len(nrow(conds))) {
    n_now <- conds$n[i]
    p_now <- conds$p[i]
    for (model in c("DGMM", "RDMM")) {
      d_all <- cost[cost$n == n_now & cost$p == p_now & cost$Model == model, , drop = FALSE]
      d <- d_all[d_all$Success %in% TRUE, , drop = FALSE]
      tab_cost[[z]] <- data.frame(
        n = n_now, p = p_now, Method = model,
        Runtime_seconds = mean_sd_string(d$Elapsed_seconds, digits = 2),
        Iterations = mean_sd_string(d$Iterations, digits = 1),
        Seconds_per_iteration = mean_sd_string(d$Seconds_per_iteration, digits = 3),
        Success_rate = sprintf("%.3f", mean(d_all$Success %in% TRUE)),
        stringsAsFactors = FALSE
      )
      z <- z + 1L
    }
  }
  tab_cost <- do.call(rbind, tab_cost)
  write.csv(tab_cost, file.path(tab_dir, "Table_computational_cost.csv"), row.names = FALSE)
  write_latex_table(tab_cost, file.path(tab_dir, "Table_computational_cost.tex"),
    "Empirical computational cost as sample size n and observed dimension p are varied while the fitted architecture is held fixed at K=(2,2) and r=(5,2). Timings are wall-clock seconds from sequential fits; entries are means with standard deviations in parentheses. Seconds per iteration separates per-iteration computational cost from differences in convergence iteration counts.",
    "tab:computational-cost", align = "c c l c c c c")

  base_p <- min(cost$p)
  # The p-scan is identified by the n value having multiple distinct p values.
  p_counts <- tapply(cost$p, cost$n, function(x) length(unique(x)))
  base_n <- as.numeric(names(p_counts)[which.max(p_counts)])
  # The n-scan is identified by the p value having multiple distinct n values.
  n_counts <- tapply(cost$n, cost$p, function(x) length(unique(x)))
  base_p <- as.numeric(names(n_counts)[which.max(n_counts)])

  nscan <- cost[cost$p == base_p & cost$Success %in% TRUE, , drop = FALSE]
  nscan$Scaling <- sprintf("Vary n (p = %d)", base_p)
  nscan$Size <- nscan$n
  pscan <- cost[cost$n == base_n & cost$Success %in% TRUE, , drop = FALSE]
  pscan$Scaling <- sprintf("Vary p (n = %d)", base_n)
  pscan$Size <- pscan$p
  cost_plot <- rbind(nscan, pscan)

  if (nrow(cost_plot)) {
    p <- ggplot2::ggplot(cost_plot,
      ggplot2::aes(x = Size, y = Elapsed_seconds, group = Model,
                   shape = Model, linetype = Model)) +
      ggplot2::stat_summary(fun = mean, geom = "line") +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 2) +
      ggplot2::facet_wrap(~Scaling, scales = "free_x", ncol = 2) +
      ggplot2::labs(x = "Problem size (see panel title)", y = "Mean wall-clock fit time (seconds)",
                    shape = "Method", linetype = "Method") +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
    save_plot_both(p, file.path(fig_dir, "Figure_computational_cost"), 9, 4.8)
  }
}

# Update the failure-rate audit to include the two revision experiments when run.
if (exists("covs")) {
  cov_models <- levels(covs$Model)
  audit <- rbind(audit,
    data.frame(Experiment = "Covariance sensitivity", Model = cov_models,
               Failure_rate = vapply(cov_models, function(m)
                 mean(!(covs$Success[covs$Model == m] %in% TRUE)), numeric(1))))
}
if (exists("cost")) {
  cost_models <- levels(cost$Model)
  audit <- rbind(audit,
    data.frame(Experiment = "Computational cost", Model = cost_models,
               Failure_rate = vapply(cost_models, function(m)
                 mean(!(cost$Success[cost$Model == m] %in% TRUE)), numeric(1))))
}
write.csv(audit, file.path(tab_dir, "simulation_failure_rates.csv"), row.names = FALSE)
cat("All tables and figures generated under results/tables and results/figures.\n")

