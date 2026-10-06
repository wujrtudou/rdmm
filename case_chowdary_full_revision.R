############################################################
## case_chowdary_full_revision.R
## Full rerun of the Chowdary case study for the RDMM paper
##
## Scope: Chowdary case study ONLY.
## - Refit DGMM and RDMM over the exact manuscript grid
## - 10 random starts per candidate architecture
## - retain the start with the largest observed-data log-likelihood
## - select architecture by BIC using the mclust sign convention
##   (larger is better; reported BIC = - internal conventional BIC)
## - recompute ARI/MR only AFTER fitting/model selection
## - save Table 7 / Table 8 data and LaTeX rows
## - reconstruct posterior latent means for the selected RDMM
## - produce the requested Chowdary visualization
##
## Expected project layout:
## C:/Users/wujrt/Desktop/rdmm-main/
##   R/
##   scripts/
##   case_chowdary_full_revision.R
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

## ==========================================================
## 0. Paths and run controls
## ==========================================================

project_root <- "C:/Users/wujrt/Desktop/rdmm-main"
output_dir <- file.path(project_root, "results", "chowdary_revision")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## TRUE = reuse completed architecture fits from checkpoint files.
## FALSE = force every architecture to be refitted from scratch.
resume_from_checkpoint <- FALSE

## Exact manuscript settings: K1 fixed at 2; K2 in {2,4};
## 12 latent-dimension configurations shown in Table 8.
K1 <- 2L
K2_values <- c(2L, 4L)
r_pairs <- data.frame(
  r1 = c(4, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11),
  r2 = c(2, 3, 2, 4, 3, 5, 2, 4, 3, 5, 4, 5)
)

## Ten starts per architecture, matching the manuscript analysis.
seeds <- seq.int(10000L, 100000L, by = 10000L)

## Fitting settings used in the existing Chowdary analysis code.
fit_control <- list(
  it = 250L,
  eps = 1e-3,
  init = "kmeans",
  init_est = "factanal",
  psi_floor = 1e-6
)

rdmm_control <- c(
  fit_control,
  list(
    method = "sem",
    nu = 10,
    estimate_nu = TRUE,
    nu_structure = "common",
    nu_bounds = c(2.05, 200),
    min_iter = 10L,
    moving_window = 5L,
    verbose = FALSE
  )
)

cat("Project root:", project_root, "\n")
cat("Output directory:", output_dir, "\n\n")

if (!file.exists(file.path(project_root, "R", "load_models.R"))) {
  stop(
    "Cannot find R/load_models.R at: ", project_root,
    "\nCheck project_root at the top of this script."
  )
}

## ==========================================================
## 1. Packages and local model code
## ==========================================================

required_packages <- c("ICGE", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

source(file.path(project_root, "R", "load_models.R"))
load_rdmm_models(project_root)

if (!exists("robustdeepgmm", mode = "function")) stop("robustdeepgmm() not loaded.")
if (!exists("deepgmm", mode = "function")) stop("deepgmm() not loaded.")
if (!exists("adjustedRandIndex", mode = "function")) stop("adjustedRandIndex() not loaded.")
if (!exists(".rdmm_estep", mode = "function")) stop(".rdmm_estep() not loaded.")

## ==========================================================
## 2. Load and prepare Chowdary exactly once
## ==========================================================

data("chowdary", package = "ICGE", envir = environment())
if (!exists("chowdary")) stop("ICGE::chowdary could not be loaded.")

truth_original <- factor(as.character(unlist(chowdary[1L, , drop = TRUE])))
X_raw <- t(as.matrix(chowdary[-1L, , drop = FALSE]))
storage.mode(X_raw) <- "double"

if (nrow(X_raw) != length(truth_original)) {
  stop("Number of expression profiles and class labels do not match.")
}
if (anyNA(X_raw) || any(!is.finite(X_raw))) stop("X contains NA/NaN/Inf.")
if (anyNA(truth_original)) stop("Class labels contain missing values.")

## Remove only zero-variance genes, then standardize each retained gene once.
variable_sd <- apply(X_raw, 2L, sd)
keep_gene <- is.finite(variable_sd) & variable_sd > 0
if (!all(keep_gene)) {
  message("Removing ", sum(!keep_gene), " zero/non-finite-variance genes.")
}
X_raw <- X_raw[, keep_gene, drop = FALSE]
X <- scale(X_raw)
X <- as.matrix(X)
storage.mode(X) <- "double"

truth <- as.integer(truth_original)
n <- nrow(X)
p <- ncol(X)

cat("Dataset: ICGE::chowdary\n")
cat("n =", n, ", p =", p, "\n")
cat("Classes:\n")
print(table(truth_original))
cat("\n")

if (n != 104L) warning("Expected n=104, obtained n=", n)
if (p != 182L) warning("Expected p=182 after preprocessing, obtained p=", p)

## Exact candidate grid in the order used in the manuscript table.
model_grid <- do.call(
  rbind,
  lapply(K2_values, function(k2) {
    data.frame(k1 = K1, k2 = k2, r_pairs, pathways = K1 * k2)
  })
)
model_grid$model_id <- seq_len(nrow(model_grid))
model_grid <- model_grid[, c("model_id", "k1", "k2", "r1", "r2", "pathways")]

stopifnot(all(p > model_grid$r1), all(model_grid$r1 > model_grid$r2))
cat("Candidate architectures:", nrow(model_grid), "\n")
print(model_grid)
cat("\n")

## ==========================================================
## 3. Utility functions
## ==========================================================

safe_scalar <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || !is.finite(x[1L])) NA_real_ else x[1L]
}

fit_loglik <- function(fit) {
  for (nm in c("loglik", "lik", "likelihood", "logLik")) {
    if (!is.null(fit[[nm]])) {
      x <- suppressWarnings(as.numeric(fit[[nm]]))
      x <- x[is.finite(x)]
      if (length(x)) return(tail(x, 1L))
    }
  }
  NA_real_
}

fit_bic_internal <- function(fit) safe_scalar(fit$bic)
fit_bic_mclust <- function(fit) {
  ## Current project code uses conventional BIC = -2 logL + q log n.
  ## Manuscript reports the mclust sign convention, so negate it.
  b <- fit_bic_internal(fit)
  if (is.finite(b)) -b else NA_real_
}

first_layer_prediction <- function(fit, n_obs) {
  if (is.null(fit$s)) stop("fit$s is missing.")
  if (is.matrix(fit$s)) {
    if (nrow(fit$s) != n_obs || ncol(fit$s) < 1L) stop("Unexpected fit$s dimensions.")
    z <- fit$s[, 1L]
  } else {
    z <- fit$s
  }
  if (length(z) != n_obs || anyNA(z)) stop("Invalid first-layer assignments.")
  as.integer(as.factor(z))
}

## Exact MR for K1=2, but written generally for small numbers of classes.
all_perms <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], all_perms(x[-i]))))
}

misclassification_rate <- function(true_label, cluster_label) {
  true_label <- as.integer(as.factor(true_label))
  cluster_label <- as.integer(as.factor(cluster_label))
  tab <- table(cluster_label, true_label)

  ## Pad if needed, then brute-force permutation. Here it is only 2 x 2.
  m <- max(nrow(tab), ncol(tab))
  padded <- matrix(0, m, m)
  padded[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  perms <- all_perms(seq_len(m))
  correct_counts <- apply(perms, 1L, function(perm) {
    sum(padded[cbind(seq_len(m), perm)])
  })
  best <- perms[which.max(correct_counts), ]
  cluster_vals <- as.integer(rownames(tab))
  class_vals <- as.integer(colnames(tab))
  mapped <- rep(NA_integer_, nrow(tab))
  for (i in seq_len(nrow(tab))) {
    if (best[i] <= ncol(tab)) mapped[i] <- class_vals[best[i]]
  }
  pred_matched <- mapped[match(cluster_label, cluster_vals)]
  accuracy <- mean(pred_matched == true_label, na.rm = TRUE)
  list(
    MR = 1 - accuracy,
    Accuracy = accuracy,
    matched_prediction = pred_matched,
    mapping = data.frame(cluster = cluster_vals, class = mapped)
  )
}

score_fit <- function(fit, method, model_row, seed, elapsed) {
  pred <- first_layer_prediction(fit, n)
  mr <- misclassification_rate(truth, pred)
  data.frame(
    method = method,
    model_id = model_row$model_id,
    k1 = model_row$k1,
    k2 = model_row$k2,
    r1 = model_row$r1,
    r2 = model_row$r2,
    pathways = model_row$pathways,
    seed = seed,
    logLik = fit_loglik(fit),
    BIC_internal = fit_bic_internal(fit),
    BIC = fit_bic_mclust(fit),
    ARI = adjustedRandIndex(pred, truth),
    MR = mr$MR,
    Accuracy = mr$Accuracy,
    nu = if (method == "RDMM" && !is.null(fit$nu)) paste(format(as.numeric(fit$nu), digits = 8), collapse = ";") else NA_character_,
    iterations = if (!is.null(fit$iterations)) safe_scalar(fit$iterations) else NA_real_,
    elapsed_seconds = elapsed,
    status = "success"
  )
}

fit_one_start <- function(method, row, seed) {
  k_now <- c(row$k1, row$k2)
  r_now <- c(row$r1, row$r2)
  t0 <- proc.time()[["elapsed"]]

  fit <- tryCatch({
    if (method == "RDMM") {
      do.call(
        robustdeepgmm,
        c(
          list(y = X, layers = 2L, k = k_now, r = r_now, scale = FALSE, seed = seed),
          rdmm_control
        )
      )
    } else {
      do.call(
        deepgmm,
        c(
          list(y = X, layers = 2L, k = k_now, r = r_now, scale = FALSE, seed = seed),
          fit_control
        )
      )
    }
  }, error = function(e) structure(list(message = conditionMessage(e)), class = "fit_error"))

  elapsed <- proc.time()[["elapsed"]] - t0
  list(fit = fit, elapsed = elapsed)
}

## ==========================================================
## 4. Fit one complete method grid with checkpointing
## ==========================================================

run_method_grid <- function(method) {
  stopifnot(method %in% c("DGMM", "RDMM"))
  method_dir <- file.path(output_dir, tolower(method))
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

  all_start_rows <- list()
  best_rows <- list()
  best_fits <- vector("list", nrow(model_grid))

  for (i in seq_len(nrow(model_grid))) {
    row <- model_grid[i, , drop = FALSE]
    key <- sprintf("%s_K2_%d_r1_%d_r2_%d", tolower(method), row$k2, row$r1, row$r2)
    checkpoint_file <- file.path(method_dir, paste0(key, ".rds"))

    cat("\n============================================================\n")
    cat(method, "architecture", i, "of", nrow(model_grid), "\n")
    cat(sprintf("K=(%d,%d), r=(%d,%d)\n", row$k1, row$k2, row$r1, row$r2))

    if (resume_from_checkpoint && file.exists(checkpoint_file)) {
      obj <- readRDS(checkpoint_file)
      cat("Loaded checkpoint:", basename(checkpoint_file), "\n")
    } else {
      start_rows <- list()
      successful_fits <- list()

      for (j in seq_along(seeds)) {
        seed_now <- seeds[j]
        cat("  start", j, "of", length(seeds), "seed", seed_now, "...")
        ans <- fit_one_start(method, row, seed_now)

        if (inherits(ans$fit, "fit_error")) {
          cat(" FAILED\n")
          start_rows[[j]] <- data.frame(
            method = method, model_id = row$model_id, k1 = row$k1, k2 = row$k2,
            r1 = row$r1, r2 = row$r2, pathways = row$pathways, seed = seed_now,
            logLik = NA_real_, BIC_internal = NA_real_, BIC = NA_real_,
            ARI = NA_real_, MR = NA_real_, Accuracy = NA_real_, nu = NA_character_,
            iterations = NA_real_, elapsed_seconds = ans$elapsed,
            status = paste0("failed: ", ans$fit$message)
          )
          next
        }

        sc <- score_fit(ans$fit, method, row, seed_now, ans$elapsed)
        start_rows[[j]] <- sc
        successful_fits[[as.character(seed_now)]] <- ans$fit
        cat(sprintf(" logLik=%.4f BIC=%.2f ARI=%.4f MR=%.4f\n",
                    sc$logLik, sc$BIC, sc$ARI, sc$MR))
      }

      starts_df <- do.call(rbind, start_rows)
      ok <- which(starts_df$status == "success" & is.finite(starts_df$logLik))
      if (!length(ok)) {
        obj <- list(starts = starts_df, best_row = NULL, best_fit = NULL)
      } else {
        ## Manuscript rule: retain the start with largest observed-data log-likelihood.
        best_idx <- ok[which.max(starts_df$logLik[ok])]
        best_seed <- starts_df$seed[best_idx]
        best_fit <- successful_fits[[as.character(best_seed)]]
        best_row <- starts_df[best_idx, , drop = FALSE]
        obj <- list(starts = starts_df, best_row = best_row, best_fit = best_fit)
      }
      saveRDS(obj, checkpoint_file)
    }

    all_start_rows[[i]] <- obj$starts
    if (!is.null(obj$best_row)) best_rows[[i]] <- obj$best_row
    best_fits[[i]] <- obj$best_fit
  }

  starts <- do.call(rbind, all_start_rows)
  best <- do.call(rbind, best_rows)
  rownames(starts) <- NULL
  rownames(best) <- NULL

  ## Architecture selection: mclust sign convention, larger BIC preferred.
  best$BIC_rank <- NA_integer_
  ok <- which(best$status == "success" & is.finite(best$BIC))
  ord <- ok[order(best$BIC[ok], decreasing = TRUE)]
  best$BIC_rank[ord] <- seq_along(ord)
  selected_idx <- ord[1L]
  selected_row <- best[selected_idx, , drop = FALSE]
  selected_model_id <- selected_row$model_id
  selected_fit <- best_fits[[selected_model_id]]

  write.csv(starts, file.path(output_dir, paste0("chowdary_", tolower(method), "_all_starts.csv")), row.names = FALSE)
  write.csv(best, file.path(output_dir, paste0("chowdary_", tolower(method), "_architecture_results.csv")), row.names = FALSE)
  saveRDS(best_fits, file.path(output_dir, paste0("chowdary_", tolower(method), "_best_fits_by_architecture.rds")))
  saveRDS(selected_fit, file.path(output_dir, paste0("chowdary_", tolower(method), "_bic_selected_fit.rds")))

  cat("\n", method, "BIC-selected architecture:\n", sep = "")
  print(selected_row)

  list(starts = starts, best = best, best_fits = best_fits,
       selected_row = selected_row, selected_fit = selected_fit)
}

## ==========================================================
## 5. Run ALL Chowdary DGMM and RDMM experiments
## ==========================================================

cat("\n\n#################### DGMM ####################\n")
dgmm_out <- run_method_grid("DGMM")

cat("\n\n#################### RDMM ####################\n")
rdmm_out <- run_method_grid("RDMM")

## ==========================================================
## 6. Rebuild manuscript tables from the NEW results
## ==========================================================

## Table 7: selected models plus published benchmark values.
benchmarks <- data.frame(
  Method = c("MCLUST", "MFA", "MCFA"),
  Architecture = c("VVI", "q=3", "q=1"),
  ARI = c(0.0657, 0.5858, 0.6800),
  MR = c(0.3462, 0.1154, 0.0865),
  stringsAsFactors = FALSE
)

selected_table <- rbind(
  benchmarks,
  data.frame(
    Method = "DGMM",
    Architecture = sprintf("(K2,r1,r2)=(%d,%d,%d)",
                           dgmm_out$selected_row$k2,
                           dgmm_out$selected_row$r1,
                           dgmm_out$selected_row$r2),
    ARI = dgmm_out$selected_row$ARI,
    MR = dgmm_out$selected_row$MR
  ),
  data.frame(
    Method = "RDMM",
    Architecture = sprintf("(K2,r1,r2)=(%d,%d,%d)",
                           rdmm_out$selected_row$k2,
                           rdmm_out$selected_row$r1,
                           rdmm_out$selected_row$r2),
    ARI = rdmm_out$selected_row$ARI,
    MR = rdmm_out$selected_row$MR
  )
)
write.csv(selected_table, file.path(output_dir, "Table7_chowdary_selected_models.csv"), row.names = FALSE)

## Table 8: matched architecture comparison.
dcols <- c("k2", "r1", "r2", "pathways", "BIC", "ARI", "MR", "logLik", "best_seed")
rcols <- c("k2", "r1", "r2", "pathways", "BIC", "ARI", "MR", "logLik", "best_seed", "nu")

d <- dgmm_out$best
r <- rdmm_out$best
## Rename 'seed' of retained start to best_seed for table clarity.
d$best_seed <- d$seed
r$best_seed <- r$seed

d2 <- d[, dcols]
r2 <- r[, rcols]
names(d2)[names(d2) %in% c("BIC", "ARI", "MR", "logLik", "best_seed")] <- paste0("DGMM_", names(d2)[names(d2) %in% c("BIC", "ARI", "MR", "logLik", "best_seed")])
names(r2)[names(r2) %in% c("BIC", "ARI", "MR", "logLik", "best_seed", "nu")] <- paste0("RDMM_", names(r2)[names(r2) %in% c("BIC", "ARI", "MR", "logLik", "best_seed", "nu")])

matched <- merge(d2, r2, by = c("k2", "r1", "r2", "pathways"), sort = FALSE)
matched <- merge(model_grid[, c("model_id", "k2", "r1", "r2")], matched,
                 by = c("k2", "r1", "r2"), sort = FALSE)
matched <- matched[order(matched$model_id), ]
matched$model_id <- NULL
write.csv(matched, file.path(output_dir, "Table8_chowdary_matched_architectures.csv"), row.names = FALSE)

## Selected RDMM nu summary.
nu_selected <- as.numeric(rdmm_out$selected_fit$nu)
writeLines(
  c(
    sprintf("Selected DGMM: K2=%d, r1=%d, r2=%d, BIC=%.6f, ARI=%.6f, MR=%.6f",
            dgmm_out$selected_row$k2, dgmm_out$selected_row$r1, dgmm_out$selected_row$r2,
            dgmm_out$selected_row$BIC, dgmm_out$selected_row$ARI, dgmm_out$selected_row$MR),
    sprintf("Selected RDMM: K2=%d, r1=%d, r2=%d, BIC=%.6f, ARI=%.6f, MR=%.6f",
            rdmm_out$selected_row$k2, rdmm_out$selected_row$r1, rdmm_out$selected_row$r2,
            rdmm_out$selected_row$BIC, rdmm_out$selected_row$ARI, rdmm_out$selected_row$MR),
    paste0("Selected RDMM nu: ", paste(format(nu_selected, digits = 8), collapse = ", "))
  ),
  file.path(output_dir, "chowdary_selected_model_summary.txt")
)

## ==========================================================
## 7. Export LaTeX rows from the newly generated numbers
## ==========================================================

fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)

latex7 <- c(
  "% Auto-generated by case_chowdary_full_revision.R",
  paste0("MCLUST & VVI & 0.0657 & 0.3462 \\\\"),
  paste0("MFA & $q=3$ & 0.5858 & 0.1154 \\\\"),
  paste0("MCFA & $q=1$ & 0.6800 & 0.0865 \\\\"),
  sprintf("DGMM & $(K_2,r_1,r_2)=(%d,%d,%d)$ & %s & %s \\\\",
          dgmm_out$selected_row$k2, dgmm_out$selected_row$r1, dgmm_out$selected_row$r2,
          fmt(dgmm_out$selected_row$ARI), fmt(dgmm_out$selected_row$MR)),
  sprintf("RDMM & $(K_2,r_1,r_2)=(%d,%d,%d)$ & \\textbf{%s} & \\textbf{%s} \\\\",
          rdmm_out$selected_row$k2, rdmm_out$selected_row$r1, rdmm_out$selected_row$r2,
          fmt(rdmm_out$selected_row$ARI), fmt(rdmm_out$selected_row$MR))
)
writeLines(latex7, file.path(output_dir, "Table7_chowdary_selected_models_rows.tex"))

latex8 <- c("% Auto-generated by case_chowdary_full_revision.R")
for (i in seq_len(nrow(matched))) {
  z <- matched[i, ]
  latex8 <- c(latex8, sprintf(
    "%d & %d & %d & %d & %.2f & %.4f & %.4f & %.2f & %.4f & %.4f \\\\",
    z$k2, z$r1, z$r2, z$pathways,
    z$DGMM_BIC, z$DGMM_ARI, z$DGMM_MR,
    z$RDMM_BIC, z$RDMM_ARI, z$RDMM_MR
  ))
}
writeLines(latex8, file.path(output_dir, "Table8_chowdary_matched_rows.tex"))

## ==========================================================
## 8. Reconstruct posterior latent means for selected RDMM
## ==========================================================

best_fit <- rdmm_out$selected_fit
best_arch <- rdmm_out$selected_row

## Strict dimension sanity checks for the selected fit.
expected_paths <- best_arch$k1 * best_arch$k2
if (nrow(best_fit$ps.y) != n || ncol(best_fit$ps.y) != expected_paths) {
  stop("Selected RDMM ps.y dimensions do not match n x (K1*K2).")
}
if (length(best_fit$ps.y.list) != 2L ||
    ncol(best_fit$ps.y.list[[1]]) != best_arch$k1 ||
    ncol(best_fit$ps.y.list[[2]]) != best_arch$k2) {
  stop("Selected RDMM layer posterior dimensions are inconsistent with K1/K2.")
}

## The fit was given the already-standardized X with scale=FALSE.
## Re-evaluate the exact conditional posterior quantities at the retained fit.
r_full <- c(p, best_fit$r)
estep <- .rdmm_estep(
  y = X,
  H = best_fit$H,
  mu = best_fit$mu,
  psi = best_fit$psi,
  w = best_fit$w,
  nu_path = best_fit$nu_path,
  paths = best_fit$paths,
  r_full = r_full,
  floor = fit_control$psi_floor
)


## Posterior mean of all latent variables, averaged over pathways:
## E[z_j | y_j] = sum_s tau_js E[z_j | y_j, S_j=s].
q <- sum(best_fit$r)
Z_post <- matrix(0, nrow = n, ncol = q)
for (s in seq_len(nrow(best_fit$paths))) {
  Z_post <- Z_post + estep$conditional_mean[[s]] * estep$tau[, s]
}

r1 <- best_fit$r[1L]
r2 <- best_fit$r[2L]
Z1_hat <- Z_post[, seq_len(r1), drop = FALSE]
Z2_hat <- Z_post[, r1 + seq_len(r2), drop = FALSE]

## Sanity checks: sample-level posterior means should not collapse to a few rows.
cat("\nPosterior latent means:\n")
cat("Z1_hat dimension:", paste(dim(Z1_hat), collapse = " x "),
    "; unique rows:", nrow(unique(as.data.frame(round(Z1_hat, 10)))), "\n")
cat("Z2_hat dimension:", paste(dim(Z2_hat), collapse = " x "),
    "; unique rows:", nrow(unique(as.data.frame(round(Z2_hat, 10)))), "\n")

## Posterior PCA coordinates. One PCA is fitted once for Z1 and reused in panels a-d.
pca_z1 <- prcomp(Z1_hat, center = TRUE, scale. = FALSE)
coord_z1 <- pca_z1$x[, 1:2, drop = FALSE]
var_z1 <- pca_z1$sdev^2 / sum(pca_z1$sdev^2)

pca_z2 <- prcomp(Z2_hat, center = TRUE, scale. = FALSE)
coord_z2 <- pca_z2$x[, 1:2, drop = FALSE]
var_z2 <- pca_z2$sdev^2 / sum(pca_z2$sdev^2)

first_prob <- best_fit$ps.y.list[[1L]]
second_prob <- best_fit$ps.y.list[[2L]]
first_cluster <- max.col(first_prob, ties.method = "first")
second_cluster <- max.col(second_prob, ties.method = "first")
first_conf <- apply(first_prob, 1L, max)
second_conf <- apply(second_prob, 1L, max)
path_cluster <- max.col(best_fit$ps.y, ties.method = "first")
path_conf <- apply(best_fit$ps.y, 1L, max)
path_entropy <- -rowSums(best_fit$ps.y * log(pmax(best_fit$ps.y, .Machine$double.xmin)))

## posterior_w is already sum_s tau_js E(W_js | y_j,S_j=s).
posterior_w <- as.numeric(estep$posterior_w)
if (length(posterior_w) != n || any(!is.finite(posterior_w)) || any(posterior_w <= 0)) {
  stop("Invalid posterior_w in selected RDMM fit.")
}
log_inv_w <- log10(1 / posterior_w)

viz <- data.frame(
  sample = seq_len(n),
  tissue = truth_original,
  z1_PC1 = coord_z1[, 1L],
  z1_PC2 = coord_z1[, 2L],
  z2_PC1 = coord_z2[, 1L],
  z2_PC2 = coord_z2[, 2L],
  first_layer = factor(first_cluster),
  first_confidence = first_conf,
  second_layer = factor(second_cluster),
  second_confidence = second_conf,
  pathway = factor(path_cluster),
  pathway_confidence = path_conf,
  pathway_entropy = path_entropy,
  posterior_w = posterior_w,
  log10_inverse_precision = log_inv_w
)
write.csv(viz, file.path(output_dir, "chowdary_selected_rdmm_visualization_data.csv"), row.names = FALSE)

save(
  best_fit, best_arch, estep, Z1_hat, Z2_hat, viz,
  truth_original, X, X_raw,
  file = file.path(output_dir, "chowdary_selected_rdmm_complete.RData")
)

## ==========================================================
## 9. Main requested visualization: fitted latent representation
## ==========================================================

library(ggplot2)
library(patchwork)

theme_case <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(size = 10)
  )

xlab_z1 <- sprintf("PC1 of posterior mean z^(1) (%.1f%%)", 100 * var_z1[1L])
ylab_z1 <- sprintf("PC2 of posterior mean z^(1) (%.1f%%)", 100 * var_z1[2L])

## (a) Primary RDMM partition. Alpha conveys layer-1 posterior confidence.
p_a <- ggplot(viz, aes(z1_PC1, z1_PC2, colour = first_layer, alpha = first_confidence)) +
  geom_point(size = 2.8) +
  scale_alpha_continuous(range = c(0.35, 1), limits = c(0, 1)) +
  labs(
    title = "(a) First-layer RDMM partition",
    subtitle = "Opacity = first-layer posterior confidence",
    x = xlab_z1, y = ylab_z1,
    colour = "First layer", alpha = "Posterior\nconfidence"
  ) + theme_case

## (b) External tissue labels; same fitted latent coordinates.
p_b <- ggplot(viz, aes(z1_PC1, z1_PC2, colour = tissue)) +
  geom_point(size = 2.8, alpha = 0.85) +
  labs(
    title = "(b) Known tissue class",
    subtitle = "Used only for post-selection evaluation",
    x = xlab_z1, y = ylab_z1,
    colour = "Tissue"
  ) + theme_case

## (c) Robustness mechanism on the same fitted latent coordinates.
## Larger log10(1/w_bar) means smaller posterior precision and stronger down-weighting.
p_c <- ggplot(viz, aes(z1_PC1, z1_PC2, colour = log10_inverse_precision)) +
  geom_point(size = 2.9, alpha = 0.9) +
  scale_colour_viridis_c(option = "C") +
  labs(
    title = "(c) Posterior robust down-weighting",
    subtitle = "Larger values indicate stronger down-weighting",
    x = xlab_z1, y = ylab_z1,
    colour = expression(log[10](1/bar(w)[j]))
  ) + theme_case

## (d) Second-layer structure on the SAME first-latent representation.
## This makes it visually explicit that layer 2 refines residual heterogeneity.
p_d <- ggplot(viz, aes(z1_PC1, z1_PC2, colour = second_layer, alpha = second_confidence)) +
  geom_point(size = 2.8) +
  scale_alpha_continuous(range = c(0.35, 1), limits = c(0, 1)) +
  labs(
    title = "(d) Second-layer latent structure",
    subtitle = "Opacity = second-layer posterior confidence",
    x = xlab_z1, y = ylab_z1,
    colour = "Second layer", alpha = "Posterior\nconfidence"
  ) + theme_case

main_figure <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = "Latent representation and posterior diagnostics for the BIC-selected RDMM",
    subtitle = sprintf(
      "Chowdary data; K=(%d,%d), r=(%d,%d); all %d samples shown",
      best_arch$k1, best_arch$k2, best_arch$r1, best_arch$r2, n
    )
  )

print(main_figure)
ggsave(file.path(output_dir, "Figure_Chowdary_RDMM_latent_posterior.pdf"),
       main_figure, width = 10.8, height = 8.4, units = "in")
ggsave(file.path(output_dir, "Figure_Chowdary_RDMM_latent_posterior.png"),
       main_figure, width = 10.8, height = 8.4, units = "in", dpi = 600)

## ==========================================================
## 10. Additional useful visualization: complete-pathway heatmap
## ==========================================================

path_labels <- apply(best_fit$paths, 1L, function(z) paste0("(", paste(z, collapse = ","), ")"))
order_samples <- order(viz$first_layer, viz$second_layer, viz$tissue, -viz$pathway_confidence)
rank_in_plot <- integer(n)
rank_in_plot[order_samples] <- seq_len(n)

heat <- data.frame(
  sample_rank = rep(rank_in_plot, times = ncol(best_fit$ps.y)),
  pathway = factor(rep(path_labels, each = n), levels = path_labels),
  responsibility = as.vector(best_fit$ps.y)
)

p_heat <- ggplot(heat, aes(pathway, sample_rank, fill = responsibility)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", limits = c(0, 1)) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(
    title = "Complete-pathway posterior responsibilities",
    subtitle = "Samples ordered by fitted first layer, second layer, and tissue class",
    x = expression("Complete pathway " * (S[1] * "," * S[2])),
    y = "Ordered samples",
    fill = expression(tau[js])
  ) +
  theme_case +
  theme(panel.grid = element_blank())

print(p_heat)
ggsave(file.path(output_dir, "Figure_Chowdary_RDMM_pathway_heatmap.pdf"),
       p_heat, width = 7.2, height = 7.0, units = "in")
ggsave(file.path(output_dir, "Figure_Chowdary_RDMM_pathway_heatmap.png"),
       p_heat, width = 7.2, height = 7.0, units = "in", dpi = 600)

## ==========================================================
## 11. Supplementary observed-data PCA: visual motivation for robustness
## ==========================================================

pca_obs <- prcomp(X, center = FALSE, scale. = FALSE)
obs_xy <- pca_obs$x[, 1:2, drop = FALSE]
obs_var <- pca_obs$sdev^2 / sum(pca_obs$sdev^2)
obs_dat <- data.frame(
  PC1 = obs_xy[, 1L], PC2 = obs_xy[, 2L],
  tissue = truth_original,
  first_layer = factor(first_cluster),
  log10_inverse_precision = log_inv_w
)

p_obs1 <- ggplot(obs_dat, aes(PC1, PC2, colour = tissue)) +
  geom_point(size = 2.6, alpha = 0.85) +
  labs(
    title = "(a) Observed-data PCA: full range",
    x = sprintf("PC1 (%.1f%%)", 100 * obs_var[1L]),
    y = sprintf("PC2 (%.1f%%)", 100 * obs_var[2L]),
    colour = "Tissue"
  ) + theme_case

## Automatic central zoom based on 5th-95th percentiles; no observations are removed.
xlim_zoom <- as.numeric(quantile(obs_dat$PC1, c(0.05, 0.95), names = FALSE))
ylim_zoom <- as.numeric(quantile(obs_dat$PC2, c(0.05, 0.95), names = FALSE))

p_obs2 <- ggplot(obs_dat, aes(PC1, PC2, colour = first_layer)) +
  geom_point(size = 2.6, alpha = 0.85) +
  coord_cartesian(xlim = xlim_zoom, ylim = ylim_zoom) +
  labs(
    title = "(b) Same PCA: central-region zoom",
    subtitle = "Plotting window only; no samples removed",
    x = sprintf("PC1 (%.1f%%)", 100 * obs_var[1L]),
    y = sprintf("PC2 (%.1f%%)", 100 * obs_var[2L]),
    colour = "First layer"
  ) + theme_case

obs_figure <- p_obs1 | p_obs2
print(obs_figure)
ggsave(file.path(output_dir, "Figure_Chowdary_observed_PCA.pdf"),
       obs_figure, width = 10.5, height = 4.8, units = "in")
ggsave(file.path(output_dir, "Figure_Chowdary_observed_PCA.png"),
       obs_figure, width = 10.5, height = 4.8, units = "in", dpi = 600)

## ==========================================================
## 12. Final diagnostic report
## ==========================================================

cat("\n\n############################################################\n")
cat("CHOWDARY CASE STUDY RERUN COMPLETE\n")
cat("############################################################\n")
cat("\nSelected DGMM:\n")
print(dgmm_out$selected_row)
cat("\nSelected RDMM:\n")
print(rdmm_out$selected_row)
cat("\nSelected RDMM degrees of freedom:\n")
print(best_fit$nu)
cat("\nSelected RDMM posterior dimensions:\n")
cat("ps.y:", paste(dim(best_fit$ps.y), collapse = " x "), "\n")
cat("layer 1 posterior:", paste(dim(best_fit$ps.y.list[[1]]), collapse = " x "), "\n")
cat("layer 2 posterior:", paste(dim(best_fit$ps.y.list[[2]]), collapse = " x "), "\n")
cat("Z1 posterior mean:", paste(dim(Z1_hat), collapse = " x "), "\n")
cat("Z2 posterior mean:", paste(dim(Z2_hat), collapse = " x "), "\n")
cat("\nOutput files are in:\n", output_dir, "\n", sep = "")
cat("\nMain manuscript outputs:\n")
cat("  Table7_chowdary_selected_models.csv\n")
cat("  Table8_chowdary_matched_architectures.csv\n")
cat("  Table7_chowdary_selected_models_rows.tex\n")
cat("  Table8_chowdary_matched_rows.tex\n")
cat("  Figure_Chowdary_RDMM_latent_posterior.pdf/png\n")
cat("\nAdditional diagnostics:\n")
cat("  Figure_Chowdary_RDMM_pathway_heatmap.pdf/png\n")
cat("  Figure_Chowdary_observed_PCA.pdf/png\n")
cat("  chowdary_selected_rdmm_complete.RData\n")
cat("  chowdary_*_all_starts.csv\n")
cat("  chowdary_*_architecture_results.csv\n")
