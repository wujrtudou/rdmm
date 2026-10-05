# Fitting, evaluation and parameter-recovery helpers.

all_permutations <- function(x) {
  x <- as.integer(x)
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], all_permutations(x[-i]))
  }))
}

best_mapping <- function(predicted, truth, groups) {
  predicted <- as.integer(predicted)
  truth <- as.integer(truth)
  perms <- all_permutations(seq_len(groups))
  err <- apply(perms, 1L, function(mp) mean(mp[predicted] != truth))
  j <- which.min(err)
  list(mapping = as.integer(perms[j, ]), error = err[j])
}

adjusted_rand_index <- function(x, y) {
  x <- as.vector(x); y <- as.vector(y)
  tab <- table(x, y)
  c2 <- function(z) z * (z - 1) / 2
  a <- sum(c2(tab))
  rp <- sum(c2(rowSums(tab)))
  cp <- sum(c2(colSums(tab)))
  tot <- c2(sum(tab))
  if (tot <= 0) return(0)
  expected <- rp * cp / tot
  max_index <- 0.5 * (rp + cp)
  denom <- max_index - expected
  if (abs(denom) < .Machine$double.eps) return(ifelse(a == max_index, 1, 0))
  (a - expected) / denom
}

misclassification_rate <- function(predicted, truth, groups) {
  best_mapping(predicted, truth, groups)$error
}

fit_dgmm_model <- function(y, r = c(5L, 2L), seed = 1L,
                           max_iter = 300L, eps = 1e-4) {
  set.seed(as.integer(seed))
  deepgmm(
    y = y, layers = 2L, k = c(2L, 2L), r = as.integer(r),
    it = as.integer(max_iter), eps = eps,
    init = "kmeans", init_est = "factanal",
    seed = as.integer(seed), scale = FALSE, psi_floor = 1e-6
  )
}

fit_rdmm_model <- function(y, r = c(5L, 2L), seed = 1L,
                           max_iter = 300L, eps = 1e-4,
                           initial_nu = 6,
                           nu_structure = "common",
                           min_iter = 100L, moving_window = 20L,
                           nu_bounds = c(2.05, 200)) {
  set.seed(as.integer(seed))
  robustdeepgmm(
    y = y, layers = 2L, k = c(2L, 2L), r = as.integer(r),
    it = as.integer(max_iter), eps = eps,
    init = "kmeans", init_est = "factanal",
    seed = as.integer(seed), scale = FALSE,
    nu = initial_nu, nu_structure = nu_structure, estimate_nu = TRUE,
    method = "sem", psi_floor = 1e-6, nu_bounds = nu_bounds,
    min_iter = as.integer(min_iter), moving_window = as.integer(moving_window),
    verbose = FALSE
  )
}

extract_path_prediction <- function(fit) {
  if (!is.null(fit$ps.y)) return(max.col(as.matrix(fit$ps.y), ties.method = "first"))
  encode_path_id(fit$s, fit$k)
}

clustering_metrics <- function(fit, sim) {
  pred1 <- as.integer(fit$s[, 1L])
  pred_path <- extract_path_prediction(fit)
  data.frame(
    FirstLayer_ARI = adjusted_rand_index(pred1, sim$layer1),
    FirstLayer_MR = misclassification_rate(pred1, sim$layer1, sim$parameters$k[1L]),
    Pathway_ARI = adjusted_rand_index(pred_path, sim$path_id),
    Pathway_MR = misclassification_rate(pred_path, sim$path_id, prod(sim$parameters$k)),
    stringsAsFactors = FALSE
  )
}

fit_summary_row <- function(fit, model, sim, elapsed) {
  cm <- clustering_metrics(fit, sim)
  nuhat <- if (model == "RDMM") mean(as.numeric(fit$nu)) else NA_real_
  loglik <- if (model == "RDMM") as.numeric(fit$loglik) else tail(as.numeric(fit$lik), 1L)
  iterations <- if (!is.null(fit$iterations)) {
    as.integer(fit$iterations)
  } else if (!is.null(fit$lik)) {
    # DGMM stores the full likelihood history but not a separate iteration count.
    max(0L, length(fit$lik) - 1L)
  } else {
    NA_integer_
  }
  data.frame(
    Model = model,
    Estimated_nu = nuhat,
    FirstLayer_ARI = cm$FirstLayer_ARI,
    FirstLayer_MR = cm$FirstLayer_MR,
    Pathway_ARI = cm$Pathway_ARI,
    Pathway_MR = cm$Pathway_MR,
    LogLik = loglik,
    BIC = if (!is.null(fit$bic)) as.numeric(fit$bic) else NA_real_,
    Iterations = iterations,
    Elapsed_seconds = elapsed,
    Success = TRUE,
    Error = "",
    stringsAsFactors = FALSE
  )
}

failed_fit_row <- function(model, msg, elapsed = NA_real_) {
  data.frame(
    Model = model, Estimated_nu = NA_real_,
    FirstLayer_ARI = NA_real_, FirstLayer_MR = NA_real_,
    Pathway_ARI = NA_real_, Pathway_MR = NA_real_,
    LogLik = NA_real_, BIC = NA_real_, Iterations = NA_integer_,
    Elapsed_seconds = elapsed, Success = FALSE, Error = as.character(msg),
    stringsAsFactors = FALSE
  )
}

fit_both_models <- function(sim, r = c(5L, 2L), seed = 1L,
                            max_iter = 300L, eps = 1e-4,
                            initial_nu = 6, nu_structure = "common",
                            min_iter = 100L, moving_window = 20L) {
  rows <- list(); fits <- list()

  t0 <- proc.time()[[3L]]
  dg <- tryCatch(fit_dgmm_model(sim$y, r, seed, max_iter, eps), error = identity)
  elapsed <- proc.time()[[3L]] - t0
  if (inherits(dg, "error")) {
    rows$DGMM <- failed_fit_row("DGMM", conditionMessage(dg), elapsed)
  } else {
    fits$DGMM <- dg
    rows$DGMM <- fit_summary_row(dg, "DGMM", sim, elapsed)
  }

  t0 <- proc.time()[[3L]]
  rd <- tryCatch(
    fit_rdmm_model(sim$y, r, seed, max_iter, eps, initial_nu,
                   nu_structure, min_iter, moving_window),
    error = identity
  )
  elapsed <- proc.time()[[3L]] - t0
  if (inherits(rd, "error")) {
    rows$RDMM <- failed_fit_row("RDMM", conditionMessage(rd), elapsed)
  } else {
    fits$RDMM <- rd
    rows$RDMM <- fit_summary_row(rd, "RDMM", sim, elapsed)
  }

  list(results = do.call(rbind, rows), fits = fits)
}

# ----- Identifiable / rotation-invariant parameter recovery -----

pred_order_for_true_labels <- function(mapping) {
  # mapping[predicted] = true. Return predicted index for true labels 1,...,K.
  vapply(seq_along(mapping), function(t) which(mapping == t)[1L], integer(1))
}

array_component_matrix <- function(x, component) {
  matrix(x[component, , ], nrow = dim(x)[2L], ncol = dim(x)[3L])
}

projection_matrix <- function(A, tol = 1e-10) {
  A <- as.matrix(A)
  qrA <- qr(A, tol = tol)
  q <- qr.Q(qrA)
  rank <- qrA$rank
  if (rank == 0L) return(matrix(0, nrow(A), nrow(A)))
  Q <- q[, seq_len(rank), drop = FALSE]
  Q %*% t(Q)
}

subspace_error <- function(Ahat, Atrue) {
  Ph <- projection_matrix(Ahat)
  Pt <- projection_matrix(Atrue)
  sqrt(sum((Ph - Pt)^2)) / sqrt(2 * max(1L, ncol(Atrue)))
}

collapse_true_path <- function(pars, a, b) {
  L1 <- pars$lambda[[1L]][[a]]
  L2 <- pars$lambda[[2L]][[b]]
  eta1 <- pars$eta[[1L]][, a]
  eta2 <- pars$eta[[2L]][, b]
  P1 <- pars$psi[[1L]][[a]]
  P2 <- pars$psi[[2L]][[b]]

  mu_z1 <- eta2
  Sigma_z1 <- P2 + L2 %*% t(L2)
  mu_y <- as.numeric(eta1 + L1 %*% mu_z1)
  Sigma_y <- P1 + L1 %*% Sigma_z1 %*% t(L1)
  list(mu = mu_y, scale = (Sigma_y + t(Sigma_y)) / 2)
}

collapse_fit_path <- function(fit, a_pred, b_pred) {
  L1 <- array_component_matrix(fit$H[[1L]], a_pred)
  L2 <- array_component_matrix(fit$H[[2L]], b_pred)
  eta1 <- fit$mu[[1L]][, a_pred]
  eta2 <- fit$mu[[2L]][, b_pred]
  P1 <- array_component_matrix(fit$psi[[1L]], a_pred)
  P2 <- array_component_matrix(fit$psi[[2L]], b_pred)

  mu_z1 <- eta2
  Sigma_z1 <- P2 + L2 %*% t(L2)
  mu_y <- as.numeric(eta1 + L1 %*% mu_z1)
  Sigma_y <- P1 + L1 %*% Sigma_z1 %*% t(L1)
  list(mu = mu_y, scale = (Sigma_y + t(Sigma_y)) / 2)
}

parameter_recovery_metrics <- function(fit, sim) {
  if (is.null(fit$s) || is.null(fit$H) || is.null(fit$mu) || is.null(fit$psi)) {
    stop("Fit object does not contain parameters required for recovery metrics.")
  }
  pars <- sim$parameters

  map1 <- best_mapping(fit$s[, 1L], sim$layer1, pars$k[1L])$mapping
  map2 <- best_mapping(fit$s[, 2L], sim$layer2, pars$k[2L])$mapping
  ord1 <- pred_order_for_true_labels(map1)
  ord2 <- pred_order_for_true_labels(map2)

  w1 <- fit$w[[1L]][ord1]
  w2 <- fit$w[[2L]][ord2]
  pi1_rmse <- sqrt(mean((w1 - pars$pi[[1L]])^2))
  pi2_rmse <- sqrt(mean((w2 - pars$pi[[2L]])^2))

  # First-layer locations live in observed coordinates and can be compared directly.
  eta1_hat <- fit$mu[[1L]][, ord1, drop = FALSE]
  eta1_rmse <- sqrt(mean((eta1_hat - pars$eta[[1L]])^2))

  # Diagonal first-layer specific variances.
  psi1_hat <- do.call(cbind, lapply(ord1, function(a) diag(array_component_matrix(fit$psi[[1L]], a))))
  psi1_true <- do.call(cbind, lapply(pars$psi[[1L]], diag))
  psi1_rmse <- sqrt(mean((psi1_hat - psi1_true)^2))

  # Rotation-invariant first-layer loading subspace error.
  loading1_err <- mean(vapply(seq_along(ord1), function(a_true) {
    a_pred <- ord1[a_true]
    subspace_error(array_component_matrix(fit$H[[1L]], a_pred),
                   pars$lambda[[1L]][[a_true]])
  }, numeric(1)))

  # Fully identifiable observed-level pathway location and scale recovery.
  path_mu_sq <- numeric(0)
  path_scale_rel <- numeric(0)
  for (a_true in seq_len(pars$k[1L])) {
    for (b_true in seq_len(pars$k[2L])) {
      tru <- collapse_true_path(pars, a_true, b_true)
      est <- collapse_fit_path(fit, ord1[a_true], ord2[b_true])
      path_mu_sq <- c(path_mu_sq, mean((est$mu - tru$mu)^2))
      denom <- sqrt(sum(tru$scale^2))
      path_scale_rel <- c(path_scale_rel,
                          sqrt(sum((est$scale - tru$scale)^2)) / max(denom, 1e-12))
    }
  }

  data.frame(
    Nu_hat = mean(as.numeric(fit$nu)),
    Nu_error = mean(as.numeric(fit$nu)) - sim$nu,
    Pi1_RMSE = pi1_rmse,
    Pi2_RMSE = pi2_rmse,
    Eta1_RMSE = eta1_rmse,
    Psi1_RMSE = psi1_rmse,
    Loading1_Subspace_Error = loading1_err,
    PathMean_RMSE = sqrt(mean(path_mu_sq)),
    PathScale_RelFrob = mean(path_scale_rel),
    stringsAsFactors = FALSE
  )
}

# ----- RDMM trace for reviewer-requested convergence plots -----
# This repeats the same internal SEM updates as robustdeepgmm(), but records
# parameter histories at every iteration. It is used only for diagnostics/figures.
fit_rdmm_trace <- function(y, r = c(5L, 2L), seed = 1L,
                           it = 300L, eps = 1e-4, initial_nu = 6,
                           min_iter = 100L, moving_window = 20L,
                           psi_floor = 1e-6, nu_bounds = c(2.05, 200)) {
  set.seed(as.integer(seed))
  y <- as.matrix(y)
  k <- c(2L, 2L)
  r_full <- c(ncol(y), as.integer(r))
  initial <- .rdmm_initialise(y, 2L, k, r_full, "kmeans", "factanal", psi_floor)

  H <- initial$H; mu <- initial$mu; psi <- initial$psi
  psi_inv <- initial$psi.inv; w <- initial$w
  paths <- .rdmm_paths(k)
  nu_path <- .rdmm_expand_nu(initial_nu, "common", paths, k)

  trace <- vector("list", as.integer(it))
  likelihood <- numeric(0)
  ratio <- Inf

  for (iteration in seq_len(as.integer(it))) {
    estep <- .rdmm_estep(y, H, mu, psi, w, nu_path, paths, r_full, psi_floor)
    draws <- .rdmm_draw_missing(estep, psi_floor)
    upd <- .rdmm_mstep_sem(y, estep, draws, k, r_full, psi_floor)
    new_nu <- .rdmm_update_nu(estep, nu_path, "common", k, nu_bounds)
    new_estep <- .rdmm_estep(y, upd$H, upd$mu, upd$psi, upd$w,
                             new_nu, paths, r_full, psi_floor)

    likelihood <- c(likelihood, new_estep$loglik)
    H <- upd$H; mu <- upd$mu; psi <- upd$psi; psi_inv <- upd$psi.inv; w <- upd$w
    nu_path <- new_nu

    trace[[iteration]] <- data.frame(
      Iteration = iteration,
      LogLik = new_estep$loglik,
      Nu = mean(new_nu),
      Pi11 = w[[1L]][1L],
      Eta111 = mu[[1L]][1L, 1L],
      stringsAsFactors = FALSE
    )

    if (length(likelihood) >= 2L * moving_window && iteration >= min_iter) {
      old_ma <- mean(likelihood[(length(likelihood) - 2L * moving_window + 1L):
                                  (length(likelihood) - moving_window)])
      new_ma <- mean(tail(likelihood, moving_window))
      ratio <- abs(new_ma - old_ma) / (abs(old_ma) + 1)
    }
    if (iteration >= min_iter && is.finite(ratio) && ratio < eps) break
  }

  do.call(rbind, trace[seq_len(iteration)])
}

read_env_int <- function(name, default) {
  z <- Sys.getenv(name, unset = "")
  if (!nzchar(z)) return(as.integer(default))
  as.integer(z)
}

read_env_num <- function(name, default) {
  z <- Sys.getenv(name, unset = "")
  if (!nzchar(z)) return(as.numeric(default))
  as.numeric(z)
}

read_env_int_vector <- function(name, default) {
  z <- Sys.getenv(name, unset = "")
  if (!nzchar(z)) return(as.integer(default))
  vals <- trimws(strsplit(z, ",", fixed = TRUE)[[1L]])
  vals <- as.integer(vals[nzchar(vals)])
  if (!length(vals) || anyNA(vals)) stop(name, " must be a comma-separated integer vector.")
  vals
}

run_parallel <- function(X, FUN, cores = 1L,
                         project_root = Sys.getenv("PROJECT_ROOT", unset = ".")) {
  cores <- max(1L, as.integer(cores))

  # Serial fallback.
  if (cores <= 1L || length(X) <= 1L) {
    return(lapply(X, FUN))
  }

  # Windows: use a PSOCK cluster because mclapply() is not available.
  if (.Platform$OS.type == "windows") {
    project_root <- normalizePath(
      project_root, winslash = "/", mustWork = TRUE
    )

    cl <- parallel::makeCluster(cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    cat(sprintf("Windows PSOCK cluster started with %d workers.\n", cores))

    # Load model and simulation code in every worker.
    init_status <- parallel::clusterCall(
      cl,
      function(root) {
        source(file.path(root, "R", "load_models.R"), local = .GlobalEnv)
        load_rdmm_models(root)

        source(
          file.path(root, "R", "simulation_parameters.R"),
          local = .GlobalEnv
        )
        source(
          file.path(root, "R", "simulation_helpers.R"),
          local = .GlobalEnv
        )

        TRUE
      },
      project_root
    )

    if (!all(vapply(init_status, isTRUE, logical(1)))) {
      stop("At least one Windows worker failed during initialization.")
    }

    # Export condition-specific objects used by the worker function
    # (e.g. nu, pars, dd, nu_true, eps_grid, max_iter).
    fun_env <- environment(FUN)
    free_vars <- codetools::findGlobals(FUN, merge = FALSE)$variables
    export_vars <- free_vars[
      vapply(
        free_vars,
        exists,
        logical(1),
        envir = fun_env,
        inherits = TRUE
      )
    ]

    if (length(export_vars)) {
      parallel::clusterExport(
        cl,
        varlist = unique(export_vars),
        envir = fun_env
      )
    }

    return(parallel::parLapply(cl, X, FUN))
  }

  # Linux/macOS: use forked workers.
  parallel::mclapply(
    X,
    FUN,
    mc.cores = cores,
    mc.preschedule = FALSE
  )
}
