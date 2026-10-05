deepgmm <- function(y, layers, k, r,
            it = 250, eps = 0.001, init = "kmeans", init_est = "factanal",
            seed = NULL, scale = TRUE, psi_floor = 1e-6) {

  # DGMM now uses exactly the same high-dimensional initialization engine as
  # RDMM and DStMM.  The shared initializer provides:
  #   * identical initial clustering rules;
  #   * minimum cluster-size repair;
  #   * factanal only when the within-group covariance can be full rank;
  #   * SVD/PCA fallback when p >= n_g or factanal fails;
  #   * diagonal residual variances bounded below by psi_floor.
  if (!exists(".rdmm_initialise", mode = "function", inherits = TRUE) ||
      !exists(".rdmm_fix_names", mode = "function", inherits = TRUE)) {
    stop("Shared DGMM/RDMM/DStMM initialization helpers are not loaded. ",
         "Load the project with load_dstmm_project(PROJECT_ROOT).")
  }

  if (is.data.frame(y)) y <- as.matrix(y)
  if (!is.matrix(y) || !is.numeric(y)) stop("y must be a numeric matrix.")
  if (anyNA(y) || any(!is.finite(y))) {
    stop("y contains missing or non-finite values.")
  }
  if (!is.finite(psi_floor) || length(psi_floor) != 1L || psi_floor <= 0) {
    stop("psi_floor must be one positive finite number.")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("The value of seed must be one finite integer-like number.")
    }
    set.seed(as.integer(seed))
  }

  if (scale) {
    ys <- base::scale(y)
    sc <- attr(ys, "scaled:scale")
    if (any(!is.finite(sc) | sc == 0)) {
      stop("Cannot scale y because at least one variable has zero/non-finite scale.")
    }
    y <- as.matrix(ys)
  }

  # Use the same canonical naming rules as RDMM/DStMM.
  fixed <- .rdmm_fix_names(init, init_est)
  init <- fixed$init
  init_est <- fixed$init_est

  # Retain the original DGMM argument checks for compatibility with the SEM
  # implementation, then construct the same r_full object used by RDMM/DStMM.
  valid_args(Y = y, layers = layers, k = k, r = r, it = it,
             eps = eps, init = init)
  numobs <- nrow(y)
  p <- ncol(y)
  k <- as.integer(k)
  r_latent <- as.integer(r)
  r_full <- c(p, r_latent)

  # Shared initialization: this is the key change.  No legacy
  # factanal -> princomp -> random-column fallback remains in DGMM.
  initial <- .rdmm_initialise(
    y = y,
    layers = as.integer(layers),
    k = k,
    r_full = r_full,
    init = init,
    init_est = init_est,
    psi_floor = psi_floor
  )

  H_list <- initial$H
  psi_list <- initial$psi
  psi_inv_list <- initial$psi.inv
  mu_list <- initial$mu
  w_list <- initial$w

  if (layers == 1L) {
    out <- deep.sem.alg.1(y, numobs, p, r_full[2L], k, H_list, psi_list,
                         psi_inv_list, mu_list, w_list, it, eps)
  } else if (layers == 2L) {
    out <- deep.sem.alg.2(y, numobs, p, r_full, k, H_list, psi_list,
                         psi_inv_list, mu_list, w_list, it, eps)
  } else if (layers == 3L) {
    out <- deep.sem.alg.3(y, numobs, p, r_full, k, H_list, psi_list,
                         psi_inv_list, mu_list, w_list, it, eps)
  } else {
    stop("layers must be 1, 2, or 3.")
  }

  out$lik <- out$likelihood
  output <- out[c("H", "w", "mu", "psi", "lik", "bic", "aic", "clc",
                  "icl_bic", "s", "h")]
  output <- c(output, list(
    k = k,
    r = r_latent,
    numobs = numobs,
    layers = as.integer(layers),
    seed = seed,
    init = init,
    init_est = init_est,
    psi_floor = psi_floor,
    initialization_engine = "shared_highdim_rdmm_v2"
  ))
  output$call <- match.call()
  class(output) <- "dgmm"

  invisible(output)
}
