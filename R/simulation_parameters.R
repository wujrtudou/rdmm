# Simulation parameters matching the manuscript.

make_manuscript_parameters <- function() {
  p <- 20L
  r1 <- 5L
  r2 <- 2L
  k1 <- 2L
  k2 <- 2L

  eta1 <- matrix(0, nrow = p, ncol = k1)
  eta1[1:3, 1L] <- -2.5
  eta1[1:3, 2L] <-  2.5

  eta2 <- matrix(0, nrow = r1, ncol = k2)
  eta2[1L, 1L] <-  1.25
  eta2[1L, 2L] <- -1.25

  a <- matrix(c(1.0, 0.9, 1.1, 0.8), ncol = 1L)
  A <- kronecker(diag(r1), a)
  D2 <- diag(c(1.10, 0.90, 1.05, 0.95, 1.00))
  lambda1 <- list(A, A %*% D2)

  lambda2 <- list(
    matrix(c(
      1.0,  0.0,
      0.8,  0.2,
      0.0,  1.0,
      0.2,  0.8,
      0.6, -0.6
    ), nrow = r1, byrow = TRUE),
    matrix(c(
       0.9,  0.1,
       0.7, -0.2,
       0.1,  0.9,
      -0.2,  0.7,
       0.5,  0.5
    ), nrow = r1, byrow = TRUE)
  )

  psi1 <- list(0.5 * diag(p), 0.5 * diag(p))
  psi2 <- list(0.3 * diag(r1), 0.3 * diag(r1))

  pi1 <- c(0.5, 0.5)
  pi2 <- c(0.5, 0.5)

  paths <- as.matrix(expand.grid(
    layer1 = seq_len(k1),
    layer2 = seq_len(k2),
    KEEP.OUT.ATTRS = FALSE
  ))

  list(
    p = p,
    r = c(r1, r2),
    k = c(k1, k2),
    pi = list(pi1, pi2),
    eta = list(eta1, eta2),
    lambda = list(lambda1, lambda2),
    psi = list(psi1, psi2),
    paths = paths
  )
}

encode_path_id <- function(path_matrix, k = c(2L, 2L)) {
  path_matrix <- as.matrix(path_matrix)
  mult <- cumprod(c(1L, head(as.integer(k), -1L)))
  as.integer(1L + rowSums((path_matrix - 1L) * matrix(
    mult, nrow(path_matrix), ncol(path_matrix), byrow = TRUE
  )))
}

rmvnorm_chol <- function(n, mean, sigma) {
  sigma <- (sigma + t(sigma)) / 2
  z <- matrix(rnorm(n * length(mean)), nrow = n)
  sweep(z %*% chol(sigma), 2L, mean, "+")
}

simulate_student_rdmm <- function(n = 1000L, nu = 4,
                                  pars = make_manuscript_parameters(),
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  n <- as.integer(n)
  p <- pars$p
  r1 <- pars$r[1L]
  r2 <- pars$r[2L]

  s1 <- sample.int(pars$k[1L], n, replace = TRUE, prob = pars$pi[[1L]])
  s2 <- sample.int(pars$k[2L], n, replace = TRUE, prob = pars$pi[[2L]])
  pathway <- cbind(s1, s2)

  # Shared precision W ~ Gamma(nu/2, nu/2); conditional variances are divided by W.
  W <- rgamma(n, shape = nu / 2, rate = nu / 2)
  inv_sqrt_W <- 1 / sqrt(W)

  z2 <- matrix(rnorm(n * r2), nrow = n, ncol = r2)
  z2 <- z2 * inv_sqrt_W

  z1 <- matrix(NA_real_, n, r1)
  for (b in seq_len(pars$k[2L])) {
    idx <- which(s2 == b)
    if (!length(idx)) next
    Lam <- pars$lambda[[2L]][[b]]
    mu <- pars$eta[[2L]][, b]
    mean_part <- sweep(z2[idx, , drop = FALSE] %*% t(Lam), 2L, mu, "+")
    eps <- rmvnorm_chol(length(idx), rep(0, r1), pars$psi[[2L]][[b]])
    eps <- eps * inv_sqrt_W[idx]
    z1[idx, ] <- mean_part + eps
  }

  y <- matrix(NA_real_, n, p)
  for (a in seq_len(pars$k[1L])) {
    idx <- which(s1 == a)
    if (!length(idx)) next
    Lam <- pars$lambda[[1L]][[a]]
    mu <- pars$eta[[1L]][, a]
    mean_part <- sweep(z1[idx, , drop = FALSE] %*% t(Lam), 2L, mu, "+")
    eps <- rmvnorm_chol(length(idx), rep(0, p), pars$psi[[1L]][[a]])
    eps <- eps * inv_sqrt_W[idx]
    y[idx, ] <- mean_part + eps
  }

  colnames(y) <- paste0("Y", seq_len(p))
  path_id <- encode_path_id(pathway, pars$k)

  list(
    y = y, z1 = z1, z2 = z2, W = W,
    layer1 = s1, layer2 = s2,
    pathway = pathway, path_id = path_id,
    nu = nu, parameters = pars
  )
}

simulate_gaussian_dgmm <- function(n = 1000L,
                                   pars = make_manuscript_parameters(),
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  n <- as.integer(n)
  p <- pars$p
  r1 <- pars$r[1L]
  r2 <- pars$r[2L]

  s1 <- sample.int(pars$k[1L], n, replace = TRUE, prob = pars$pi[[1L]])
  s2 <- sample.int(pars$k[2L], n, replace = TRUE, prob = pars$pi[[2L]])
  pathway <- cbind(s1, s2)

  z2 <- matrix(rnorm(n * r2), nrow = n, ncol = r2)

  z1 <- matrix(NA_real_, n, r1)
  for (b in seq_len(pars$k[2L])) {
    idx <- which(s2 == b)
    if (!length(idx)) next
    Lam <- pars$lambda[[2L]][[b]]
    mu <- pars$eta[[2L]][, b]
    mean_part <- sweep(z2[idx, , drop = FALSE] %*% t(Lam), 2L, mu, "+")
    eps <- rmvnorm_chol(length(idx), rep(0, r1), pars$psi[[2L]][[b]])
    z1[idx, ] <- mean_part + eps
  }

  y <- matrix(NA_real_, n, p)
  for (a in seq_len(pars$k[1L])) {
    idx <- which(s1 == a)
    if (!length(idx)) next
    Lam <- pars$lambda[[1L]][[a]]
    mu <- pars$eta[[1L]][, a]
    mean_part <- sweep(z1[idx, , drop = FALSE] %*% t(Lam), 2L, mu, "+")
    eps <- rmvnorm_chol(length(idx), rep(0, p), pars$psi[[1L]][[a]])
    y[idx, ] <- mean_part + eps
  }

  colnames(y) <- paste0("Y", seq_len(p))
  path_id <- encode_path_id(pathway, pars$k)

  list(
    y = y, z1 = z1, z2 = z2,
    layer1 = s1, layer2 = s2,
    pathway = pathway, path_id = path_id,
    parameters = pars
  )
}

make_nested_t_contamination <- function(y, eps_grid = c(0, .1, .2, .3, .4, .5),
                                        df = 4, scale_multiplier = 2,
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)

  ord <- sample.int(n)
  # One chi-square draw per observation gives a multivariate t perturbation.
  denom <- sqrt(rchisq(n, df = df) / df)
  G <- matrix(rnorm(n * p, sd = sqrt(scale_multiplier)), nrow = n, ncol = p)
  delta <- G / denom

  out <- vector("list", length(eps_grid))
  names(out) <- sprintf("%.2f", eps_grid)
  for (i in seq_along(eps_grid)) {
    eps <- eps_grid[i]
    m <- floor(eps * n)
    yc <- y
    if (m > 0L) {
      idx <- ord[seq_len(m)]
      yc[idx, ] <- yc[idx, , drop = FALSE] + delta[idx, , drop = FALSE]
    } else {
      idx <- integer(0)
    }
    out[[i]] <- list(y = yc, contaminated = idx)
  }

  list(levels = out, order = ord, perturbation = delta, eps_grid = eps_grid)
}

# Create a controlled violation of the diagonal residual-covariance assumption.
# Marginal residual variances are kept identical to the manuscript setting,
# while AR(1) residual correlation is introduced within each layer.
make_correlated_residual_parameters <- function(rho1 = 0, rho2 = rho1,
                                                pars = make_manuscript_parameters()) {
  if (length(rho1) != 1L || length(rho2) != 1L ||
      !is.finite(rho1) || !is.finite(rho2) ||
      abs(rho1) >= 1 || abs(rho2) >= 1) {
    stop("rho1 and rho2 must be finite scalars strictly between -1 and 1.")
  }

  ar1_cor <- function(d, rho) {
    idx <- seq_len(d)
    rho ^ abs(outer(idx, idx, "-"))
  }

  out <- pars
  for (a in seq_len(out$k[1L])) {
    base <- out$psi[[1L]][[a]]
    sdv <- sqrt(pmax(diag(base), 0))
    D <- diag(sdv, length(sdv))
    out$psi[[1L]][[a]] <- D %*% ar1_cor(length(sdv), rho1) %*% D
  }
  for (b in seq_len(out$k[2L])) {
    base <- out$psi[[2L]][[b]]
    sdv <- sqrt(pmax(diag(base), 0))
    D <- diag(sdv, length(sdv))
    out$psi[[2L]][[b]] <- D %*% ar1_cor(length(sdv), rho2) %*% D
  }
  out
}

# Parameters for the computational-cost benchmark when p is varied.
# The two-layer architecture and latent dimensions are held fixed at the
# manuscript values (K = (2,2), r = (5,2)).  For p > 20, the manuscript's
# 20-variable first-layer loading/location pattern is repeated so that changing
# p changes the amount of observed information without changing the fitted
# architecture.  The benchmark uses p >= 20 by default.
make_runtime_parameters <- function(p = 20L,
                                    pars = make_manuscript_parameters()) {
  p <- as.integer(p)
  if (length(p) != 1L || !is.finite(p) || p < pars$p) {
    stop("For the runtime benchmark, p must be an integer >= ", pars$p, ".")
  }

  if (p == pars$p) return(pars)

  out <- pars
  idx <- rep(seq_len(pars$p), length.out = p)
  out$p <- p

  # Repeat the manuscript signal pattern in blocks of 20 variables.
  out$eta[[1L]] <- pars$eta[[1L]][idx, , drop = FALSE]
  out$lambda[[1L]] <- lapply(pars$lambda[[1L]], function(L) {
    L[idx, , drop = FALSE]
  })

  # Preserve the manuscript marginal residual variance (0.5) for every
  # observed variable.  The runtime study is not a covariance-misspecification
  # experiment, so Psi remains diagonal here.
  out$psi[[1L]] <- lapply(seq_len(out$k[1L]), function(a) {
    v <- mean(diag(pars$psi[[1L]][[a]]))
    v * diag(p)
  })

  out
}

