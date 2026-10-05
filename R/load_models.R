# Load the original DGMM and RDMM implementations used by the manuscript.
# Run from the project root or pass project_root explicitly.

load_rdmm_models <- function(project_root = ".") {
  required <- c("mvtnorm", "corpcor")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Run scripts/00_install_packages.R first.")
  }

  assign("rmvnorm", mvtnorm::rmvnorm, envir = .GlobalEnv)
  assign("dmvnorm", mvtnorm::dmvnorm, envir = .GlobalEnv)
  assign("is.positive.definite", corpcor::is.positive.definite, envir = .GlobalEnv)
  assign("make.positive.definite", corpcor::make.positive.definite, envir = .GlobalEnv)

  core <- c(
    "misc.R", "ginv.R", "adjustedRandIndex.R", "initial_clustering.R",
    "fix_para.R", "factanal_para.R", "ppca_para.R", "valid_args.R",
    "compute_est.R", "compute_lik.R", "deep.sem.alg.1.R",
    "deep.sem.alg.2.R", "deep.sem.alg.3.R"
  )
  for (f in core) source(file.path(project_root, "R", "dgmm", f), local = .GlobalEnv)

  # robustdeepgmm.R defines the shared high-dimensional initializer used by both methods.
  source(file.path(project_root, "R", "robustdeepgmm.R"), local = .GlobalEnv)
  source(file.path(project_root, "R", "dgmm", "deepgmm.R"), local = .GlobalEnv)
  source(file.path(project_root, "R", "dgmm", "print.R"), local = .GlobalEnv)

  invisible(TRUE)
}
