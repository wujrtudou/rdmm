# RDMM Simulation

This repository contains the R implementation and reproducibility code for a **Robust Deep Mixture Model (RDMM)** and its comparison with the **Deep Gaussian Mixture Model (DGMM)**. RDMM extends the Gaussian deep mixture formulation with pathway-wise Student-*t* tails, allowing the model to adapt to heavy-tailed observations and contamination while retaining the layered latent-variable structure of DGMM.

The repository includes the main simulation studies, reviewer-requested sensitivity analyses, publication tables and figures, and a rerun of the Chowdary gene-expression case study.

## Repository structure

```text
.
├── R/
│   ├── robustdeepgmm.R          # RDMM implementation
│   ├── dgmm/                    # DGMM implementation and supporting routines
│   ├── load_models.R            # Loads DGMM/RDMM code and dependencies
│   ├── simulation_parameters.R  # Data-generating models
│   ├── simulation_helpers.R     # Fitting, metrics, parallelization, recovery
│   └── output_helpers.R         # Aggregation, tables, figures
├── scripts/
│   ├── 00_install_packages.R
│   ├── 01_sim1_parameter_recovery.R
│   ├── 02_dimension_misspecification.R
│   ├── 03_gaussian_contamination.R
│   ├── 04_simulation_visuals_and_convergence.R
│   ├── 05_make_tables_and_figures.R
│   ├── 06_validate_outputs.R
│   ├── 07_covariance_sensitivity.R
│   ├── 08_computational_cost.R
│   ├── run_all.R
│   └── run_revision_experiments.R
├── results/
│   ├── raw/                     # Raw simulation/benchmark results
│   ├── tables/                  # CSV and LaTeX summary tables
│   └── figures/                 # PNG and PDF figures
├── chowdary_revision/           # Bundled Chowdary case-study outputs
├── case_chowdary.R              # Chowdary case-study rerun
├── simulation_multiworker.R     # Example multi-core Windows run script
└── RUN_REVISION_ANALYSES.md     # Notes for the added revision analyses
```

## Requirements

Use a recent version of **R**. The core simulations require:

- `mvtnorm`
- `corpcor`
- `ggplot2`

Install the core dependencies from the project root with:

```bash
Rscript scripts/00_install_packages.R
```

The Chowdary case study additionally requires:

- `ICGE` (for the `chowdary` dataset)
- `patchwork`

If you use `init = "mclust"` when calling the model directly, install `mclust` as well.

## Quick start

Run commands from the repository root. The scripts use the `PROJECT_ROOT` environment variable and otherwise default to the current working directory.

To load the model implementations in an interactive R session:

```r
source("R/load_models.R")
load_rdmm_models(".")
```

This makes `deepgmm()` and `robustdeepgmm()` available.

A minimal RDMM fit to a numeric matrix `X` is:

```r
source("R/load_models.R")
load_rdmm_models(".")

fit <- robustdeepgmm(
  y = X,
  layers = 2,
  k = c(2, 2),
  r = c(5, 2),
  it = 250,
  eps = 1e-3,
  init = "kmeans",
  init_est = "factanal",
  seed = 1,
  scale = TRUE,
  nu = 10,
  nu_structure = "common",
  estimate_nu = TRUE,
  method = "sem"
)
```

The latent dimensions must satisfy `p > r[1] > ... > r[layers] >= 1`, where `p` is the number of observed variables.

## Reproducing the simulation studies

### Full convenience run

The convenience runner executes the main simulations, diagnostics, covariance-sensitivity analysis, computational-cost benchmark, and then rebuilds tables and figures:

```bash
Rscript scripts/run_all.R
```

The default Monte Carlo settings are publication scale (`N_REP = 500` for the main experiments and `COST_N_REP = 20` for the timing benchmark), so a full run can be computationally expensive.

After a run, validate the generated raw outputs with:

```bash
Rscript scripts/06_validate_outputs.R
```

### Individual experiments

| Script | Analysis |
|---|---|
| `01_sim1_parameter_recovery.R` | Student-*t* tail-weight experiment and RDMM parameter recovery |
| `02_dimension_misspecification.R` | Sensitivity to under/over-specifying latent dimensions |
| `03_gaussian_contamination.R` | Gaussian DGMM data with nested multivariate-*t* contamination |
| `04_simulation_visuals_and_convergence.R` | Simulated-data visualization and multiple-start convergence traces |
| `07_covariance_sensitivity.R` | Sensitivity to correlated AR(1) residual covariance when fitting diagonal residual covariance |
| `08_computational_cost.R` | Sequential runtime/scalability benchmark over sample size and dimension |
| `05_make_tables_and_figures.R` | Aggregates raw shards and regenerates publication tables/figures |
| `06_validate_outputs.R` | Checks row counts and reports fit-failure rates |

For example, to run only Simulation 1:

```r
Sys.setenv(PROJECT_ROOT = getwd(), N_REP = 500, CORES = 4)
source("scripts/01_sim1_parameter_recovery.R")
```

## Parallel runs and sharding

The Monte Carlo experiments support both multi-core execution and explicit replication shards. Useful environment variables are:

- `PROJECT_ROOT`: repository root; defaults to `.`
- `N_REP`: total/default number of Monte Carlo replications; default `500`
- `REP_START`, `REP_END`: inclusive replication range for a shard
- `CORES`: number of workers; default `1`
- `COND`: run a single indexed design condition where supported
- `MAX_ITER`: maximum fitting iterations; default `300`

Example: split a 500-replication experiment into two jobs:

```r
# Job 1
Sys.setenv(PROJECT_ROOT = getwd(), N_REP = 500,
           REP_START = 1, REP_END = 250, CORES = 8)
source("scripts/01_sim1_parameter_recovery.R")

# Job 2
Sys.setenv(PROJECT_ROOT = getwd(), N_REP = 500,
           REP_START = 251, REP_END = 500, CORES = 8)
source("scripts/01_sim1_parameter_recovery.R")
```

Each job writes a separate CSV shard under `results/raw/`. `scripts/05_make_tables_and_figures.R` reads all matching shards and removes duplicate condition/replication/model rows before summarizing them.

Parallel execution uses a PSOCK cluster on Windows and `mclapply()` on Unix-like systems.

## Revision analyses

The two added reviewer-driven analyses can be run together with:

```bash
Rscript scripts/run_revision_experiments.R
```

See [`RUN_REVISION_ANALYSES.md`](RUN_REVISION_ANALYSES.md) for additional details and quick-test commands.

### Residual-covariance sensitivity

`scripts/07_covariance_sensitivity.R` generates data with AR(1) residual correlation at both layers using

```text
rho = 0.0, 0.2, 0.4, 0.6
```

while fitting both RDMM and DGMM under the diagonal-residual-covariance specification. The marginal residual variances are kept fixed.

Main outputs include:

```text
results/tables/Table_covariance_sensitivity.csv
results/tables/Table_covariance_sensitivity_recovery.csv
results/figures/Figure_covariance_sensitivity_clustering.png
```

### Computational cost

`scripts/08_computational_cost.R` holds `K = (2, 2)` and `r = (5, 2)` fixed and benchmarks:

- `n = 250, 500, 1000, 2000` at `p = 20`
- `p = 20, 40, 80` at `n = 1000`

Timing is intentionally **sequential**. Do not parallelize the benchmark if wall-clock comparisons are intended for reporting, because CPU contention would make the timing results difficult to interpret.

The benchmark accepts separate `COST_*` controls, including `COST_N_REP`, `COST_REP_START`, `COST_REP_END`, `COST_N_GRID`, `COST_P_GRID`, `COST_COND`, `COST_MAX_ITER`, `COST_MIN_ITER`, and `COST_MOVING_WINDOW`.

## Generating tables and figures

Once the required raw CSV files exist, regenerate publication outputs with:

```bash
Rscript scripts/05_make_tables_and_figures.R
```

The script writes summary tables in both CSV and LaTeX form to `results/tables/` and figures in PNG and PDF form to `results/figures/`.

Because the aggregation code reads all matching raw shards, you can run simulation jobs separately and build the final summaries afterward.

## Chowdary case study

`case_chowdary.R` reruns the Chowdary gene-expression analysis across the manuscript architecture grid for both DGMM and RDMM, using 10 random starts per candidate architecture and selecting valid fits by final observed-data likelihood within an architecture and BIC across architectures. It also creates selected-model tables and latent-space/PCA visualizations.

Before running it, install the additional dependencies:

```r
install.packages(c("ICGE", "patchwork"))
```

Then **edit `project_root` near the top of `case_chowdary.R`**. The current script contains a machine-specific Windows path:

```r
project_root <- "C:/Users/wujrt/Desktop/rdmm-main"
```

Replace it with the path to your local checkout. The current script writes newly generated case-study files to:

```text
results/chowdary_revision/
```

The archive also contains previously generated Chowdary outputs under the top-level `chowdary_revision/` directory.

Run the case study with:

```bash
Rscript case_chowdary.R
```

## Output conventions

Raw experiment files contain fit-level results such as:

- estimated degrees of freedom for RDMM;
- first-layer and pathway adjusted Rand index (ARI);
- first-layer and pathway misclassification rate (MR);
- log-likelihood and BIC;
- iteration count and elapsed fit time;
- fit success/failure status and error messages.

Parameter-recovery experiments additionally report identifiable or rotation-invariant recovery summaries such as mixing-proportion RMSE, observed-location RMSE, residual-variance RMSE, loading-subspace error, and pathway-level location/scale recovery.

Random seeds are generated deterministically from the condition and replication number so that simulation shards can be reproduced independently.

## Notes

- Run scripts from the project root unless you explicitly set `PROJECT_ROOT`.
- Publication-scale runs are computationally intensive; use replication shards for the main Monte Carlo studies when appropriate.
- Keep the computational-cost benchmark sequential for interpretable timings.
- `simulation_multiworker.R` is an example configured with a machine-specific Windows path. Edit both `setwd()` and `PROJECT_ROOT` before using it on another machine.
- Existing files in `results/` provide examples of the expected output layout and naming conventions.

