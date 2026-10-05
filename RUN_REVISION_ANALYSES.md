# Reviewer revision analyses

This project now contains two additional analyses aimed directly at Reviewer 1,
Comments 5 and 6.

## 1. Residual-covariance sensitivity

Script: `scripts/07_covariance_sensitivity.R`

Data are generated with AR(1) residual covariance at both layers, with
`rho = 0, 0.2, 0.4, 0.6`, while preserving the manuscript marginal residual
variances.  RDMM and DGMM are still fitted with the manuscript's diagonal-Psi
specification.  This measures sensitivity to violation of the diagonal
residual-covariance assumption without changing the fitted model.

Publication-scale default: 500 replications per rho value.  As with the other
Monte Carlo experiments, this script can be sharded using `REP_START`,
`REP_END`, `COND`, and `CORES`.

Example quick test:

```r
Sys.setenv(PROJECT_ROOT = getwd(), N_REP = 5, REP_END = 5, CORES = 1)
source("scripts/07_covariance_sensitivity.R")
```

Main outputs after running `scripts/05_make_tables_and_figures.R`:

- `results/tables/Table_covariance_sensitivity.csv/.tex`
- `results/tables/Table_covariance_sensitivity_recovery.csv/.tex`
- `results/figures/Figure_covariance_sensitivity_clustering.png/.pdf`

## 2. Computational cost / scalability

Script: `scripts/08_computational_cost.R`

The benchmark holds `K=(2,2)` and `r=(5,2)` fixed and varies:

- `n = 250, 500, 1000, 2000` at `p=20`;
- `p = 20, 40, 80` at `n=1000`.

The shared baseline `(n,p)=(1000,20)` is run only once.  Each condition is timed
for DGMM and RDMM separately.  The output records total wall-clock fit time,
iteration count, seconds per iteration, clustering metrics, and success status.
A system-information text file is also saved so the hardware/software context
can be stated in the paper.

**Important:** runtime benchmarking is intentionally sequential.  Do not run
replications in parallel when collecting publication timing numbers because CPU
contention makes per-fit wall-clock comparisons difficult to interpret.

Default: 20 repetitions per condition.  A quick smoke test is:

```r
Sys.setenv(PROJECT_ROOT = getwd(), COST_N_REP = 2, COST_REP_END = 2,
           COST_N_GRID = "250,1000", COST_P_GRID = "20,40")
source("scripts/08_computational_cost.R")
source("scripts/05_make_tables_and_figures.R")
```

For the full benchmark:

```r
Sys.setenv(PROJECT_ROOT = getwd(), COST_N_REP = 20)
source("scripts/08_computational_cost.R")
source("scripts/05_make_tables_and_figures.R")
```

Optional environment variables:

- `COST_N_REP`, `COST_REP_START`, `COST_REP_END`
- `COST_N_GRID` (comma-separated integer vector)
- `COST_P_GRID` (comma-separated integer vector, default values are >=20)
- `COST_BASE_N`, `COST_BASE_P`
- `COST_COND` for one benchmark condition only
- `COST_MAX_ITER`, `COST_MIN_ITER`, `COST_MOVING_WINDOW`

Main outputs:

- `results/raw/computational_cost_*.csv`
- `results/raw/computational_cost_system_info.txt`
- `results/tables/Table_computational_cost.csv/.tex`
- `results/figures/Figure_computational_cost.png/.pdf`

## Recommended workflow

Run the large Monte Carlo simulations and covariance sensitivity as parallel or
sharded jobs if desired.  Run the computational-cost benchmark as a separate,
quiet, sequential job on the machine whose specifications you will report.
Finally run `scripts/05_make_tables_and_figures.R` and
`scripts/06_validate_outputs.R`.
