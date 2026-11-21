# Efficient Input Uncertainty Quantification via Green Simulation Using Sample Path Likelihood Ratios — code

Code that reproduces Table 1 (Section 5) of

> Feng, B. M., & Song, E. (2019). Efficient input uncertainty quantification via green simulation
> using sample path likelihood ratios. In *Proceedings of the 2019 Winter Simulation Conference*,
> 3693–3704. IEEE. https://doi.org/10.1109/WSC40007.2019.9004835

If you use this code, please cite the paper.

## Quick start

```r
# R >= 4.3; from this folder:
install.packages(c("matrixStats", "future", "future.apply"))
source("reproduce_paper.R")   # or: Rscript reproduce_paper.R
```

By default (`RERUN_SIMULATIONS <- FALSE`) this computes Table 1 from the saved per-macro results of
the November 2026 rerun (`results/rerun/`, a few MB) in seconds and checks it against the paper.
Set `RERUN_SIMULATIONS <- TRUE` to run the 10,000 macro replications again (see [Run times](#run-times)).

## The experiment

An M/M/1/10 queue with arrival rate 0.9 and service rate 1. Both rates are unknown and estimated
from m = 50 "real-world" observations each. The quantity of interest is the steady-state expected
time in system (TIS), and the goal is a 90% confidence interval (CI) for it that accounts for the
input uncertainty. All methods are bootstrap percentile CIs: estimate the TIS at B bootstrapped
parameter pairs (B = 1,000 or 100) and take the 5% and 95% sample quantiles. They differ in how the
TIS at each bootstrapped parameter is estimated:

| Method (Table 1) | Estimate of the TIS at a bootstrapped parameter | Replications |
|---|---|---|
| Oracle | closed form (`mm1k_true_tis()`) | 0 |
| SMC | average of R replications at that parameter (Algorithm 1) | B·R = 1,000 |
| ILR | likelihood ratio estimator from R replications at ONE sampling parameter (Algorithm 2) | R = 1,000 or 100 |
| MLR | mixture likelihood ratio estimator from one replication at each of R sampling parameters | R = 1,000 or 100 |

Each replication serves 100 customers (Lindley's equation) and reports the average time in system
of customers 51–100, after a warm-up of 50 customers (the run length of the 2019 experiments;
`scripts/table1_simulate.R` takes another length as its fourth argument); the likelihood ratios use
all simulated interarrival and service times,
including those of the warm-up and of rejected arrivals. ILR samples at min{θ̂, 2·min_b θ_b} for
each rate; MLR samples R parameters with replacement from the bootstrapped rates below
2·min_b θ_b (the conditions of Proposition 1). Coverage is judged against the closed-form TIS at
the true rates. Everything is repeated over 10,000 macro replications, each with new "real-world"
data.

In one macro run the B = 100 panel uses the first 100 of the 1,000 bootstrapped parameters, and the
R = 100 cases of ILR and MLR use the first 100 of the 1,000 replications.

## Folder structure

```
reproduce_paper.R          entry point: Table 1 + verification
verify_against_paper.R     checks coverage and width of every Table 1 entry
R/
  mm1k_queue.R             the queue: closed-form TIS, one replication (Lindley's
                           equation), path log-densities, rate estimation
  iuq_methods.R            Oracle, SMC, ILR and MLR estimates; percentile CI
  tables.R                 Table 1 from the per-macro results
  utilities.R              loading the library, reproducible parallel macro runs, timing
scripts/table1_simulate.R  the simulation (writes results/<folder>/table1_per_macro.csv.gz)
results/rerun/             per-macro results of the November 2026 rerun
tests/run_tests.R          the rewritten functions against the 2019 code (original_2019/)
original_2019/             the code of June 2019 that implemented the methods (reference)
verify/                    the paper's printed values and the verification report
tables/                    output
```

## Code notation → paper notation

| Code | Paper |
|---|---|
| `lambda`, `mu` | arrival and service rates, θ |
| `lambda_hat`, `mu_hat` | θ̂, estimated from the m = 50 observations |
| `lambda_boot`, `mu_boot` | bootstrapped parameters θ_1, …, θ_B |
| `lambda_sampling`, `mu_sampling` | sampling parameters θ̃_s |
| `output` of `queue_replication()` | simulation output g(Z) |
| `n_boot`, `n_rep`, `n_inner` | B, R (green simulation), R (SMC) |
| `ALPHA = 0.1` | 1 − α = 0.9 |

## Reproducibility

Every macro run k has its own L'Ecuyer-CMRG random-number stream derived from one seed
(`run_macros()`, via `future.apply`), so the results do not depend on the number of cores. The
paper's 2019 runs used other random numbers, and their results were not kept, so a rerun agrees with
the paper within Monte Carlo error, not digit for digit.

## Verification

`verify_against_paper.R [results_dir]` compares the coverage and width of all 12 rows of Table 1
(24 numbers) with the values computed from a results folder:
z = (x − v)/√(se_x² + se_v² + h²/3), where se_v = se_x·√(K_x/10,000) is the Monte Carlo error of the
paper's value and h is half a unit in the last printed digit; |z| > 2.5 is flagged.

**Rerun (November 2026, `results/rerun/`, 10,000 macro replications, new random numbers):** all
24 entries agree with the paper within Monte Carlo error (largest |z| = 2.13; about 0.3 values above
2.5 would be expected by chance). `verify/verification_rerun.csv` has every entry.
`tests/run_tests.R` shows that the rewritten functions reproduce the 2019 code.

## Run times

10,000 macro replications took 27 minutes on 24 cores (Linux server, 2 × Xeon Gold 6126). Average
time per macro replication:

| B | Method | R | Simulation (s) | Likelihood ratios (s) |
|---|---|---|---|---|
| 1,000 | SMC | 1 | 1.07 | — |
| 1,000 | MLR | 1,000 / 100 | 0.46 / 0.05 | 0.079 / 0.004 |
| 1,000 | ILR | 1,000 / 100 | 0.46 / 0.05 | 0.038 / 0.004 |
| 100 | SMC | 10 | 0.52 | — |
| 100 | MLR | 1,000 / 100 | 0.46 / 0.05 | 0.046 / 0.001 |
| 100 | ILR | 1,000 / 100 | 0.46 / 0.05 | 0.004 / 0.001 |

The likelihood-ratio times are far below the run-time columns of Table 1 because the mixture
denominator is computed once per replication (see Change log); the paper's run times describe the
2019 implementation.

## Change log

* **2019-06** Code for the paper (`original_2019/`): the queue simulator and the four methods.
  The script with the exact Table 1 settings was not kept; the settings are taken from Section 5
  of the paper and from the run calls in `original_2019/runmacroreps.R` (replications of 100
  customers with 50 warm-up customers), and the rerun reproduces Table 1 with them.
* **2026-11** Maintained version (this folder):
  * `scripts/table1_simulate.R` runs Table 1 with the settings of Section 5; reproducible per-macro
    random-number streams, parallel macro runs, per-macro results instead of a saved workspace;
  * readable rewrite (descriptive names, plain loops, comments in the paper's notation);
    `tests/run_tests.R` shows that the queue simulation, the bootstrap and the Oracle/SMC intervals
    are identical to the 2019 code on the same random numbers, and the ILR/MLR estimates agree to
    rounding;
  * faster likelihood ratios: the mixture density in the denominator does not depend on the target
    parameter, so it is computed once per replication (about R² + B·R density evaluations instead
    of B·R², the same estimator); simulated interarrival times are stored in a preallocated vector;
  * `sample()` on a single allowed rate would have drawn from 1:x; the MLR sampling now picks by
    position (no effect on the results).

## License

MIT (see [LICENSE](LICENSE)).
