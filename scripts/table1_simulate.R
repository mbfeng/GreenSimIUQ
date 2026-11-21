################################################################################
# table1_simulate.R                                                            #
# Simulation for Table 1 in Section 5 of the paper                             #
# M/M/1/10 queue, arrival rate 0.9, service rate 1, m = 50 real-world          #
# observations of each, replications of r = 100 customers (first r0 = 50 are   #
# warm-up), 90% CIs, simulation budget C = 1000, B = 1000 and B = 100          #
# Writes <out_dir>/table1_per_macro.csv.gz: one row per macro run, B and       #
# method, with the CI limits, coverage indicator, width and run times          #
# Usage (from the Code/ folder):                                               #
#   Rscript scripts/table1_simulate.R [n_macro] [n_cores] [out_dir] [r]        #
# -- n_macro: number of macro replications (default 10000, as in the paper)    #
# -- n_cores: number of parallel workers (default 1)                           #
# -- out_dir: output folder (default results/rerun)                            #
# -- r: customers per replication, including warm-up (default 100)             #
################################################################################

#-- command-line arguments
args <- commandArgs(trailingOnly = TRUE)
n_macro <- if (length(args) >= 1) as.integer(args[1]) else 10000
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1
out_dir <- if (length(args) >= 3) args[3] else file.path("results", "rerun")
r_customers <- if (length(args) >= 4) as.integer(args[4]) else 100
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#-- load the library and set up the parallel workers
source(file.path("R", "utilities.R"))
source_library("R")
setup_parallel(n_cores)

#-- experiment parameters (Section 5 of the paper)
LAMBDA <- 0.9 # true arrival rate
MU <- 1 # true service rate
K <- 10 # system capacity
N_OBS <- 50 # m, real-world observations of each input
R0 <- 50 # warm-up customers per replication
R_CUSTOMERS <- r_customers # customers per replication (incl. warm-up)
ALPHA <- 0.1 # 90% confidence intervals
BUDGET <- 1000 # C, total number of replications
B_VALUES <- c(1000, 100) # numbers of bootstrapped parameters (B)
R_GREEN <- c(1000, 100) # numbers of replications for ILR and MLR (R)
SEED <- 2019 # seed of the experiment

TRUE_TIS <- mm1k_true_tis(LAMBDA, MU, K) # for the coverage indicator

ci_row <- function(n_boot, method, n_rep, estimates, seconds_simulation,
                   seconds_lr) {
  # One row of results: the percentile CI from the estimates and its
  # coverage of the true TIS
  # Inputs:
  # -- n_boot, method, n_rep: B, method name, R
  # -- estimates: estimated TIS at the B bootstrapped parameters
  # -- seconds_simulation, seconds_lr: run times
  # Output:
  # -- data frame with one row

  ci <- percentile_ci(estimates, ALPHA)
  row <- data.frame(B = n_boot, method = method, R = n_rep,
                    lower = ci[1], upper = ci[2],
                    covered = as.numeric(ci[1] <= TRUE_TIS & TRUE_TIS <= ci[2]),
                    width = ci[2] - ci[1],
                    seconds_simulation = seconds_simulation,
                    seconds_lr = seconds_lr)

  return(row)
}

one_macro_run <- function(k) {
  # One macro replication: new real-world data, bootstrap, and all methods
  # Input:
  # -- k: index of the macro replication
  # Output:
  # -- data frame, one row per (B, method, R)
  # Design within one macro run:
  # -- the B = 100 panel uses the first 100 of the 1000 bootstrapped parameters
  # -- ILR and MLR with R = 100 use the first 100 of the 1000 replications
  #    (a 10%-budget version of the same experiment)

  #-- 1. real-world data -> estimated rates
  lambda_hat <- estimate_rates(LAMBDA, N_OBS, 1)
  mu_hat <- estimate_rates(MU, N_OBS, 1)

  #-- 2. parametric bootstrap
  lambda_boot_all <- estimate_rates(lambda_hat, N_OBS, max(B_VALUES))
  mu_boot_all <- estimate_rates(mu_hat, N_OBS, max(B_VALUES))

  rows <- list() # placeholder for the results
  for (n_boot in B_VALUES) {
    lambda_boot <- lambda_boot_all[1:n_boot]
    mu_boot <- mu_boot_all[1:n_boot]

    #-- Oracle
    start <- Sys.time()
    estimates <- oracle_estimates(lambda_boot, mu_boot, K)
    rows[[length(rows) + 1]] <- ci_row(n_boot, "Oracle", NA, estimates, NA,
                                       seconds_since(start))

    #-- SMC: R = C / B replications at every bootstrapped parameter
    n_inner <- BUDGET / n_boot
    start <- Sys.time()
    estimates <- smc_estimates(lambda_boot, mu_boot, n_inner, K, R0,
                               R_CUSTOMERS)
    rows[[length(rows) + 1]] <- ci_row(n_boot, "SMC", n_inner, estimates,
                                       seconds_since(start), NA)

    #-- ILR: max(R_GREEN) replications at one sampling parameter
    n_max <- max(R_GREEN)
    lambda_sampling <- rep(ilr_sampling_rate(lambda_hat, lambda_boot), n_max)
    mu_sampling <- rep(ilr_sampling_rate(mu_hat, mu_boot), n_max)
    start <- Sys.time()
    paths <- queue_replications(lambda_sampling, mu_sampling, K, R0,
                                R_CUSTOMERS)
    seconds_per_path <- seconds_since(start) / n_max
    for (n_rep in R_GREEN) {
      first <- 1:n_rep # the first n_rep replications
      start <- Sys.time()
      estimates <- green_estimates(paths[first, ], lambda_sampling[first],
                                   mu_sampling[first], lambda_boot, mu_boot)
      rows[[length(rows) + 1]] <- ci_row(n_boot, "ILR", n_rep, estimates,
                                         seconds_per_path * n_rep,
                                         seconds_since(start))
    }

    #-- MLR: one replication at each of max(R_GREEN) sampling parameters
    lambda_sampling <- mlr_sampling_rates(lambda_boot, n_max)
    mu_sampling <- mlr_sampling_rates(mu_boot, n_max)
    start <- Sys.time()
    paths <- queue_replications(lambda_sampling, mu_sampling, K, R0,
                                R_CUSTOMERS)
    seconds_per_path <- seconds_since(start) / n_max
    for (n_rep in R_GREEN) {
      first <- 1:n_rep # the first n_rep sampling parameters and replications
      start <- Sys.time()
      estimates <- green_estimates(paths[first, ], lambda_sampling[first],
                                   mu_sampling[first], lambda_boot, mu_boot)
      rows[[length(rows) + 1]] <- ci_row(n_boot, "MLR", n_rep, estimates,
                                         seconds_per_path * n_rep,
                                         seconds_since(start))
    }
  }
  result <- stack_rows(rows)
  result$macro <- k

  return(result)
}

#-- run the macro replications and save the per-macro results
start <- Sys.time()
results <- run_macros(n_macro, one_macro_run, seed = SEED)
message(sprintf("Table 1: %d macro runs in %.0f s", n_macro,
                seconds_since(start)))
write.csv(stack_rows(results),
          gzfile(file.path(out_dir, "table1_per_macro.csv.gz")),
          row.names = FALSE)
