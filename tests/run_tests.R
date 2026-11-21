################################################################################
# run_tests.R -- the rewritten code against the original 2019 code             #
#                                                                              #
#   Rscript tests/run_tests.R          (about a minute)                        #
#                                                                              #
# Runs the new functions (R/) and the 2019 functions (original_2019/) on the   #
# same random numbers and checks that they agree: exactly for the queue        #
# simulation, the bootstrap and the Oracle/SMC intervals, and up to rounding   #
# (1e-10 relative) for the likelihood-ratio estimators, which are computed in  #
# a different (faster) order.                                                  #
################################################################################

source(file.path("R", "utilities.R"))
source_library("R")
original <- new.env()
sys.source(file.path("original_2019", "MM1kFun.R"), envir = original)
sys.source(file.path("original_2019", "GSFun.R"), envir = original)

check <- function(name, new_value, old_value, tolerance = 0) {
  # Compare a value of the new code with the 2019 code and print the result
  # Inputs:
  # -- name: description of the check
  # -- new_value, old_value: values from the new and the 2019 code
  # -- tolerance: 0 = must be identical, otherwise relative tolerance
  # Output:
  # -- TRUE if the check passed

  if (tolerance == 0) {
    ok <- identical(unname(new_value), unname(old_value))
  } else {
    ok <- isTRUE(all.equal(unname(new_value), unname(old_value),
                           tolerance = tolerance))
  }
  cat(sprintf("%-8s %s\n", if (ok) "ok" else "FAILED", name))

  return(ok)
}

seed <- function(s) {
  # Fix the random numbers (R's default generator, as in 2019)

  set.seed(s, kind = "Mersenne-Twister", normal.kind = "Inversion")
}

#-- test settings
n_failed <- 0 # number of failed checks
k <- 10 # system capacity
r0 <- 50
r <- 200
alpha <- 0.1

# Queue simulation, at moderate and heavy load (many rejected arrivals)
for (setting in list(c(0.9, 1), c(1.6, 1), c(0.5, 1))) {
  for (s in 1:3) {
    seed(s)
    new <- queue_replication(setting[1], setting[2], k, r0, r)
    seed(s)
    old <- original$MM1kTIS_Lindley(setting[1], setting[2], k, r0, r)
    n_failed <- n_failed + !check(sprintf("queue replication, rates %s, seed %d",
                  paste(setting, collapse = "/"), s),
          c(new$output, new$n_interarrival, new$sum_interarrival,
            new$n_service, new$sum_service),
          c(old$Ybar, old$n_Arr, old$sum_Arr, old$n_Sr, old$sum_Sr))
  }
}
n_failed <- n_failed + !check("closed-form TIS", mm1k_true_tis(0.9, 1, k),
      original$MM1kTIS_Calc(0.9, 1, k))

# Bootstrap
seed(10)
new <- estimate_rates(0.9, 50, 20)
seed(10)
old <- original$MM1kParam_Bootstrap(0.9, 50, 20)
n_failed <- n_failed + !check("bootstrap", new, old)
lambda_boot <- new
seed(11)
mu_boot <- estimate_rates(1, 50, 20)

# Oracle and SMC intervals
n_failed <- n_failed + !check("Oracle CI", percentile_ci(oracle_estimates(lambda_boot, mu_boot, k), alpha),
      original$MM1kIUQ_Qtl_Calc(lambda_boot, mu_boot, k, alpha)$CI)
seed(12)
new <- percentile_ci(smc_estimates(lambda_boot, mu_boot, 3, k, r0, r), alpha)
seed(12)
old <- original$MM1kIUQ_Qtl_Nested(lambda_boot, mu_boot, k, r0, r, alpha, 3)$CI
n_failed <- n_failed + !check("SMC CI", new, old)

# ILR: R replications at one sampling parameter
n_rep <- 40
lambda_s <- ilr_sampling_rate(0.92, lambda_boot)
mu_s <- ilr_sampling_rate(1.03, mu_boot)
seed(13)
paths <- queue_replications(rep(lambda_s, n_rep), rep(mu_s, n_rep), k, r0, r)
new <- green_estimates(paths, rep(lambda_s, n_rep), rep(mu_s, n_rep),
                       lambda_boot, mu_boot)
seed(13)
old <- original$MM1kIUQ_Qtl_MIS(lambda_boot, mu_boot, lambda_s, mu_s,
                                k, r0, r, alpha, n_rep)
n_failed <- n_failed + !check("ILR estimates", new, old$TIS_B, tolerance = 1e-10)
n_failed <- n_failed + !check("ILR CI", percentile_ci(new, alpha), old$CI, tolerance = 1e-10)

# MLR: one replication at each of R sampling parameters
seed(14)
lambda_s <- mlr_sampling_rates(lambda_boot, n_rep)
mu_s <- mlr_sampling_rates(mu_boot, n_rep)
seed(14)
old_lambda_s <- original$MM1kParam_ResampThres(matrix(lambda_boot), n_rep)
old_mu_s <- original$MM1kParam_ResampThres(matrix(mu_boot), n_rep)
n_failed <- n_failed + !check("MLR sampling parameters", c(lambda_s, mu_s),
      c(old_lambda_s, old_mu_s))
seed(15)
paths <- queue_replications(lambda_s, mu_s, k, r0, r)
new <- green_estimates(paths, lambda_s, mu_s, lambda_boot, mu_boot)
seed(15)
old <- original$MM1kIUQ_Qtl_MIS(lambda_boot, mu_boot, lambda_s, mu_s,
                                k, r0, r, alpha, 1)
n_failed <- n_failed + !check("MLR estimates", new, old$TIS_B, tolerance = 1e-10)
n_failed <- n_failed + !check("MLR CI", percentile_ci(new, alpha), old$CI, tolerance = 1e-10)

cat(sprintf("\n%d failed\n", n_failed))
if (n_failed > 0) {
  quit(status = 1)
}
