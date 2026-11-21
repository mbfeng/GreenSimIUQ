################################################################################
# iuq_methods.R                                                                #
# Input-uncertainty confidence intervals (CIs) compared in Section 5 and       #
# Table 1 of the paper: bootstrap percentile CIs from estimates of the         #
# conditional mean mu(theta_b) at B bootstrapped parameters                    #
# -- Oracle: Algorithm 1 without inner simulation, mu(theta_b) in closed form  #
# -- SMC: Algorithm 1 with inner simulation, R reps at every theta_b           #
# -- ILR: Algorithm 2, individual likelihood ratio, R reps at ONE sampling     #
#    parameter, reused for every theta_b                                       #
# -- MLR: Algorithm 2, mixture likelihood ratio, one rep at each of R sampling #
#    parameters, reused for every theta_b                                      #
# Paper notation -> code:                                                      #
# -- bootstrapped parameters theta_b -> lambda_boot, mu_boot                   #
# -- sampling parameters theta~_s -> lambda_sampling, mu_sampling              #
# -- significance level alpha (1 - alpha = 0.9) -> alpha                       #
################################################################################

# ===== Percentile confidence interval
percentile_ci <- function(estimates, alpha) {
  # Bootstrap percentile CI
  # Inputs:
  # -- estimates: estimated conditional means at the bootstrapped parameters
  # -- alpha: significance level, the CI has level 1 - alpha
  # Output:
  # -- vector of 2 values, the alpha/2 and 1-alpha/2 sample quantiles
  #    (quantile type 4, linear interpolation of the empirical cdf, as in the
  #    2019 code)

  ci <- quantile(estimates, c(alpha / 2, 1 - alpha / 2), type = 4,
                 names = FALSE)

  return(ci)
}

# ===== Algorithm 1: Oracle and SMC
oracle_estimates <- function(lambda_boot, mu_boot, k) {
  # Oracle: closed-form TIS at every bootstrapped parameter
  # Inputs:
  # -- lambda_boot, mu_boot: vectors of length B, bootstrapped rates
  # -- k: system capacity
  # Output:
  # -- vector of length B, TIS at the bootstrapped parameters

  return(mm1k_true_tis(lambda_boot, mu_boot, k))
}

smc_estimates <- function(lambda_boot, mu_boot, n_inner, k, r0, r) {
  # SMC: average of n_inner replications at every bootstrapped parameter
  # Inputs:
  # -- lambda_boot, mu_boot: vectors of length B, bootstrapped rates
  # -- n_inner: number of replications per bootstrapped parameter (R)
  # -- k, r0, r: queue capacity, warm-up and total customers
  # Output:
  # -- vector of length B, estimated TIS at the bootstrapped parameters

  n_boot <- length(lambda_boot) # number of bootstrapped parameters (B)
  estimates <- numeric(n_boot) # placeholder for the estimates
  for (b in 1:n_boot) {
    lambda_b <- rep(lambda_boot[b], n_inner) # same parameter for all reps
    mu_b <- rep(mu_boot[b], n_inner)
    paths <- queue_replications(lambda_b, mu_b, k, r0, r)
    estimates[b] <- mean(paths$output)
  }

  return(estimates)
}

# ===== Algorithm 2: sampling parameters
ilr_sampling_rate <- function(rate_hat, rate_boot) {
  # Sampling parameter of ILR, min{theta^, 2 min_b theta_b} (Section 5), which
  # makes the variance of the likelihood ratios finite (Proposition 1);
  # applied to the arrival and the service rate separately
  # Inputs:
  # -- rate_hat: rate estimated from the real-world data (theta^)
  # -- rate_boot: vector of bootstrapped rates
  # Output:
  # -- the sampling rate (as in the 2019 code, 1e-6 is subtracted from
  #    2 min_b theta_b to make the bound strict)

  return(min(rate_hat, 2 * min(rate_boot) - 1e-6))
}

mlr_sampling_rates <- function(rate_boot, n_sampling) {
  # Sampling parameters of MLR: random draws with replacement from the
  # bootstrapped rates below 2 min_b theta_b (Section 5)
  # Inputs:
  # -- rate_boot: vector of bootstrapped rates
  # -- n_sampling: number of sampling rates (R)
  # Output:
  # -- vector of n_sampling sampling rates

  allowed <- rate_boot[rate_boot < 2 * min(rate_boot) - 1e-6] # feasible rates
  # sample() would treat a single number x as 1:x, so pick by position
  picked <- sample.int(length(allowed), size = n_sampling, replace = TRUE)

  return(allowed[picked])
}

# ===== Algorithm 2: likelihood ratio estimators
green_estimates <- function(paths, lambda_sampling, mu_sampling,
                            lambda_boot, mu_boot) {
  # Likelihood ratio (green simulation) estimates of mu(theta_b) at every
  # bootstrapped parameter, reusing the replications at the sampling
  # parameters
  # Inputs:
  # -- paths: data frame of R replications (queue_replications()),
  #    replication i simulated at (lambda_sampling[i], mu_sampling[i])
  # -- lambda_sampling, mu_sampling: vectors of length R, sampling rates
  # -- lambda_boot, mu_boot: vectors of length B, bootstrapped rates
  # Output:
  # -- vector of length B, estimated TIS at the bootstrapped parameters
  # formula: estimate_b = (1/R) sum_i g(Z_i) f_b(Z_i) / ((1/R) sum_s f_s(Z_i))
  # -- MLR: the denominator is the mixture of the R sampling densities
  #    (balance heuristic)
  # -- ILR: all R sampling parameters are equal, the mixture is that density
  # Efficiency: the denominator does not depend on b, so it is computed once
  # per replication (on the log scale, with the log-sum-exp trick against
  # overflow): about R^2 + B R density evaluations instead of B R^2 in the
  # 2019 code; same estimator up to rounding

  n_paths <- nrow(paths) # number of replications (R)

  #-- log mixture density (1/R) sum_s f_s(Z_i) of every replication i
  is_ilr <- length(unique(lambda_sampling)) == 1 &&
    length(unique(mu_sampling)) == 1 # one sampling parameter
  if (is_ilr) {
    log_mixture <- as.vector(log_path_density(paths, lambda_sampling[1],
                                              mu_sampling[1]))
  } else {
    log_sampling <- log_path_density(paths, lambda_sampling, mu_sampling) # [i, s]
    largest <- matrixStats::rowMaxs(log_sampling) # for numerical stability
    log_mixture <- largest + log(rowMeans(exp(log_sampling - largest)))
  }

  #-- likelihood ratios and estimates for all bootstrapped parameters
  log_target <- log_path_density(paths, lambda_boot, mu_boot) # [i, b]
  likelihood_ratio <- exp(log_target - log_mixture) # [i, b] = f_b / mixture
  estimates <- as.vector(paths$output %*% likelihood_ratio) / n_paths

  return(estimates)
}
