################################################################################
# mm1k_queue.R                                                                 #
# The M/M/1/k queue of Section 5 in the paper                                  #
# -- exponential interarrival times with rate lambda                           #
# -- exponential service times with rate mu                                    #
# -- one server, room for k customers (an arrival that finds k customers is    #
#    rejected)                                                                 #
# Performance measure: steady-state expected time in system (TIS)              #
# Paper notation -> code:                                                      #
# -- input parameters theta = (arrival rate, service rate) -> (lambda, mu)     #
# -- real-world observations, m = 50 of each -> n_obs                          #
# -- one simulation replication, sample path Z -> queue_replication()          #
# -- simulation output g(Z): average TIS of the customers after warm-up        #
################################################################################

# ===== Closed-form solution
mm1k_true_tis <- function(lambda, mu, k) {
  # Steady-state expected time in system of an M/M/1/k queue (the "oracle"
  # in the paper), vectorized over lambda and mu
  # Inputs:
  # -- lambda: interarrival rate(s)
  # -- mu: service rate(s)
  # -- k: system capacity
  # Output:
  # -- long-run mean time in system, one value per (lambda, mu)

  rho <- lambda / mu # traffic intensity
  number_in_system <- rho / (1 - rho) -
    (k + 1) * rho^(k + 1) / (1 - rho^(k + 1)) # mean number in system
  prob_full <- rho^k / ((1 - rho^(k + 1)) / (1 - rho)) # P(k customers in system)
  tis <- number_in_system / (lambda * (1 - prob_full)) # Little's law

  return(tis)
}

# ===== Simulation
queue_replication <- function(lambda, mu, k, r0, r) {
  # One simulation replication of the M/M/1/k queue, started empty, using
  # Lindley's equation W_i = max(0, W_{i-1} + S_{i-1} - A_i)
  # Inputs:
  # -- lambda: interarrival rate
  # -- mu: service rate
  # -- k: system capacity
  # -- r0: number of warm-up customers
  # -- r: total number of customers served (including warm-up)
  # Output: a list of 5 items
  # -- output: average time in system of customers r0+1, ..., r, i.e., g(Z)
  # -- n_interarrival, sum_interarrival: number and sum of ALL simulated
  #    interarrival times (including warm-up and rejected arrivals)
  # -- n_service, sum_service: number and sum of all simulated service times
  # Note: for exponential inputs these four numbers are sufficient statistics
  # of the sample path's likelihood (see log_path_density()). The random
  # numbers are drawn in the same order as in the 2019 code
  # (original_2019/MM1kFun.R, MM1kTIS_Lindley): first all r service times,
  # then the interarrival times one by one.

  service <- rexp(r, mu) # service times of the r customers
  wait <- numeric(r) # placeholder, waiting times
  time_in_system <- numeric(r) # placeholder, times in system
  departure <- numeric(r) # placeholder, departure times
  interarrival <- numeric(2 * r) # placeholder, all simulated interarrival times
  n_interarrival <- 0 # number of simulated interarrival times so far

  #-- customer 1 arrives to an empty system
  gap <- rexp(1, lambda) # first interarrival time
  n_interarrival <- n_interarrival + 1
  interarrival[n_interarrival] <- gap
  clock <- gap # time of the latest accepted arrival
  time_in_system[1] <- service[1] # no waiting for the first customer
  departure[1] <- clock + service[1]
  n_in_system <- 1 # customers in system at the latest arrival

  #-- customers 2, ..., r
  for (i in 2:r) {
    since_last <- 0 # time since the latest accepted arrival
    repeat { # until an arrival is accepted
      gap <- rexp(1, lambda) # time to the next arrival
      n_interarrival <- n_interarrival + 1
      if (n_interarrival > length(interarrival)) {
        interarrival <- c(interarrival, numeric(r)) # extend the placeholder
      }
      interarrival[n_interarrival] <- gap
      arrival_gap <- since_last + gap # time since the latest accepted arrival
      if (clock + arrival_gap > departure[i - 1]) {
        n_ahead <- 0 # the arrival finds the system empty
        break
      }
      # n_ahead: number of the n_in_system customers still in the system
      for (n_ahead in n_in_system:1) {
        if (clock + arrival_gap < departure[i - n_ahead]) {
          break
        }
      }
      if (n_ahead < k) {
        break # arrival accepted
      }
      since_last <- arrival_gap # arrival rejected: the system is full
    }
    n_in_system <- n_ahead + 1 # including the new customer
    clock <- clock + arrival_gap # advance the clock to the new arrival
    wait[i] <- max(0, wait[i - 1] + service[i - 1] - arrival_gap) # Lindley
    time_in_system[i] <- wait[i] + service[i]
    departure[i] <- clock + time_in_system[i]
  }

  after_warmup <- time_in_system[-(1:r0)] # drop the warm-up customers

  return(list(output = mean(after_warmup),
              n_interarrival = n_interarrival,
              sum_interarrival = sum(interarrival[1:n_interarrival]),
              n_service = r,
              sum_service = sum(service)))
}

queue_replications <- function(lambda, mu, k, r0, r) {
  # Independent replications at the parameters (lambda[i], mu[i])
  # Inputs:
  # -- lambda, mu: vectors of rates, one entry per replication
  # -- k, r0, r: as in queue_replication()
  # Output:
  # -- data frame, one row per replication, columns as the list returned by
  #    queue_replication()

  n <- length(lambda) # number of replications
  result <- data.frame(output = numeric(n), n_interarrival = numeric(n),
                       sum_interarrival = numeric(n), n_service = numeric(n),
                       sum_service = numeric(n)) # placeholder
  for (i in 1:n) {
    replication <- queue_replication(lambda[i], mu[i], k, r0, r)
    result$output[i] <- replication$output
    result$n_interarrival[i] <- replication$n_interarrival
    result$sum_interarrival[i] <- replication$sum_interarrival
    result$n_service[i] <- replication$n_service
    result$sum_service[i] <- replication$sum_service
  }

  return(result)
}

# ===== Likelihoods and parameter estimation
log_path_density <- function(paths, lambda, mu) {
  # Log density of each simulated sample path under the rates (lambda, mu)
  # Inputs:
  # -- paths: data frame from queue_replications(), one row per path
  # -- lambda, mu: vectors of rates, one entry per parameter
  # Output:
  # -- matrix, [i, j] = log density of path i under parameter j
  # formula: n_service log(mu) - mu sum_service
  #          + n_interarrival log(lambda) - lambda sum_interarrival
  # (outer products for efficiency: all paths and parameters at once)

  log_density <- (outer(paths$n_service, log(mu)) -
                    outer(paths$sum_service, mu)) +
    (outer(paths$n_interarrival, log(lambda)) -
       outer(paths$sum_interarrival, lambda))

  return(log_density)
}

estimate_rates <- function(rate, n_obs, n) {
  # Maximum likelihood estimates of an exponential rate, each from n_obs
  # observations, used for the "real-world" data and for the parametric
  # bootstrap
  # Inputs:
  # -- rate: true (or estimated) rate of the exponential distribution
  # -- n_obs: number of observations per estimate (m in the paper)
  # -- n: number of estimates
  # Output:
  # -- vector of n estimated rates, 1 / (sample mean)

  observations <- matrix(rexp(n_obs * n, rate), n_obs, n) # column = data set

  return(1 / colMeans(observations))
}
