################################################################################
# tables.R                                                                     #
# Table 1 of the paper from the per-macro results: for every number of         #
# bootstraps B, method and R, the empirical coverage probability and average   #
# width of the CIs, with standard errors sd / sqrt(K) (K = number of macro     #
# runs), and the average run time per macro run                                #
################################################################################

table1 <- function(results_dir) {
  # Table 1 from a results folder
  # Input:
  # -- results_dir: folder with table1_per_macro.csv.gz
  # Output:
  # -- data frame, one row per (B, method, R) in the order of the paper

  per_macro <- read.csv(file.path(results_dir, "table1_per_macro.csv.gz"))
  per_macro$R[is.na(per_macro$R)] <- 0 # Oracle has no replications

  rows <- list() # placeholder, one data frame per table row
  for (n_boot in sort(unique(per_macro$B), decreasing = TRUE)) {
    for (method in c("Oracle", "SMC", "MLR", "ILR")) {
      runs_method <- per_macro[per_macro$B == n_boot &
                                 per_macro$method == method, ]
      for (n_rep in sort(unique(runs_method$R), decreasing = TRUE)) {
        runs <- runs_method[runs_method$R == n_rep, ] # all macro runs
        rows[[length(rows) + 1]] <- data.frame(
          B = n_boot, method = method, R = n_rep, n_macro = nrow(runs),
          coverage = mean(runs$covered), # empirical coverage probability
          coverage_se = standard_error(runs$covered),
          width = mean(runs$width), # average CI width
          width_se = standard_error(runs$width),
          seconds_simulation = mean(runs$seconds_simulation),
          seconds_lr = mean(runs$seconds_lr)
        )
      }
    }
  }

  return(stack_rows(rows))
}
