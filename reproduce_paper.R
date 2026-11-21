################################################################################
# reproduce_paper.R                                                            #
#                                                                              #
# Reproduces Table 1 of                                                        #
#   Feng, B. M. and Song, E. (2019). Efficient input uncertainty               #
#   quantification via green simulation using sample path likelihood ratios.   #
#   Proceedings of the 2019 Winter Simulation Conference, 3693-3704.           #
#                                                                              #
# Run from this folder:   Rscript reproduce_paper.R                            #
#                                                                              #
# By default (RERUN_SIMULATIONS = FALSE) Table 1 is computed from the saved    #
# per-macro results of the November 2026 rerun (results/rerun/) and checked    #
# against the paper (verify_against_paper.R). The results of the paper's own   #
# 2019 runs were not kept.                                                     #
#                                                                              #
# With RERUN_SIMULATIONS = TRUE the 10,000 macro replications are run again    #
# (about 1 hour on 24 cores, see README) into results/rerun_new/.              #
################################################################################

RERUN_SIMULATIONS <- FALSE
N_CORES <- max(1, parallel::detectCores() - 2) # used only for reruns

source(file.path("R", "utilities.R"))
source_library("R")
rscript <- file.path(R.home("bin"), "Rscript")

results_dir <- file.path("results", "rerun")
if (RERUN_SIMULATIONS) {
  results_dir <- file.path("results", "rerun_new")
  command <- paste(shQuote(rscript), "scripts/table1_simulate.R", 10000,
                   N_CORES, results_dir)
  message("Running: ", command)
  if (system(command) != 0) {
    stop("table1_simulate.R failed")
  }
}

dir.create("tables", showWarnings = FALSE)
t1 <- table1(results_dir)
write.csv(t1, file.path("tables", "Table1_Coverage_Width_Runtime.csv"),
          row.names = FALSE)
cat("\n==== Table 1 ====\n")
shown <- t1
shown$coverage <- sprintf("%.2f%% (%.2f%%)", 100 * t1$coverage,
                          100 * t1$coverage_se)
shown$width <- sprintf("%.2f (%.3f)", t1$width, t1$width_se)
print(shown[, c("B", "method", "R", "coverage", "width",
                "seconds_simulation", "seconds_lr")], row.names = FALSE)

system(paste(shQuote(rscript), "verify_against_paper.R", results_dir))
