################################################################################
# verify_against_paper.R                                                       #
#                                                                              #
# Compares the coverage probabilities and CI widths of Table 1 (the paper's    #
# printed values in verify/paper_table1.csv) with those computed from a        #
# results folder.                                                              #
#                                                                              #
#   Rscript verify_against_paper.R [results_dir]     (default results/rerun)   #
#                                                                              #
# The reruns use different random numbers from the paper's runs, so each       #
# entry must agree within Monte Carlo error:                                   #
#   z = (x - v) / sqrt(se_x^2 + se_v^2 + h^2 / 3),  |z| <= 2.5,                #
# where x is the computed value with standard error se_x, v the paper's        #
# value, se_v = se_x * sqrt(K_x / 10000) the Monte Carlo error of the paper's  #
# value (10,000 macro runs), and h half a unit in the last printed digit.      #
# Run times depend on the computer and are not compared.                       #
################################################################################

args <- commandArgs(trailingOnly = TRUE)
results_dir <- file.path("results", "rerun")
if (length(args) >= 1) {
  results_dir <- args[1]
}

source(file.path("R", "utilities.R"))
source_library("R")

N_MACRO_PAPER <- 10000
paper <- read.csv(file.path("verify", "paper_table1.csv"),
                  colClasses = c("numeric", "character", "numeric",
                                 "character", "character"))
computed <- table1(results_dir)

half_unit <- function(printed) {
  # Half a unit in the last printed digit (the rounding error of the paper)
  # Input:
  # -- printed: the paper's value as text, e.g. "89.97"
  # Output:
  # -- half a unit in the last digit, e.g. 0.005

  decimals <- 0 # number of printed decimals
  if (grepl("\\.", printed)) {
    decimals <- nchar(sub(".*\\.", "", printed))
  }

  return(0.5 * 10^(-decimals))
}

checks <- list()
for (i in seq_len(nrow(paper))) {
  entry <- paper[i, ]
  row <- computed[computed$B == entry$B & computed$method == entry$method &
                    computed$R == entry$R, ]
  for (quantity in c("coverage", "width")) {
    printed <- entry[[quantity]]
    scale <- if (quantity == "coverage") 100 else 1 # coverage printed in %
    x <- row[[quantity]] * scale
    se_x <- row[[paste0(quantity, "_se")]] * scale
    se_paper <- se_x * sqrt(row$n_macro / N_MACRO_PAPER)
    z <- (x - as.numeric(printed)) /
      sqrt(se_x^2 + se_paper^2 + half_unit(printed)^2 / 3)
    checks[[length(checks) + 1]] <- data.frame(
      B = entry$B, method = entry$method, R = entry$R, quantity = quantity,
      paper = printed, computed = signif(x, 5), se = signif(se_x, 3),
      z = round(z, 2), status = if (abs(z) <= 2.5) "ok" else "DIFFERENT")
  }
}
report <- stack_rows(checks)
label <- basename(normalizePath(results_dir))
write.csv(report, file.path("verify", paste0("verification_", label, ".csv")),
          row.names = FALSE)
cat(sprintf("\nVerification of %s: %d checks\n", results_dir, nrow(report)))
print(table(report$status))
print(report, row.names = FALSE)
cat(sprintf(paste("\nlargest |z| = %.2f; expected number of |z| > 2.5 by",
                  "chance (no real difference): %.1f\n"),
            max(abs(report$z)), 0.0124 * nrow(report)))
