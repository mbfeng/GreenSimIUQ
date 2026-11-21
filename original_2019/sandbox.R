rm(list = ls())
source("MM1kFun.R")
source("GSFun.R")

# # simulator parameters
# lambda = 4.5
# mu = 5
# k = 10
# alpha = 0.05
# r0 = 50
# r = 200
# 
# # design of experiment
# n_samp = 5e1    # sample size of observations (also bootstrapped sample size)
# n_B = 1e1       # number of bootstrapped input distributions
# n_R_Inner = 1e2 # number of inner sample paths (for nested sim)
# n_R_GS = 1e2    # total number of sample paths in green simulation
# n_S = 1e1       # number of sampling input distributions (divisor of n_R_GS)
# n_macro = 1e2   # number of macro replications
# seed = NULL     # random seed for reproducibility
# types = c("NoInner", "WithInner", "ILR", "MLR")
# selfnorm = F

IUQ1 <- MM1k_RunMacros_MIS(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                               n_samp = 5e1, n_B = 1e2, n_R_GS = 1e2, n_S = 1,
                               alpha = 0.05, n_macro = 1e3, seed = NULL, 
                               type = "Quantile")
    
IUQ2 <- MM1k_RunMacros_MIS(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                           n_samp = 5e1, n_B = 1e2, n_R_GS = 1e2, n_S = 1,
                           alpha = 0.05, n_macro = 1e3, seed = NULL, 
                           type = "Shrinkage")

IUQ3 <- MM1k_RunMacros_MIS(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                           n_samp = 5e1, n_B = 1e2, n_R_GS = 1e3, n_S = 1,
                           alpha = 0.05, n_macro = 1e3, seed = NULL, 
                           type = "Quantile")

IUQ4 <- MM1k_RunMacros_MIS(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                           n_samp = 5e1, n_B = 1e2, n_R_GS = 1e3, n_S = 1,
                           alpha = 0.05, n_macro = 1e3, seed = NULL, 
                           type = "Shrinkage")

IUQ5 <- MM1k_RunMacros_MIS(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                           n_samp = 5e1, n_B = 1e2, n_R_GS = 1e2, n_S = 1e1,
                           alpha = 0.05, n_macro = 1e3, seed = NULL, 
                           type = "Quantile")

IUQ6 <- MM1k_RunMacros_MIS(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                           n_samp = 5e1, n_B = 1e2, n_R_GS = 1e2, n_S = 1e1,
                           alpha = 0.05, n_macro = 1e3, seed = NULL, 
                           type = "Shrinkage")

IUQ7 <- MM1k_RunMacros_MIS(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                           n_samp = 5e1, n_B = 1e2, n_R_GS = 1e3, n_S = 1e1,
                           alpha = 0.05, n_macro = 1e3, seed = NULL, 
                           type = "Quantile")

IUQ8 <- MM1k_RunMacros_MIS(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                           n_samp = 5e1, n_B = 1e2, n_R_GS = 1e3, n_S = 1e1,
                           alpha = 0.05, n_macro = 1e3, seed = NULL, 
                           type = "Shrinkage")
    
save.image(file = "IUQ_Quantile_Shrinkage.RData")
    