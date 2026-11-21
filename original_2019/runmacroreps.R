rm(list = ls())
source("MM1kFun.R")
source("GSFun.R")

MM1k_RunMacros <- function(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                           n_samp = 5e1, n_B = 1e1, n_S = 1e1,
                           n_R_Inner = 1e2, n_R_GS = 1e2, 
                           alpha = 0.05, n_macro = 1e2, seed = NULL, 
                           types = c("NoInner", "WithInner", "ILR", "MLR"),
                           selfnorm = F){
    # run macro reps to calculate IUQ CIs for the MM1k example using 
    # different methods as specified in "types"
    # Inputs:
    # -- lambda: real interarrival rate in the MM1k queue
    # -- mu: real service rate in the MM1k queue
    # -- k: system capacity in the MM1k queue
    # -- alpha: confidence level of the IUQ confidence interval
    # -- r0: warm-up observations
    # -- r: total observations (including warm-up)
    # -- n_samp: no. of real-world/bootstrap service & interarrival time obs.
    # -- n_B: no. of bootstrapped distributions
    # -- n_S: no. of sampling distributions, should be a divisor of n_R_GS
    # -- n_R_Inner: no. of replications per bootstrapped distribution
    # -- n_R_GS: no. of replications in total for GS estimators
    # -- n_macro: no. of macro replications
    # -- seed(optional): random seed for reproducibility
    
    # design of experiment (DOE) parameters 
    DOE <- list(lambda = lambda, mu = mu, k = k, alpha = alpha, r0 = r0, r = r, 
                n_samp = n_samp, n_B = n_B, n_R_Inner = n_R_Inner, 
                n_R_GS = n_R_GS, n_macro = n_macro, seed = seed, types = types,
                selfnorm = selfnorm)
    
    # IUQ results
    IUQ <- vector("list", length(types))
    names(IUQ) <- types
    
    #---- Analytical solution for time-in-system (TIS), for benchmarking
    TIS_True <- MM1kTIS_Calc(lambda, mu, k) # analytical solution
    
    #---- Presimulate common parameters for all methods
    set.seed(seed)  # set seed for reproducibility
    lambda_C <- MM1kParam_Bootstrap(lambda, n_samp, n_macro)  # arrival rates est. from real-world obs.
    mu_C <- MM1kParam_Bootstrap(mu, n_samp, n_macro)  # service rates est. from real-world obs.
    
    lambda_B <- matrix(0, n_B, n_macro)    # bootstrapped arrival rates
    mu_B <- matrix(0, n_B, n_macro)        # bootstrapped service rates
    
    for(i in 1:n_macro) # bootstrapped params in each macro rep
    {
        lambda_B[, i] <- MM1kParam_Bootstrap(lambda_C[i], n_samp, n_B)
        mu_B[, i] <- MM1kParam_Bootstrap(mu_C[i], n_samp, n_B)
    }

        #---- IUQ without inner simulation
    if("NoInner" %in% types){
        print("Running NoInner.")
        zero_vec <- rep(0, n_macro)
        IUQ_Temp <- data.frame(CI_LB = zero_vec, CI_UB = zero_vec)
        # TIS_Temp <- matrix(0, n_B, n_macro)  # placeholder for results
        Time_Temp <- rep(0, n_macro)         # placeholder for runtimes
        
        start_time <- Sys.time()
        for(i in 1:n_macro) # each macro rep uses a different real-world sample
        {
            if(i %% 10 == 0){   print(sprintf("Running the %i-th macro rep.", i))  }
            start_time_per_macro <- Sys.time()
            
            Out_i <- MM1kIUQ_Qtl_Calc(lambda_B[, i], mu_B[, i], k, alpha)
            IUQ_Temp[i, ] <- Out_i$CI
            # TIS_Temp[, i] <- Out_i$TIS_Calc
            
            end_time_per_macro <- Sys.time()
            Time_Temp[i] <- difftime(end_time_per_macro, start_time_per_macro,
                                     units = "secs")
        }
        
        IUQ_Temp$CI_width <- IUQ_Temp$CI_UB - IUQ_Temp$CI_LB
        IUQ_Temp$coverage <- (IUQ_Temp$CI_LB<=TIS_True)*(IUQ_Temp$CI_UB>=TIS_True) # coverage indicator
        IUQ_Temp$RunTime <- Time_Temp
        Avg_Temp <- colMeans(IUQ_Temp)
        end_time <- Sys.time()
        
        print(sprintf("NoInner completed in %5.3f seconds.", 
                      difftime(end_time, start_time,units = "secs")))
        print(sprintf("Empirical coverage is %5.3f, CI width is %5.3f.",
                      Avg_Temp["coverage"], Avg_Temp["CI_width"]))
        
        IUQ$NoInner <- list(Avg_NoInner = Avg_Temp, IUQ_NoInner = IUQ_Temp)
    }
    
    #---- IUQ with inner simulation
    if("WithInner" %in% types){
        print("Running WithInner.")
        zero_vec <- rep(0, n_macro)
        IUQ_Temp <- data.frame(CI_LB = zero_vec, CI_UB = zero_vec)
        # TIS_Temp <- matrix(0, n_B, n_macro)  # placeholder for results
        Time_Temp <- rep(0, n_macro)         # placeholder for runtimes
        
        start_time <- Sys.time()
        for(i in 1:n_macro) # each macro rep uses a different real-world sample
        {
            if(i %% 10 == 0){   print(sprintf("Running the %i-th macro rep.", i))  }
            start_time_per_macro <- Sys.time()
            
            Out_i <- MM1kIUQ_Qtl_Nested(lambda_B[, i], mu_B[, i], k, r0, r, alpha, n_R_Inner)
            IUQ_Temp[i, ] <- Out_i$CI
            # TIS_Temp[, i] <- Out_i$TIS_B
            
            end_time_per_macro <- Sys.time()
            Time_Temp[i] <- difftime(end_time_per_macro, start_time_per_macro,
                                     units = "secs")
        }
        
        IUQ_Temp$CI_width <- IUQ_Temp$CI_UB - IUQ_Temp$CI_LB
        IUQ_Temp$coverage <- (IUQ_Temp$CI_LB<=TIS_True)*(IUQ_Temp$CI_UB>=TIS_True) # coverage indicator
        IUQ_Temp$RunTime <- Time_Temp
        Avg_Temp <- colMeans(IUQ_Temp)
        end_time <- Sys.time()
        
        print(sprintf("WithInner completed in %5.3f seconds.", 
                      difftime(end_time, start_time,units = "secs")))
        print(sprintf("Empirical coverage is %5.3f, CI width is %5.3f.",
                      Avg_Temp["coverage"], Avg_Temp["CI_width"]))
        
        IUQ$WithInner <- list(Avg_WithInner = Avg_Temp, IUQ_WithInner = IUQ_Temp)
    }
    
    if("ILR" %in% types){
        print("Running ILR.")
        if(selfnorm){
            print("Self-normalization applied.")
        }
        zero_vec <- rep(0, n_macro)
        IUQ_Temp <- data.frame(CI_LB = zero_vec, CI_UB = zero_vec)
        # TIS_Temp <- matrix(0, n_B, n_macro)  # placeholder for results
        Time_Temp <- rep(0, n_macro)         # placeholder for runtimes
        
        # parameters repeated from lambda_C and mu_C, but limit to thresholds
        # the thresholds make sure the ILR variance is finite
        lambda_S <- MM1kParam_Thres(lambda_C, lambda_B)
        mu_S <- MM1kParam_Thres(mu_C, mu_B)

        start_time <- Sys.time()
        for(i in 1:n_macro) # each macro rep uses a different real-world sample
        {
            if(i %% 10 == 0){   print(sprintf("Running the %i-th macro rep.", i))  }
            start_time_per_macro <- Sys.time()
            
            Out_i <- MM1kIUQ_Qtl_MIS(lambda_B[, i], mu_B[, i], 
                                     lambda_S[i], mu_S[i],
                                     k, r0, r, alpha, n_R_GS,
                                     selfnorm = selfnorm)
            IUQ_Temp[i, ] <- Out_i$CI
            # TIS_Temp[, i] <- Out_i$TIS_B
            
            end_time_per_macro <- Sys.time()
            Time_Temp[i] <- difftime(end_time_per_macro, start_time_per_macro,
                                     units = "secs")
        }
        
        IUQ_Temp$CI_width <- IUQ_Temp$CI_UB - IUQ_Temp$CI_LB
        IUQ_Temp$coverage <- (IUQ_Temp$CI_LB<=TIS_True)*(IUQ_Temp$CI_UB>=TIS_True) # coverage indicator
        IUQ_Temp$RunTime <- Time_Temp
        Avg_Temp <- colMeans(IUQ_Temp)
        end_time <- Sys.time()
        
        print(sprintf("ILR completed in %5.3f seconds.", 
                      difftime(end_time, start_time,units = "secs")))
        print(sprintf("Empirical coverage is %5.3f, CI width is %5.3f.",
                      Avg_Temp["coverage"], Avg_Temp["CI_width"]))
        
        IUQ$ILR <- list(Avg_ILR = Avg_Temp, IUQ_ILR = IUQ_Temp)
    }
    
    if("MLR" %in% types){
        print("Running MLR.")
        if(selfnorm){
            print("Self-normalization applied.")
        }
        zero_vec <- rep(0, n_macro)
        IUQ_Temp <- data.frame(CI_LB = zero_vec, CI_UB = zero_vec)
        # TIS_Temp <- matrix(0, n_B, n_macro)  # placeholder for results
        Time_Temp <- rep(0, n_macro)         # placeholder for runtimes
        
        # parameters resampled from lambda_B and mu_B, subject to thresholds
        lambda_S <- MM1kParam_ResampThres(lambda_B, n_S)
        mu_S <- MM1kParam_ResampThres(mu_B, n_S)
        n_R_MLR <- n_R_GS / n_S
        
        
        start_time <- Sys.time()
        for(i in 1:n_macro) # each macro rep uses a different real-world sample
        {
            if(i %% 10 == 0){   print(sprintf("Running the %i-th macro rep.", i))  }
            start_time_per_macro <- Sys.time()
            
            Out_i <- MM1kIUQ_Qtl_MIS(lambda_B[, i], mu_B[, i], 
                                     lambda_S[, i], mu_S[, i],
                                     k, r0, r, alpha, n_R_MLR,
                                     selfnorm = selfnorm)
            IUQ_Temp[i, ] <- Out_i$CI
            # TIS_Temp[, i] <- Out_i$TIS_B
            
            end_time_per_macro <- Sys.time()
            Time_Temp[i] <- difftime(end_time_per_macro, start_time_per_macro,
                                     units = "secs")
        }
        
        IUQ_Temp$CI_width <- IUQ_Temp$CI_UB - IUQ_Temp$CI_LB
        IUQ_Temp$coverage <- (IUQ_Temp$CI_LB<=TIS_True)*(IUQ_Temp$CI_UB>=TIS_True) # coverage indicator
        IUQ_Temp$RunTime <- Time_Temp
        Avg_Temp <- colMeans(IUQ_Temp)
        end_time <- Sys.time()
        
        print(sprintf("MLR completed in %5.3f seconds.", 
                      difftime(end_time, start_time,units = "secs")))
        print(sprintf("Empirical coverage is %5.3f, CI width is %5.3f.",
                      Avg_Temp["coverage"], Avg_Temp["CI_width"]))
        
        IUQ$MLR <- list(Avg_MLR = Avg_Temp, IUQ_MLR = IUQ_Temp)
    }
    
    filename <- sprintf("%slambda%smu%sOuter%sInner%sGS%sSamp%s.RData",
                        format(Sys.Date(),"%Y%m%d"),
                        lambda, mu, n_B, n_R_Inner, n_R_GS, n_S)

    save(DOE, IUQ, file = filename)
    
    return(list(DOE, IUQ))
}

# Experiment 1: Fixed simulation budget: n_B * n_R_Inner = 1e3
# -- 1. n_B = 1e3, n_R_Inner = 1e0
# -- 2. n_B = 1e2, n_R_Inner = 1e1
# -- 3. Add n_R_GS = 1e3 to both above so computational budget is comparable
# -- 4. Add n_R_GS = 1e2 to both above so computational budget is less
Outer1e2GS1e2 <- MM1k_RunMacros(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                                n_samp = 5e1, n_B = 1e1, n_S = 1e1,
                                n_R_Inner = 1e2, n_R_GS = 1e2, 
                                alpha = 0.05, n_macro = 1e2, seed = NULL, 
                                types = c("ILR", "MLR"), selfnorm = F)
rm(Outer1e2GS1e2)

Outer1e3GS1e2 <- MM1k_runmacro(lambda = 4.5, mu = 5, k = 10, alpha = 0.05,
                               r0 = 50, r = 100,
                               n_samp = 5e1, n_B = 1e3, n_R_Inner = 0,
                               n_R_GS = 1e2, n_macro = 1e3,
                               seed = 0, types = c("ILR", "MLR"))
rm(Outer1e3GS1e2)

Outer1e2Inner1e1GS1e3 <- MM1k_runmacro(lambda = 4.5, mu = 5, k = 10, alpha = 0.05,
                                       r0 = 50, r = 100,
                                       n_samp = 5e1, n_B = 1e2, n_R_Inner = 1e1,
                                       n_R_GS = 1e3, n_macro = 1e3,
                                       seed = 0, types = c("NoInner", "WithInner", "ILR", "MLR"))
rm(Outer1e2Inner1e1GS1e3)

Outer1e3Inner1e0GS1e3 <- MM1k_runmacro(lambda = 4.5, mu = 5, k = 10, alpha = 0.05,
                                       r0 = 50, r = 100,
                                       n_samp = 5e1, n_B = 1e3, n_R_Inner = 1e0,
                                       n_R_GS = 1e3, n_macro = 1e3,
                                       seed = 0, types = c("NoInner", "WithInner", "ILR", "MLR"))
rm(Outer1e3Inner1e0GS1e3)
