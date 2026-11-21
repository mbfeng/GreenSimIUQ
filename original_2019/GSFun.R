MM1kIUQ_Qtl_Calc <- function(lambda_B, mu_B, k, alpha){
    # Calculates sample quantile confidence interval for time-in-system (TIS) 
    # in MM1k queue using analytical formula
    # Inputs:
    # -- lambda_B: vector of length n_B, bootstrapped interarrival rates
    # -- mu_B: vector of lenght n_B, bootstrapped service rates
    # -- k: integer, system capacity
    # -- alpha: confidence level of the IUQ confidence interval
    # Outputs: a list of 2 items
    # -- CI: confidence interval, a vector of 2 values
    # -- TIS_B: the estimated mean for the n_B bootstrapped parameters
    
    n_B <- length(lambda_B) # no. of bootstrapped params
    TIS_B <- rep(0, n_B) # placeholder for TIS of given params
    
    # inner simulations
    for (b in 1:n_B)
    { 
        TIS_B[b] <- MM1kTIS_Calc(lambda_B[b], mu_B[b], k) 
    }
    
    # quantiles based on linear interpolation of the empirical cdf
    CI <- quantile(TIS_B, probs = c(alpha/2, 1-alpha/2), type = 4)
    
    return(list(CI = CI, TIS_B = TIS_B))
}

MM1kIUQ_Qtl_Nested <- function(lambda_B, mu_B, k, r0, r, alpha, n_R){
    # Calculates sample quantile confidence interval for time-in-system (TIS) 
    # in MM1k queue using nested sim: Lindley equation as inner sim
    # Inputs:
    # -- lambda_B: vector of length n_B, bootstrapped interarrival rates
    # -- mu_B: vector of lenght n_B, bootstrapped service rates
    # -- k: integer, system capacity
    # -- alpha: confidence level of the IUQ confidence interval
    # -- n_R: number of i.i.d. sample paths
    # Outputs: a list of 2 items
    # -- CI: confidence interval, a vector of 2 values
    # -- TIS_B: the estimated mean for the n_B bootstrapped parameters
    
    n_B <- length(lambda_B)     # number of bootstrapped distributions
    TIS_B <- rep(0, n_B)   # placeholder for TIS of given params
    inner_sims <- rep(0, n_R)   # placeholder for inner sim. outputs
    
    for(b in 1:n_B){
        lambda_b = lambda_B[b]
        mu_b = mu_B[b]
        
        for(i in 1:n_R){
            inner_sims[i] <- MM1kTIS_Lindley(lambda_b, mu_b, k, r0, r)$Ybar
        }
        TIS_B[b] <- mean(inner_sims)
    }
    
    # quantiles based on linear interpolation of the empirical cdf
    CI <- quantile(TIS_B, probs = c(alpha/2, 1-alpha/2), type = 4)
    
    return(list(CI = CI, TIS_B = TIS_B))
}

MM1kIUQ_Qtl_MIS <- function(lambda_B, mu_B, lambda_S, mu_S,
                            k, r0, r, alpha, n_R, selfnorm = F){
    # Calculates sample quantile confidence interval for time-in-system (TIS) 
    # in MM1k queue using multiple importance sampling (MIS)
    # Inputs:
    # -- lambda_B: vector of length n_B, bootstrapped interarrival rates
    # -- mu_B: vector of lenght n_B, bootstrapped service rates
    # -- lambda_S: vector of length n_S, sampling interarrival rate(s)
    # -- mu_S: vector of lenght n_S, sampling service rate(s)
    # -- k: integer, system capacity
    # -- alpha: confidence level of the IUQ confidence interval
    # -- n_R: number of i.i.d. sample paths
    # -- selfnorm (optional): indicator for self-normalization, default True
    # Outputs: a list of 2 items
    # -- CI: confidence interval, a vector of 2 values
    # -- TIS_B: the estimated mean for the n_B bootstrapped parameters
    
    n_B <- length(lambda_B)     # number of bootstrapped distributions
    n_S <- length(lambda_S)     # number of bootstrapped distributions
    TIS_B <- rep(0, n_B)   # placeholder for TIS of given params
    
    lambda_S <- rep(lambda_S, n_R)
    mu_S <- rep(mu_S, n_R)
    n_R_S <- n_R * n_S
    
    Ybar <- rep(0, n_R_S)
    n_Arr <- rep(0, n_R_S)
    sum_Arr <- rep(0, n_R_S)
    n_Sr <- rep(0, n_R_S)
    sum_Sr <- rep(0, n_R_S)
    
    # simulates the inputs/outputs for likelihood ratio calculations
    for(i in 1:n_R_S){
        Out_i <- MM1kTIS_Lindley(lambda_S[i], mu_S[i], k, r0, r)
        Ybar[i] <- Out_i$Ybar
        n_Arr[i] <- Out_i$n_Arr
        sum_Arr[i] <- Out_i$sum_Arr
        n_Sr[i] <- Out_i$n_Sr  
        sum_Sr[i] <- Out_i$sum_Sr
    }
    
    # calculate the exponents needed for the denominator of the LR
    LR_denom <- matrix(0, nrow = n_R_S, ncol = n_S)
    for(s in 1:n_S){
        LR_denom[, s] <- MM1kPdfExp(n_Arr, sum_Arr, n_Sr, sum_Sr, 
                                    lambda_S[s], mu_S[s])
    }
    
    # estimate the TIS for each pair of bootstrapped params
    for (b in 1:n_B)
    {
        lambda_b <- lambda_B[b]
        mu_b <- mu_B[b]
        
        LR_num <- MM1kPdfExp(n_Arr, sum_Arr, n_Sr, sum_Sr,
                             lambda_b, mu_b)
        LR <- 1/rowMeans(exp(LR_denom - LR_num))
        
        TIS_B[b] <- mean(Ybar*LR)  # MIS estimator
        if(selfnorm){
            # self-normalization to ensure average sample LR = 1
            TIS_B[b] <- TIS_B[b]/mean(LR)
        }
    }
    
    # quantiles based on linear interpolation of the empirical cdf
    CI <- quantile(TIS_B, probs = c(alpha/2, 1-alpha/2), type = 4)
    
    return(list(CI = CI, TIS_B = TIS_B))
}

MM1kIUQ_Shr_MIS <- function(lambda_B, mu_B, lambda_S, mu_S, lambda_C, mu_C,
                            k, r0, r, alpha, n_R, selfnorm = F){
    # Calculates shrinkage confidence interval for time-in-system (TIS) 
    # in MM1k queue using multiple importance sampling (MIS)
    # Inputs:
    # -- lambda_B: vector of length n_B, bootstrapped interarrival rates
    # -- mu_B: vector of lenght n_B, bootstrapped service rates
    # -- lambda_S: vector of length n_S, sampling interarrival rate(s)
    # -- mu_S: vector of lenght n_S, sampling service rate(s)
    # -- lambda_C: constant, estimated interarrival rate(s) from data
    # -- mu_C: constant, estimated service rate(s) from data
    # -- k: integer, system capacity
    # -- alpha: confidence level of the IUQ confidence interval
    # -- n_R: number of i.i.d. sample paths
    # -- selfnorm (optional): indicator for self-normalization, default True
    # Outputs: a list of 2 items
    # -- CI: confidence interval, a vector of 2 values
    # -- TIS_B: the estimated mean for the n_B bootstrapped parameters
    
    n_B <- length(lambda_B)     # number of bootstrapped distributions
    n_S <- length(lambda_S)     # number of bootstrapped distributions
    TIS_B <- rep(0, n_B)        # placeholder for TIS of given params
    S2_B  <- rep(0, n_B)        # sample var of TIS of given params
    
    lambda_S <- rep(lambda_S, n_R)
    mu_S <- rep(mu_S, n_R)
    n_R_S <- n_R * n_S
    
    Ybar <- rep(0, n_R_S)
    n_Arr <- rep(0, n_R_S)
    sum_Arr <- rep(0, n_R_S)
    n_Sr <- rep(0, n_R_S)
    sum_Sr <- rep(0, n_R_S)
    
    # simulates the inputs/outputs for likelihood ratio calculations
    for(i in 1:n_R_S){
        Out_i <- MM1kTIS_Lindley(lambda_S[i], mu_S[i], k, r0, r)
        Ybar[i] <- Out_i$Ybar
        n_Arr[i] <- Out_i$n_Arr
        sum_Arr[i] <- Out_i$sum_Arr
        n_Sr[i] <- Out_i$n_Sr  
        sum_Sr[i] <- Out_i$sum_Sr
    }
    
    # calculate the exponents needed for the denominator of the LR
    LR_denom <- matrix(0, nrow = n_R_S, ncol = n_S)
    for(s in 1:n_S){
        LR_denom[, s] <- MM1kPdfExp(n_Arr, sum_Arr, n_Sr, sum_Sr, 
                                    lambda_S[s], mu_S[s])
    }
    
    # estimate the TIS for each pair of bootstrapped params
    for (b in 1:n_B)
    {
        lambda_b <- lambda_B[b]
        mu_b <- mu_B[b]
        
        LR_num <- MM1kPdfExp(n_Arr, sum_Arr, n_Sr, sum_Sr,
                             lambda_b, mu_b)
        LR <- 1/rowMeans(exp(LR_denom - LR_num))
        YLR <- Ybar*LR
        
        TIS_B[b] <- mean(YLR)   # MIS estimator
        S2_B[b] <- var(YLR)     # sample variance of MIS estimator
        if(selfnorm){
            # self-normalization to ensure average sample LR = 1
            TIS_B[b] <- TIS_B[b]/mean(LR)
        }
    }
    
    Vhat <- var(TIS_B)
    sighatI2 <- max(Vhat - mean(S2_B)/n_R_S, 0)
    
    # estimate the TIS for theta_C, param based on observed data
    LR_num <- MM1kPdfExp(n_Arr, sum_Arr, n_Sr, sum_Sr,
                         lambda_C, mu_C)
    LR <- 1/rowMeans(exp(LR_denom - LR_num))
    YLR <- Ybar*LR
    TIS_C <- mean(YLR)   # MIS estimator for muhat(thetahat)
    S2_C <- var(YLR)     # sample variance of MIS estimator
    
    # quantiles based on linear interpolation of the empirical cdf
    Z <- qnorm(1 - alpha/2)
    CI <- c(TIS_C-Z*sqrt(sighatI2 + S2_C/n_R_S), TIS_C+Z*sqrt(sighatI2 + S2_C/n_R_S))
    
    return(list(CI = CI, TIS_B = TIS_B))
}

MM1k_RunMacros_MIS <- function(lambda = 4.5, mu = 5, k = 10, r0 = 50, r = 200,
                               n_samp = 5e1, n_B = 1e1, n_R_GS = 1e2, n_S = 1e1,
                               alpha = 0.05, n_macro = 1e2, seed = NULL, 
                               type = "Quantile"){
    
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
    
    if(n_S == 1){
        print(sprintf("Running ILR with %s", type))
        lambda_S <- matrix(MM1kParam_Thres(lambda_C, lambda_B), ncol = n_macro)
        mu_S <- matrix(MM1kParam_Thres(mu_C, mu_B), ncol = n_macro)
    }
    else{
        print(sprintf("Running MLR with %s", type))
        lambda_S <- MM1kParam_ResampThres(lambda_B, n_S)
        mu_S <- MM1kParam_ResampThres(mu_B, n_S)
        n_R_GS <- n_R_GS / n_S
    }
    
    zero_vec <- rep(0, n_macro)
    IUQ_Temp <- data.frame(CI_LB = zero_vec, CI_UB = zero_vec)
    # TIS_Temp <- matrix(0, n_B, n_macro)  # placeholder for results
    Time_Temp <- rep(0, n_macro)         # placeholder for runtimes
    
    start_time <- Sys.time()
    for(i in 1:n_macro) # each macro rep uses a different real-world sample
    {
        if(i %% 10 == 0){   print(sprintf("Running the %i-th macro rep.", i))  }
        start_time_per_macro <- Sys.time()
        
        if(type == "Quantile"){
            Out_i <- MM1kIUQ_Qtl_MIS(lambda_B[, i], mu_B[, i],
                                     lambda_S[, i], mu_S[, i],
                                     k, r0, r, alpha, n_R_GS,
                                     selfnorm = F)
        }else if(type == "Shrinkage"){
            Out_i <- MM1kIUQ_Shr_MIS(lambda_B[, i], mu_B[, i], 
                                     lambda_S[, i], mu_S[, i],
                                     lambda_C[i], mu_C[i],
                                     k, r0, r, alpha, n_R_GS,
                                     selfnorm = F)
        }
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
    
    return(list(Avg = Avg_Temp, IUQ = IUQ_Temp))
}

