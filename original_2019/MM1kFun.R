#############
# MM1kFun.R #
#############
# Notations
# -- lambda: mean interarrival rate (exponential distribution)
# -- mu: mean service rate (exponential distribution)
# -- k: system capacity
# -- r0: warm-up observations
# -- r: total observations (including warm-up)

MM1kTIS_Calc <- function(lambda, mu, k){
    # Calculates the steady-state time-in-system (TIS) in an M/M/1/k queue with 
    # interarrival rate lambda, service rate mu, and capacity k.
    # Inputs:
    # -- lambda: mean interarrival rate (exponential distribution)
    # -- mu: mean service rate (exponential distribution)
    # -- k: system capacity
    # Outputs:
    # -- TIS: constant, long-run mean time-in-system
    
    # Steady-State statistics
    rho <- lambda/mu
    NIS <- rho/(1-rho) - (k+1)*rho^(k+1)/(1-rho^(k+1))   # number in system
    pk <- rho^k/((1-rho^(k+1))/(1-rho))                  # (Eunhye) ???
    TIS <- NIS/(lambda*(1-pk))                           # time in system
    
    return(TIS)
}

MM1kTIS_Lindley <- function(lambda, mu, k, r0, r){
    # Simulates the steady-state time-in-system (TIS) in an M/M/1/k queue using 
    # Lindley's equation with interarrival rate lambda, service rate mu,
    # capacity k, number of warm-up obs r0, and total number of obs r>r0.
    # Inputs:
    # -- lambda: mean interarrival rate (exponential distribution)
    # -- mu: mean service rate (exponential distribution)
    # -- k: system capacity
    # -- r0: warm-up observations
    # -- r: total observations (including warm-up)
    # Outputs: a list of 7 items
    # -- data: simulated TIS after warmup
    # -- Ybar: estimated mean TIS, sample average of data
    # -- S2: estimated variance of TIS, sample variance of data
    # -- n_Arr: number of simulated interarrival times, including warm ups
    # -- sum_Arr: sum of simulated interarrival times, including warm ups
    # -- n_Sr: number of simulated service times
    # -- sum_Sr: sum of simulated service times
    
    # placeholders for results, note that only r cumstomers receive services
    Sr <- rexp(r, mu)       # service times for r customers
    W <- matrix(0, r, 1)    # wait time for r custerms
    TIS <- matrix(0, r, 1)  # time-in-system for r customers
    D <- matrix(0, r, 1)    # departure time for r customers
    Arr <- matrix(0, r, 1)  # arrival times for r customers
    
    prevAi <- 0     # system initialization, no previous arrival
    departed <- 0   # system initialization, no previous departure
    
    cur_arr <- rexp(1, lambda)  # first arrival, enters the system immediately
    clock <- cur_arr            # event clock: arrivals, services, & departures
    All_arr <- cur_arr          # all arrivals
    Arr[1, ] <- cur_arr         # arrivals of serviced customers
    W[1, ] <- 0                 # wait times of serviced customers
    TIS[1, ] <- Sr[1]           # time-in-system of serviced customers
    D[1, ] <- clock + Sr[1]     # departure times of serviced customers
    NIS <- 1                    # number of customers in system
    
    for(i in 2:r)
    {
        repeat  # repeat until 1 arrival is accepted
        {
            cur_arr <- rexp(1, lambda)      # time-to-next arrival
            All_arr <- c(All_arr, cur_arr)  # time-to-next arrival
            Ai <- prevAi + cur_arr          # arrival time
            if(clock + Ai > D[(i-1), ])     # i-th arrival finds system empty
            {  j = 0; break; }              # empty system
            
            for(j in NIS:1)
            {
                # count (j-1) departures in time period of Ai
                if((clock + Ai) < D[(i-j), ]){   break;  }
            }
            
            if(j < k){  break;  }
            else{   prevAi <- Ai; }  # arrival is rejected due to full system
        }   
        
        NIS <- j + 1            # increment the number of customers in system
        clock <- clock + Ai     # increment clock to next arrival
        W[i, ] <- max(0, W[(i-1), ] + Sr[(i-1)] - Ai)   # Lindley's equation
        TIS[i, ] <- W[i, ] + Sr[i ] # time in system for the i-th customer
        D[i, ] <- clock + TIS[i, ]  # time of i-th departure
        departed <- departed + (j-1)
        prevAi <- 0             # reset rejected arrivals
        
        Arr[i, ] <- clock
    }
    TIS.after.warmup <- TIS[-(1:r0)]
    
    return(list(data = TIS.after.warmup,
                Ybar = mean(TIS.after.warmup), S2 = var(TIS.after.warmup),
                n_Arr = length(All_arr), sum_Arr = sum(All_arr),
                n_Sr = length(Sr), sum_Sr = sum(Sr)))
}

MM1kPdfExp <- function(n_Arr, sum_Arr, n_Sr, sum_Sr, lambda, mu){
    # calculates the sum of the exponents of joint pdf of interarrival time Arr
    # and service times Sr using parameters lambda & mu
    # Inputs (all inputs have the same length):
    # -- n_Arr: vector of the number of interarrival times in each replication
    # -- sum_Arr: vector of SUMs OF interarrival times
    # -- n_Sr: vector of the number of service times in each replication
    # -- sum_Sr : vector of SUMs OF service times
    # -- lambda: contant of interarrival rate
    # -- mu: constant of service rate
    
    # exponent of joint density is considered for numerical stability
    return((n_Sr*log(mu)-mu*sum_Sr) + (n_Arr*log(lambda)-lambda*sum_Arr))
}

MM1kParam_Bootstrap <- function(theta, n_samp, n_out){
    # Bootstrap n_Out parameters, each of which is estimated from n_samp 
    # observations of exponential random variables with rate theta
    # Inputs:
    # -- theta:  constant, rate parameter of exponential distribution
    # -- n_samp: no. of obs., i.e., samples, for one parameter estimation
    # -- n_out:  no. of parameters to be estimated
    # Outputs:
    # -- theta_B: estimated parameters, vector of length n_out
    
    sim_times <- matrix(rexp(n_samp*n_out, theta), n_samp, n_out)
    
    return(theta_out = 1/colMeans(sim_times))
}

MM1kParam_Thres <- function(theta_C, theta_B){
    # convert n_macro param vector theta_C into n_macro param vector whose 
    # entries are the minimum of theta_C and twice the column minimums of theta_B
    # Inputs:
    # -- theta_C: n_macro param vecotr
    # -- theta_B: n_B-by-n_macro parameter matrix
    # Output:
    # -- theta_thres: n_macro param vector whose entries are the minimum of 
    #                 theta_C twice the column minimums of theta_B
    
    n_macro <- length(theta_C) # no. of macro reps
    n_B <- length(theta_B) # no. of macro reps
    thres <- 2*apply(theta_B, 2, min) - 1e-6    # desired thresholds
    
    return(theta_thres = pmin(theta_C, thres))
}

MM1kParam_Resamp <- function(theta_B, n_R){
    # resample from each column the n_B-by-n_macro param matrix theta_B for form 
    # an n_R-by-n_macro parameter matrix
    # Inputs:
    # -- theta_B: n_B-by-n_macro parameter matrix
    # -- n_R:     no. of sample paths per macro rep.
    # Output:
    # -- n_B-by-n_macro param matrix whose columns are random samples of columns
    #    of theta_B
    
    theta_out <- apply(theta_B, MARGIN = 2, function(x) sample(x, replace = TRUE, size = n_R))
    return(theta_out)
}

MM1kParam_ResampThres <- function(theta_B, n_S){
    # resample from each column of the n_B-by-n_macro param matrix theta_B to 
    # form an n_S-by-n_macro parameter matrix, with thresholds that make sure 
    # the likelihood ratio estimator's variance is finite
    # Inputs:
    # -- theta_C: n_macro param vecotr
    # -- theta_B: n_B-by-n_macro parameter matrix
    # -- n_S:     no. of paris of sampling params per macro rep.
    # Output:
    # -- n_S-by-n_macro param matrix whose columns are random samples of columns
    #    of theta_B within the desired thresholds
    
    n_macro <- dim(theta_B)[2] # no. of macro reps
    thres <- 2*apply(theta_B, 2, min) - 1e-6    # desired thresholds
    
    theta_out <- matrix(0, n_S, n_macro)
    for(i in 1:n_macro){
        ind <- theta_B[, i] < thres[i]  # indices for "feasible" params
        samp_i <- theta_B[ind, i]       # feasible params
        theta_out[, i] <- sample(x = samp_i, size = n_S, replace = T)  # sample with replacement
    }
    
    return(theta_out)
}

############ Codes below this point have not been heavily tested

MM1kTIS_SS <- function(lambda, mu, k, n){
    # Simulates the steady-state time in system (TIS) using its limiting
    # distribution, i.e., Erlang(x, mu), where x follows s-s NIS distribution
    # for M/M/1/(k-1)
    
    k <- k-1
    rho <- lambda/mu
    if(rho == 1){
        p0 <- 1/(k+1)
    }
    else{
        p0 <- (1-rho)/(1-rho^(k+1))
    }
    
    # inverse cdf
    temp <- NULL
    u <- runif(n, 0, 1)
    if(rho==1){
        temp <- ceiling((k+1)*u)-1
    }
    else if(rho<1){
        temp <- floor((log(1-((1-rho)/p0)*u))/log(rho))
    }
    else{
        temp <- floor((log(1+((rho-1)/p0)*u))/log(rho))
    }
    
    temp <- temp + 1 # add yourself to # waiting in front of you
    temp <- cbind(temp)
    
    # Erlang distribution is Gamma with the integer shape
    TIS <- rgamma(length(temp), shape = temp, rate = mu)
    
    return(list(data = TIS, Ybar = mean(TIS), S2 = var(TIS)))
}

MM1kTIS_Lindley_SS <-function(lambda, mu, k, r){
    # Simulates the steady-state time in system (TIS) using Lindley's equation
    # for fast initialization, we sample the initial # in system from the truncated geometric 
    # distribution, which is the steady-state distribution of NIS.
    
    # lambda: arrival rate (Poisson process)
    # mu: mean service rate (exponential distribution)
    # k: system capacity
    # r: total observations 
    
    # s-s NIS distribution
    rho = lambda/mu
    ssprob = 0:k
    ssprob = (1-rho)*rho^ssprob/(1-rho^(k+1))
    ssprob_rescaled = ssprob[2:(k+1)]/sum(ssprob[2:(k+1)]) 
    ssprob_rescaled = cumsum(ssprob_rescaled) # we sample from this cdf > excluding the system being empty
    
    # first, we sample the initial NIS at the time of the first arrival (conditional on NIS>=1)
    initial_prob <- runif(1)
    initial_NIS <- min(which(ssprob_rescaled>=initial_prob))-1 # number in system before the first arrival
    
    # placeholders for results, note that only r cumstomers receive services
    Sr <- rexp(r + initial_NIS, mu)       # service times for r customers
    W <- matrix(0,r + initial_NIS,1)      # wait time for r custerms
    TIS <- matrix(0,r + initial_NIS,1)    # time-in-system for r customers
    D <- matrix(0,r + initial_NIS,1)      # departure time for r customers
    
    prevAi <- 0     # system initialization, no previous arrival
    departed <- 0   # system initialization, no previous departure
    
    clock <- 0 # first arrival happens at clock = 0
    
    if(initial_NIS>0)
    {
        D[1, ] <- Sr[1]
        if(initial_NIS>1)
        {
            for(x in 2:initial_NIS)
            {
                D[x, ] = D[x-1, ] + Sr[x]
            }  
        }
        W[initial_NIS+1, ] = D[initial_NIS, ]
        TIS[initial_NIS+1, ] = W[initial_NIS+1, ] + Sr[initial_NIS+1]
        D[initial_NIS+1, ] = TIS[initial_NIS+1, ]
        NIS = initial_NIS + 1
    }  
    else
    {
        # initially the system is empty
        TIS[1, ] = Sr[1]            
        D[1, ] = clock + Sr[1]      
        NIS = 1
    } 
    
    All_arr <- NULL          # all arrivals
    
    for(i in (initial_NIS+2):(r + initial_NIS))
    {
        repeat
        {
            cur_arr <- rexp(1, lambda)       # time-to-next arrival
            All_arr <- c(All_arr, cur_arr)  # time-to-next arrival
            Ai <- prevAi + cur_arr          # arrival time
            if(clock+Ai > D[(i-1), ])  # this is the case where the arrival finds the system empty.
            {  j=0; break;}
            
            for(j in NIS:1)
            {
                # count how many leave in time period of Ai
                if((clock + Ai) < D[(i-j), ]){   break;  }
            }
            
            if(j < k){  break;  }
            else{   prevAi <- Ai;   }
        }
        
        NIS <- j+1
        clock <- clock + Ai
        W[i, ] <- max(0, W[(i-1), ]+Sr[(i-1)]-Ai)
        TIS[i, ] <- W[i, ] + Sr[i]
        D[i, ] <- clock + TIS[i, ]
        departed <- departed + (j-1)
        prevAi <- 0
    }
    data = TIS[(initial_NIS+1):(initial_NIS + r)]
    
    return(list(data = data,
                Ybar = mean(data),
                S2 = var(data),
                Sr = Sr, All_arr = All_arr))
}