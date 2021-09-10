
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



# Gibbs sampler for restricted MS-VAR models
                   
# Arguments:
#-----------------------------------------------------------------------------
# N: number of iterations for the Gibbs sampler


# h: number of periods to forecast


# priors:
# List containing:
    # alpha (M,M) - priors for the transition matrix, see the function G.hidden.markov.dirichlet.priors()
    # Beta.bar (K*(1+K*p), 1) - priors for the VAR coeffs, mean
    # V.Beta.bar.m1 (K*(1+K*p), K*(1+K*p)) - priors for the VAR coeffs, variance
    # S.mean (number) - mean of lognorm standard deviation coeffs
    # S.sd (positive number) - SD of lognorm standard deviation coeffs
    # w (d, 1) - Hidden Markov Chain, transition probabilities

# restrictions
# List containing: 
    # Beta.0 ((1+K*p), K) - Restrictions for AR coefficients - regime independent
    # Beta.St ((1+K*p), K, M) - Restrictions AR coefficients - regime dependent
    # heteroskedasticity: TRUE or FALSE
    # Sigma.0 (K, K) - Restrictions for Sigma coefficients - regime independent
    # Sigm.St (K, K, M) - Restrictions Sigma coefficients - regime dependent   
    # M (M, M) - Restrictions for transition probabilities 
    # dj (M, 1) - Restrictions for transition probabilities
# About restrictions: for Beta.0, Beta.St, Sigma.0, Sigma.St, enter 0 if the parameter is not to be estimated, or enter 1 otherwise.


# starting.values: for startup of the algorithm
# List containing:
    # Y (T, K) - data
    # U (TT, K, M) - residuals
    # xi (M, TT) - estimated regimes, not used for initialization, i.e. a (M*TT) matrix containing NAs will do
    # Beta.0 ((1+K*p), K) - AR coefficients - regime independent
    # Beta.St ((1+K*p), K, M) - AR coefficients - regime dependent
    # Sigma.0 (K, K) - variance covariance matrices - regime independent
    # Sigma.St (K, K, M) - variance covariance matrices - regime dependent  
    # PR_TR (M, M) - transition probabilities
	# w (d,1) - restricted transition probabilities


# debug: debug mode, print informations


# print.iterations: number for how often to print the results of the iteration


# Returns:
#-----------------------------------------------------------------------------
# last.draws: last draws from the Gibbs sampler, can be used to launch iterations for a new sampler
# List containing:
    # Y - data
    # u - residuals
    # xi (M, TT) - draws for regime
    # Beta - draws for intercepts and AR coefficients
    # Sigma - draws for variance
    # S - draws for standard deviations
    # R - draws for correlations
    # PR_TR - draws for transition probabilities
	# w - draws for transition probabilities


# posteriors: All draws from all iterations
# List containing:
    # Beta ((1+K*p)*K,M,N)
    # Sigma (K*K,M,N)
    # S (K,M,N)
    # R (K*(K-1)/2,M,N)
    # PR_TR (M*M,N)
	# w (d,N)
    # S.t (TT,N) - draws of regimes
    # F.S.t (h,N) - forecasts of regimes, horizon h
    # F.Y (h,K,N) - forecasts of series, horizon h                 


  
G.MSBVAR.griddy <- function(N, h, priors, restrictions, starting.values, debug=FALSE, print.iterations=50) {

    # parameters
    #-----------------------------------------------------------------------------
    aux <- starting.values
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]
    d   <- sum(restrictions$dj) 

    # List of vectorized posteriors for each iteration
    #----------------------------------------------------------------------------- 
    posteriors <- list(
        Beta    = array(NA, c(K*(1+K*p),M,N)),
        Sigma   = array(NA, c(K*K,M,N)),
        S       = array(NA, c(K,M,N)),
        R       = array(NA, c(K*(K-1)/2,M,N)),  # only upper-triangular elements 
        PR_TR   = array(NA, c(M*M,N)),
        w   	= array(NA, c(d,N)),
        S.t     = array(NA, c(TT,N)),
        F.S.t   = array(NA, c(h,N)),
        F.Y     = array(NA, c(h,K,N))
        )
        
    # N iterations of the Gibbs
    #----------------------------------------------------------------------------- 
    for(iteration in 1:N){
        
        if((iteration %% 100)==0) cat(" ",iteration)

        # Filtering and smoothing step
        if(debug) print("Filtering and smoothing step")
        aux <- G.filtering.smoothing(aux)

        # Hidden Markov Chain step
        if(debug) print("Hidden Markov Chain step") 
        aux <- G.hidden.markov.chain(aux, priors, restrictions)

        # Standard deviations step
        if(debug) print("Standard deviation step")
        tmp <- try(G.standard.deviation(aux, priors, restrictions))
		if(!inherits(tmp, "try-error")){
	    	aux <- tmp
    	}else{
			 cat(iteration,"Standard deviation error\n")
		}

        # Correlations step
        if(debug) print("Correlation step")
        aux <- G.correlation(aux, priors, restrictions) 

        # Regression step
        if(debug) print("Regression step")
        aux <- G.regression.Beta.INWP(aux, priors, restrictions) 

        # Forecasting 
        if(debug) print("Forecasting step")
        posteriors  <- G.forecast.MSVAR(aux, posteriors, h, iteration)

        # Save posteriors as vectors
        for(regime in 1:M){
            posteriors$Beta[,regime,iteration]  <- matrix(aux$Beta[,,regime], ncol=1)
            posteriors$S[,regime,iteration]     <- matrix(aux$S[,regime], ncol=1)
            R   <- aux$R[,,regime]
            posteriors$R[,regime,iteration]     <- matrix(R[upper.tri(R,diag=FALSE)],ncol=1) 
            posteriors$Sigma[,regime,iteration] <- matrix(aux$Sigma[,,regime], ncol=1)    
        }
        posteriors$PR_TR[,iteration]    <- matrix(aux$PR_TR, ncol=1)
        posteriors$w[,iteration]    	<- matrix(aux$w, ncol=1)
        posteriors$S.t[,iteration]      <- matrix(max.col(t(aux$xi)), ncol=1)

        # Print iteration results
        if((iteration %% print.iterations)==0){
            cat("\n---------------------------------------------------------- \nIteration:",iteration,"/", N,"\n") 
            cat("Count transitions\n")
            print(count.regime.transitions(aux$xi))
            cat("Transition probabilities\n")
            print(aux$PR_TR)
            ### cat("S\n") 
            ### print(aux$S)
            ### cat("R\n") 
            ### print(aux$R)
            cat("Sigma\n")
            print(aux$Sigma)
            cat("Intercepts + AR coefficients\n")
            print(aux$Beta)
        }
    }


    # Output
    #-----------------------------------------------------------------------------
    output  <- list(
                last.draws  = aux,
                posteriors  = posteriors
               )

    return(output)
                    
}
