
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



# Gibbs sampler for restricted VAR models
                   
# Arguments:
#-----------------------------------------------------------------------------
# N: number of iterations for the Gibbs sampler


# h: number of periods to forecast


# priors:
# List containing:
    # Beta.bar (K*(1+K*p)), 1) - priors for the VAR coeffs, mean
    # V.Beta.bar.m1 (K*(1+K*p), K*(1+K*p)) - priors for the VAR coeffs, variance
    # S.bar (K, K) - priors for the cov. matrix, scale 
    # nu.bar (positive number) - priors for the cov. matrix, degrees of freedom


# restrictions
# List containing: 
    # Beta ((1+K*p), K) - Restrictions for AR coefficients
    # Sigma (K, K) - Restrictions for covariance matrix
# About restrictions: for each parameter, enter 0 if the parameter is not to be estimated, or enter 1 otherwise.


# starting.values: for startup of the algorithm
# List containing:
    # Y (T,K) - data
    # U (TT,K), with TT=T-p - residuals
    # Beta ((K*p+1),K) - AR coefficients
    # Sigma (K,K) - variance covariance matrices


# debug: debug mode, print informations


# print.iterations: number for how often to print the results of the iteration 


# Returns:
#-----------------------------------------------------------------------------
# last.draws: last draws from the Gibbs sampler, can be used to launch iterations for a new sampler
# List containing:
    # Y - data
    # U - residuals
    # Beta - draws for intercepts and AR coefficients
    # Sigma - draws for variance


# posteriors: All draws from all iterations
# List containing:
    # Beta ((K*p+1)*K,N)
    # Sigma (K*K,N)
    # F.Y (h,K,N) - forecasts of series, horizon h 


  
G.BVAR.Griddy <- function(N, h, priors, starting.values, restrictions, debug=FALSE, print.iterations=100) {

    # parameters
    #-----------------------------------------------------------------------------
    aux <- starting.values
    p   <- dim(aux$Y)[1] - dim(aux$U)[1]
    TT  <- dim(aux$Y)[1] - p
    K   <- dim(aux$Y)[2] 


    # List of vectorized posteriors for each iteration
    #----------------------------------------------------------------------------- 
    posteriors <- list(
        Beta    = array(NA, c((K*p+1)*K,N)),
        Sigma   = array(NA, c(K*K,N)),
        S       = array(NA, c(K,N)),
        R       = array(NA, c(K*(K-1)/2,N)),  # only upper-triangular elements
        F.Y     = array(NA, c(h,K,N))
        )

    # N iterations of the Gibbs
    #----------------------------------------------------------------------------- 
    for(iteration in 1:N){

        # Standard deviations step
        if(debug) print("Standard deviation step")
        aux <- G.standard.deviation(aux, priors, restrictions)  

        # Correlations step
        if(debug) print("Correlation step")
        aux <- G.correlation(aux, priors, restrictions)    

        # Regression step 
        if(debug) print("Regression step")
        aux <- G.regression.BVAR.INWP(aux, priors, restrictions) 

        # Forecasting 
        if(debug) print("Forecasting")
        posteriors  <- G.forecast.BVAR(aux, posteriors, h, iteration)

        # Save posteriors as vectors
        posteriors$Beta[,iteration] <- matrix(aux$Beta[,], ncol=1) 
        posteriors$Sigma[,iteration]<- matrix(aux$Sigma[,], ncol=1)  
        posteriors$S[,iteration]<- matrix(aux$S[], ncol=1)  
        posteriors$R[,iteration]<- matrix(aux$R[upper.tri(aux$R[,],diag=FALSE)],ncol=1)  

        # Print iteration results
        if((debug) && ((iteration %% print.iterations)==0)){
            cat("---------------------------------------------------------- \nIteration:",iteration,"/", N,"\n")  
            cat("Intercepts + AR coefficients\n")
            print(aux$Beta)
            cat("Sigma\n")
            print(aux$Sigma)
        }else if((iteration %% print.iterations)==0) cat(" ",iteration)  
    }


    # Output
    #-----------------------------------------------------------------------------
    output  <- list(
                last.draws  = aux,
                posteriors  = posteriors
               )

    return(output)
                    
}



#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------




G.BVAR.Adaptive.Grid <- function(N, h, priors, starting.values, G.grid, restrictions, debug=FALSE, print.iterations=100) {

    # parameters
    #-----------------------------------------------------------------------------
    aux <- starting.values
    p   <- dim(aux$Y)[1] - dim(aux$U)[1]
    TT  <- dim(aux$Y)[1] - p
    K   <- dim(aux$Y)[2]


    # List of vectorized posteriors for each iteration
    #----------------------------------------------------------------------------- 
    posteriors <- list(
        Beta    = array(NA, c((K*p+1)*K,N)),
        Sigma   = array(NA, c(K*K,N)),
        S       = array(NA, c(K,N)),
        R       = array(NA, c(K*(K-1)/2,N)),  # only upper-triangular elements 
        F.Y     = array(NA, c(h,K,N))
        )

    # N iterations of the Gibbs
    #----------------------------------------------------------------------------- 
    for(iteration in 1:N){

        # Standard deviations step
        if(debug) print("Standard deviation step")
        aux <- G.standard.deviation.adaptive.grid(aux, priors, restrictions, G.grid)  

        # Correlations step
        if(debug) print("Correlation step")
        aux <- G.correlation.adaptive.grid(aux, priors, restrictions, G.grid) 

        # Regression step 
        if(debug) print("Regression step")
        aux <- G.regression.BVAR.INWP(aux, priors, restrictions) 

        # Forecasting 
        if(debug) print("Forecasting")
        posteriors  <- G.forecast.BVAR(aux, posteriors, h, iteration)

        # Save posteriors as vectors
        posteriors$Beta[,iteration] <- matrix(aux$Beta, ncol=1) 
        posteriors$Sigma[,iteration]<- matrix(aux$Sigma, ncol=1)  
        posteriors$S[,iteration]<- matrix(aux$S, ncol=1)  
        posteriors$R[,iteration]<- matrix(aux$R[upper.tri(aux$R[,],diag=FALSE)],ncol=1)  

        # Print iteration results
        if((debug) && ((iteration %% print.iterations)==0)){
            cat("---------------------------------------------------------- \nIteration:",iteration,"/", N,"\n")  
            cat("Intercepts + AR coefficients\n")
            print(aux$Beta)
            cat("Sigma\n")
            print(aux$Sigma)
        }else if((iteration %% print.iterations)==0) cat(" ",iteration)  
    }


    # Output
    #-----------------------------------------------------------------------------
    output  <- list(
                last.draws  = aux,
                posteriors  = posteriors
               )

    return(output)
                    
}               
