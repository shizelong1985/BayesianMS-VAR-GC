
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



# Modified harmonic Mean for restricted MS-VAR models
# Geweke (1999)
                   
# Arguments:
#-----------------------------------------------------------------------------
# Gibbs.output: output from G.BVAR()
 
# priors:
# see header of G-Restricted-BVAR.R

# restrictions:
# see header of G-Restricted-BVAR.R

# alpha: truncation for truncated normal distribution

# debug: debug mode, print informations

# print.iterations: number for how often to print the results of the iteration


# Returns:
#-----------------------------------------------------------------------------
# MHM
# log.prior
# log.likelihood
# log.truncated.normal


G.MHM.BVAR.Geweke   <- function(Gibbs.output, priors, restrictions, alpha, debug=FALSE, print.iterations=100){

    #----------------------------------------------------------------------------- 
    # Setup constants 
    #-----------------------------------------------------------------------------
    aux         <- Gibbs.output$last.draws
    posteriors  <- Gibbs.output$posteriors
    p           <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T           <- dim(aux$Y)[1]
    TT          <- T - p
    K           <- dim(aux$Y)[2] 
    N           <- dim(posteriors$Beta)[2]

    # Y and X matrices
    X   <- matrix(data=1, nrow=TT, ncol=1)	# A column of ones for the intercept
    if(p > 0){
        for(lag in 1:p){
            Y.m.j   <- aux$Y[(p+1-lag):(T-lag),]
            X       <- cbind(X,Y.m.j)
        }
    }
    Y <- matrix(aux$Y[(p+1):T,],ncol=K)  



    #-----------------------------------------------------------------------------
    # Storage

    #-----------------------------------------------------------------------------
    # Likelihood
    log.likelihood  <- array(NA, c(N))
    # Priors
    log.prior.S     <- array(NA, c(K,N))
    log.prior.R     <- array(NA, N)
    log.prior.Beta  <- array(NA, N)
    # Truncated normal
    log.truncated.normal    <- array(NA, c(N))


#---------------------------------------------------------------------------------------------------
# Determination of the truncation set
#---------------------------------------------------------------------------------------------------
    output.matrix       <- BVAR.posteriors.matrix(Gibbs.output, restrictions)
    posteriors.vector   <- output.matrix$posteriors.vector
    number.parameters   <- output.matrix$number.parameters

	########################################
	# Changed by TW
	########################################
    # Cancelled:
    # Truncation parameters as quantiles: alpha percentiles set as bounds
#    truncation.interval <- array(NA, c(number.parameters,2))               
#    for(i in 1:number.parameters){
#        truncation.interval[i,] <- quantile(posteriors.vector[i,], probs=c(alpha/2,1-alpha/2))
#    }

    # Mean and variance of the posteriors vector
    posterior.mean  <- apply(posteriors.vector, 1, mean)
    posterior.var   <- cov(t(posteriors.vector))


	# Added:
	quadratic.form	= array(NA, N)
	for(iteration in 1:N){
		        
        #-------------------------------------------------------------------------
        par.x       <- as.vector(posteriors.vector[,iteration])
        quadratic.form[iteration]	= t(par.x - posterior.mean) %*% solve(posterior.var) %*% (par.x - posterior.mean)
	}
	
	# A vector of indices for the posterior draws for which the truncation condition holds [TW]
    indices 	= which(quadratic.form <= qchisq(1-alpha, df = number.parameters) )    
    N.truncated	= length(indices)

        
#---------------------------------------------------------------------------------------------------
# Loop over not dismissed iterations of the Gibbs sampler
#---------------------------------------------------------------------------------------------------


    #-----------------------------------------------------------------------------------------------
    # 1/4: Likelihood 
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n1/4: Likelihood, loop over not dismissed iterations of the Gibbs sampler")
    for(iteration in 1:N){

        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        # Residuals calculated from Beta
        aux$Beta     <- matrix(posteriors$Beta[,iteration], ncol=K)
        aux$Sigma    <- matrix(posteriors$Sigma[,iteration], ncol=K)
        res.tmp      <- Y - X %*% aux$Beta
        aux$U        <- matrix(res.tmp, ncol=K)

        # Densities
        eta.tmp         <- dmvnorm(matrix(aux$U,ncol=K), sigma=matrix(aux$Sigma,ncol=K)) 
        eta.t  <- matrix(eta.tmp, nrow=1)

        # Log likelihood 
        log.likelihood[iteration]   <- sum( log(eta.t ))
    }                                                                                    

# 	log.likelihood	= log.likelihood[indices]

    #-----------------------------------------------------------------------------------------------
    # 2/4: Priors
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n2/4: Priors, loop over not dismissed iterations of the Gibbs sampler")
    
    # Standard deviations
    #-----------------------------------------------------------------------------
    if(debug) cat("\nPriors, standard deviations\n")
    for(iteration in 1:N){
        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        S.iteration <- posteriors$S[,iteration]
        for(j in 1:K){
             # Prior for Si, log normal
            log.prior.S[j,iteration] <- dlnorm(S.iteration[j], meanlog=priors$S.mean, sdlog=priors$S.sd, log=TRUE) 
        }
    }      
    # Sum for all iterations
#     log.prior.S			= log.prior.S[,indices]
    log.prior.S.joint   = colSums(log.prior.S)


    # Correlations
    #-----------------------------------------------------------------------------
    if(debug) cat("\nPriors, correlations\n")
    for(iteration in 1:N){
        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        # Homoskedasticity
        R   <-  diag(K)
        R[upper.tri(R, diag=FALSE)] <- posteriors$R[,iteration]
        R[lower.tri(R, diag=FALSE)] <- posteriors$R[,iteration]
        
        joint.density.R <- det(R)^((K*(K-1)/2)-1)
        for(i in 1:K){
            # Principal submatrices
            joint.density.R <- joint.density.R * det(matrix(R[-i,-i],nrow=K-1))^(-(K+1)/2)
        }
        log.prior.R[iteration]    <- log(joint.density.R)

       
    }      
    # Sum for all iterations
#     log.prior.R			= log.prior.R[indices]
    log.prior.R.joint   <- log.prior.R


    # Intercepts and autoregressive coefficients
    #-----------------------------------------------------------------------------
    if(debug) cat("\nPriors, intercepts and autoregressive coefficients\n") 
    for(iteration in 1:N){
        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        par.x       <- posteriors$Beta[,iteration]
        par.mean    <- priors$Beta.bar
        par.sigma   <- solve(priors$V.Beta.bar.m1)
        ldensity    <- dmvnorm(x=par.x, mean=par.mean, sigma=par.sigma, log=TRUE)

        log.prior.Beta[iteration]   <- ldensity
    }      
    # Sum for all iterations
#     log.prior.Beta			= log.prior.Beta[indices]
    log.prior.Beta.joint    <-  log.prior.Beta


    #-----------------------------------------------------------------------------------------------
    # 3/4: Truncated multivariate normal
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n3/4: Truncated multivariate normal, loop over not dismissed iterations of the Gibbs sampler")
    
#     for(iteration in 1:indices){
#         if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 
#         
#         # Density from truncated multivariate normal
#         #-------------------------------------------------------------------------
#         par.x       <- as.vector(posteriors.vector[,iteration])
#         #par.mean    <- as.vector(posterior.mean)
#         #par.sigma   <- posterior.var
# 
#         #par.lower   <- truncation.interval[,1]
#         #par.upper   <- truncation.interval[,2]
#         #log.truncated.normal[iteration] <- dtmvnorm(x=par.x, mean=par.mean, sigma=par.sigma, lower=par.lower, upper=par.upper,log=TRUE)
#         
#         log.truncated.normal[iteration] = -log(1-alpha) - (number.parameters/2)*log(2*pi) - 0.5*log(det(posterior.var)) - 0.5*quadratic.form[iteration]
#         
#     }
    par.mean    <- as.vector(posterior.mean)
    par.sigma   <- posterior.var
    log.truncated.normal = apply(posteriors.vector,2,dmvnorm,mean=par.mean, sigma=par.sigma,log=TRUE) -log(1-alpha)
    log.truncated.normal[-(indices)] = -Inf

# 	log.truncated.normal = log.truncated.normal[indices]

    #-----------------------------------------------------------------------------------------------
    # 4/4: Modied Harmonic Mean
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n4/4: Computation of the MHM")
    # Equation (15), Sims, Wagoner, Zha (2008)
    priors.joint        <- log.prior.S.joint + log.prior.R.joint + log.prior.Beta.joint
    iterations.product  <- log.truncated.normal - log.likelihood - priors.joint
    
#     iteration.product   <- array(NA, N.truncated)
#     for(iteration in 1:N.truncated){
#         iteration.truncated             <- log.truncated.normal[iteration]
#         iteration.likelihood            <- log.likelihood[iteration]
#         iteration.priors                <- log.prior.S.joint[iteration] + log.prior.R.joint[iteration] + log.prior.Beta.joint[iteration]
#         iteration.product[iteration]    <- iteration.truncated - iteration.likelihood - iteration.priors
#     }
    
    # Work with logs only
    # 1/3: Summation using formula http://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation.2Fsubtraction
    finite.elements <- iterations.product[which(is.finite(iterations.product))] # Finite elements only

    acceptance.rate <- length(finite.elements)/N * 100

    finite.elements <- sort(finite.elements, decreasing=TRUE) # Sorting by descending order allows to avoid infinite exponents due to the big difference         
    N1              <- length(finite.elements)
    tmp.sum <- 0
    for(iteration in 2:N1){
        tmp.sum <- tmp.sum + exp(finite.elements[iteration] - finite.elements[1])
    }
    log.sum <- finite.elements[1] + log(1+tmp.sum)
    # 2/3: divide by N, number Gibbs iterations
    log.mean<- log.sum - log(N)
    # 3/3: invert, yields MHM
    MHM     <- - log.mean

    #-----------------------------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------------------------
    ### log.prior   <- list(log.prior.PR_TR = log.prior.PR_TR.joint,
                        ### log.prior.S     = log.prior.S.joint,
                        ### log.prior.R     = log.prior.R.joint,
                        ### log.prior.Beta  = log.prior.Beta.joint
                        ### )

    ### output  <- list(MHM                 = MHM,
                    ### log.prior           = log.prior,
                    ### log.likelihood      = log.likelihood,
                    ### log.truncated.normal= log.truncated.normal
                    ### )

    cat("\n-----------------------------\nMHM:", MHM, ", alpha:", alpha, ", acceptance rate:", acceptance.rate)

    output  <- list(    MHM = MHM,
                        acceptance.rate = acceptance.rate,
                        alpha = alpha   )

    return(output)

}



    



#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------




BVAR.posteriors.matrix <- function(Gibbs.output, restrictions)
{
    #----------------------------------------------------------------------------- 
    # Setup constants 
    #-----------------------------------------------------------------------------
    aux         <- Gibbs.output$last.draws
    posteriors  <- Gibbs.output$posteriors
    M           <- dim(aux$PR_TR)[1]
    p           <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T           <- dim(aux$Y)[1]
    TT          <- T - p
    K           <- dim(aux$Y)[2] 
    N           <- dim(posteriors$Beta)[2]
 
    # Vectorize all the posteriors for all Gibbs iterations: number of parameters * N
    # S
    number.parameters   <-sum(diag(restrictions$Sigma))
    cat("S:", number.parameters, "\n")
    # R
    number.parameters   <- number.parameters + sum(restrictions$Sigma[upper.tri(restrictions$Sigma, diag=FALSE)])
    cat("R:", number.parameters, "\n")
    # Beta
    number.parameters   <- number.parameters + sum(restrictions$Beta) 
    cat("Beta:", number.parameters, "\n")


    # Posterior vectors
    posteriors.vector   <- array(NA, c(number.parameters,N))

    for(iteration in 1:N){
        # S
        s.restriction   <- diag(restrictions$Sigma)
        s.iteration     <- posteriors$S[,iteration]
        iteration.vector<- s.iteration[which(s.restriction != 0)]

        # R
        r.restriction   <- restrictions$Sigma[upper.tri(restrictions$Sigma, diag=FALSE)]
        r.iteration     <- posteriors$R[,iteration]
        iteration.vector    <- append(iteration.vector, r.iteration[which(r.restriction != 0)] ) 
        
        # Beta
        beta        <- matrix(posteriors$Beta[,iteration],ncol=K)
        b.iteration <-beta[which(restrictions$Beta != 0)]
        iteration.vector    <- append(iteration.vector, b.iteration)   


        # Store iteration vector into posteriors matrix
        posteriors.vector[,iteration]   <- iteration.vector    
    }

    #-----------------------------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------------------------
    output  <- list(    posteriors.vector = posteriors.vector,
                        number.parameters = number.parameters   )

    return(output)      

}           
