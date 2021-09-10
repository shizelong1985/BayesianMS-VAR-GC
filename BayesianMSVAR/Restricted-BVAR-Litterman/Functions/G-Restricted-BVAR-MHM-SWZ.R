
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



# Modified harmonic Mean for restricted MS-VAR models
# SWZ (2008)
                   
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


G.MHM.BVAR.SWZ  <- function(Gibbs.output, priors, restrictions, debug=FALSE, print.iterations=100){

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
    log.likelihood  <- array(NA, N)
    # Priors
    log.prior.S     <- array(NA, c(K,N))
    log.prior.R     <- array(NA, N)
    log.prior.Beta  <- array(NA, N)
    log.prior.joint <- array(NA, N)
    # Truncated normal
    log.truncated.normal    <- array(NA, N)


#---------------------------------------------------------------------------------------------------
# Creation of the posteriors vector
#---------------------------------------------------------------------------------------------------
    output.matrix       <- BVAR.posteriors.matrix(Gibbs.output, restrictions)
    posteriors.vector   <- output.matrix$posteriors.vector
    number.parameters   <- output.matrix$number.parameters
        
    ### # Truncation parameters as quantiles: alpha percentiles set as bounds
    ### truncation.interval <- array(NA, c(number.parameters,2))               
    ### for(i in 1:number.parameters){
        ### truncation.interval[i,] <- quantile(posteriors.vector[i,], probs=c(alpha/2,1-alpha/2))
    ### }

    ### # Mean and variance of the posteriors vector
    ### posterior.mean  <- apply(posteriors.vector, 1, mean)
    ### posterior.var   <- cov(t(posteriors.vector))
                  
   
#---------------------------------------------------------------------------------------------------
# Loop over all iterations of the Gibbs sampler
#---------------------------------------------------------------------------------------------------


    #-----------------------------------------------------------------------------------------------
    # 1/4: Likelihood 
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n1/4: Likelihood, loop over all iterations of the Gibbs sampler")
    for(iteration in 1:N){

        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        # Residuals calculated from Beta
        aux$Beta    <- matrix(posteriors$Beta[,iteration], ncol=K)
        res.tmp     <- Y - X %*% aux$Beta
        aux$U       <- matrix(res.tmp, ncol=K)

        # Densities
        eta.tmp         <- dmvnorm(matrix(aux$U,ncol=K), sigma=matrix(aux$Sigma,ncol=K)) 
        eta.t  <- matrix(eta.tmp, nrow=1)

        # Log likelihood 
        log.likelihood[iteration]   <- sum( log(eta.t ))
    }                                                                                    



    #-----------------------------------------------------------------------------------------------
    # 2/4: Priors
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n2/4: Priors, loop over all iterations of the Gibbs sampler")
    
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
    log.prior.S.joint   <- colSums(log.prior.S)


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
    log.prior.Beta.joint    <-  log.prior.Beta


    # Joint priors
    #-----------------------------------------------------------------------------
    for(iteration in 1:N){
        log.prior.joint[iteration]  <- log.prior.S.joint[iteration] + log.prior.R.joint[iteration] + log.prior.Beta.joint[iteration] 
    }


    #-----------------------------------------------------------------------------------------------
    # 3/4: Truncated elliptical density
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n3/4: Truncated elliptical density, loop over all iterations of the Gibbs sampler")

    # Posterior mode
    Theta.hat   <- array(NA, c(number.parameters,1))
    for(i in 1:number.parameters){
        Theta.hat[i,] <- quantile(posteriors.vector[i,], probs=.5)
    }           
    # Scaling matrix
    minus       <- posteriors.vector - Theta.hat %*% matrix(1, 1, N)
    Sigma.hat   <- array(0, c(number.parameters,number.parameters))
    for(iteration in 1:N){
        Sigma.hat   <- Sigma.hat + tcrossprod(as.matrix(minus[,iteration]))
    }
    Sigma.hat   <- Sigma.hat / N
    # Square root of Sigma.hat (http://realizationsinbiostatistics.blogspot.de/2008/08/matrix-square-roots-in-r_18.html)
    S.eig   <- eigen(Sigma.hat)
    S.hat   <- S.eig$vectors %*% diag(sqrt(S.eig$values)) %*% solve(S.eig$vectors)
    # Loop for r.i
    r.i     <- array(NA, N)
    Sigma.hat.inv   <- solve(Sigma.hat)
    for(iteration in 1:N){
        r.i[iteration]  <- sqrt( t(as.matrix(minus[,iteration])) %*% Sigma.hat.inv %*% as.matrix(minus[,iteration]))
    }
    # Hyperparameters for f
    c.1 <- quantile(r.i, probs=0.01, names=FALSE)
    c.10<- quantile(r.i, probs=0.1, names=FALSE)
    c.90<- quantile(r.i, probs=0.9, names=FALSE)
    a   <- c.1
    v   <- log(1/9)/log(c.10/c.90)
    b   <- c.90/0.9^(1/v)

    # Step 1: kernel, L (lower bound truncation on the posterior density kernel itself)
    kernel  <- log.likelihood + log.prior.joint
    ### kmin    <- min(kernel)
    ### kmax    <- max(kernel)
    L       <- quantile(kernel, probs=.1, names=FALSE) # a good choice of L is a value such that 90% of draws from the posterior distribution lie in theta_L.

    # Step 2: densities from elliptical distribution
    g.Theta     <- array(NA, N)
    k           <- number.parameters
    constant    <- gamma(k/2) / (2*pi^(k/2) * abs(det(S.hat)))
    for(iteration in 1:N){
        f.r                 <- v * r.i[iteration]^(v-1) / ( b^v - a^v)
        g.Theta[iteration]  <- constant * f.r / r.i[iteration]^(k-1)
    }
    q.L.hat <- 0.9


    #-----------------------------------------------------------------------------------------------
    # 4/4: Modied Harmonic Mean
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n4/4: Computation of the MHM")
    # Equation (15), Sims, Wagoner, Zha (2008)
    iteration.product   <- array(NA, N)
    for(iteration in 1:N){
        iteration.truncated             <- log.truncated.normal[iteration]
        iteration.likelihood            <- log.likelihood[iteration]
        iteration.priors                <- log.prior.S.joint[iteration] + log.prior.R.joint[iteration] + log.prior.Beta.joint[iteration]
        iteration.product[iteration]    <- iteration.truncated - iteration.likelihood - iteration.priors
    }
    
    # Work with logs only
    # 1/3: Summation using formula http://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation.2Fsubtraction
    finite.elements <- iteration.product[which(is.finite(iteration.product))] # Finite elements only

    acceptance.rate <- length(finite.elements)/length(iteration.product) * 100

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
