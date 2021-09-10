
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



# Modified harmonic Mean for restricted MS-VAR models
                   
# Arguments:
#-----------------------------------------------------------------------------
# Gibbs.output: output from G.MSBVAR()
 
# priors:
# see header of G-Restricted-MSBVAR.R

# restrictions:
# see header of G-Restricted-MSBVAR.R

# alpha: truncation for truncated normal distribution 
 
# debug: debug mode, print informations

# print.iterations: number for how often to print the results of the iteration


# Returns:
#-----------------------------------------------------------------------------
# MHM
# log.prior
# log.likelihood
# log.truncated.normal


G.MHM.MSVAR <- function(Gibbs.output, priors, restrictions, alpha, debug=FALSE, print.iterations=100){

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
    N           <- dim(posteriors$Beta)[3]
    Q           <- length(restrictions$dj)

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
    log.prior.PR_TR <- array(NA, c(Q,N))
    if(sum(restrictions$Sigma.St) == 0){
        # Homoskedasticity
        log.prior.S <- array(NA, c(K,1,N))
        log.prior.R <- array(NA, c(1,N))
    }else{
        # Heteroskedasticity
        log.prior.S <- array(NA, c(K,M,N))
        log.prior.R <- array(NA, c(M,N))
    }
    log.prior.Beta  <- array(NA, c(M,N))
    # Truncated normal
    log.truncated.normal    <- array(NA, c(N))


#---------------------------------------------------------------------------------------------------
# Determination of the truncation set
#---------------------------------------------------------------------------------------------------
    output.matrix       <- MSVAR.posteriors.matrix(Gibbs.output, restrictions)
    posteriors.vector   <- output.matrix$posteriors.vector
    number.parameters   <- output.matrix$number.parameters
        
    # Truncation parameters as quantiles: alpha percentiles set as bounds
    truncation.interval <- array(NA, c(number.parameters,2))               
    for(i in 1:number.parameters){
        truncation.interval[i,] <- quantile(posteriors.vector[i,], probs=c(alpha/2,1-alpha/2))
    }

    # Mean and variance of the posteriors vector
    posterior.mean  <- apply(posteriors.vector, 1, mean)
    posterior.var   <- cov(t(posteriors.vector))
                  
   
#---------------------------------------------------------------------------------------------------
# Loop over all iterations of the Gibbs sampler
#---------------------------------------------------------------------------------------------------


    #-----------------------------------------------------------------------------------------------
    # 1/4: Likelihood 
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n1/4: Likelihood, loop over all iterations of the Gibbs sampler")
    for(iteration in 1:N){

        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        # Residuals calculated from Beta.star: depending only on Beta.star, they can be computed now
        for(regime in 1:M){
            # Residuals: from Beta star, used for additional Gibbs runs
            aux$Beta[,,regime]  <- matrix(posteriors$Beta[,regime,iteration], ncol=K)
            res.tmp             <- Y - X %*% aux$Beta[,,regime]
            aux$U[,,regime]     <- matrix(res.tmp, ncol=K)
        }   

        # Densities for each regimes
        eta.t 	<- matrix(NA, M, TT)
        for(regime in 1:M){
            # Densities:
            eta.tmp         <- dmvnorm(matrix(aux$U[,,regime],ncol=K), sigma=matrix(aux$Sigma[,,regime],ncol=K)) 
            eta.t[regime,]  <- matrix(eta.tmp, nrow=1)
        }

        # Rebuild drawn regime history
        for(regime in 1:M){
            aux$xi[regime,] <- 0
            aux$xi[regime,which(posteriors$S.t[,iteration]==regime)] <- 1
        }

        # Log likelihood 
        log.likelihood[iteration]   <- sum( log(colSums(eta.t * aux$xi)))
    }                                                                                    



    #-----------------------------------------------------------------------------------------------
    # 2/4: Priors
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n2/4: Priors, loop over all iterations of the Gibbs sampler")
    
    # Transition probabilities
    #-----------------------------------------------------------------------------
    if(debug) cat("\nPriors, transition probabilities\n")
    for(iteration in 1:N){
        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        ### PR_TR.vec           <- matrix(matrix(posteriors$PR_TR[,iteration], nrow=M,byrow=TRUE),ncol=1) 
        ### w                   <- solve(restrictions$M) %*% PR_TR.vec 
        w                   <- posteriors$w[,iteration]

        for(j in 1:Q){
            index.begin <- sum(restrictions$dj[1:j-1]) + 1
            dj          <- restrictions$dj[j]
            index.end   <- index.begin + dj - 1

            log.prior.PR_TR[j,iteration]    <- log(ddirichlet(x=w[index.begin:index.end], alpha=priors$w[index.begin:index.end]))
        } 
    }
    # Sum for all iterations
    log.prior.PR_TR.joint   <- colSums(log.prior.PR_TR)


    # Standard deviations
    #-----------------------------------------------------------------------------
    if(debug) cat("\nPriors, standard deviations\n")
    for(iteration in 1:N){
        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        if(sum(restrictions$Sigma.St) == 0){
            # Homoskedasticity
            S.iteration <- posteriors$S[,1,iteration]
            for(j in 1:K){
                 # Prior for Si, log normal
                log.prior.S[j,1,iteration] <- dlnorm(S.iteration[j], meanlog=priors$S.mean, sdlog=priors$S.sd, log=TRUE) 
            }
        }else{
            # Heteroskedasticity
            for(regime in 1:M){
                S.iteration <- posteriors$S[,regime,iteration]
                for(j in 1:K){
                    # Prior for Si, log normal
                    log.prior.S[j,regime,iteration] <- dlnorm(S.iteration[j], meanlog=priors$S.mean, sdlog=priors$S.sd, log=TRUE) 
                }   
            }
        }
    }      
    # Sum for all iterations
    log.prior.S.joint   <- colSums(colSums(log.prior.S))


    # Correlations
    #-----------------------------------------------------------------------------
    if(debug) cat("\nPriors, correlations\n")
    for(iteration in 1:N){
        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        if(sum(restrictions$Sigma.St) == 0){
            # Homoskedasticity
            R   <-  diag(K)
            R[upper.tri(R, diag=FALSE)] <- posteriors$R[,1,iteration]
            R[lower.tri(R, diag=FALSE)] <- posteriors$R[,1,iteration]
            
            joint.density.R <- det(R)^((K*(K-1)/2)-1)
            for(i in 1:K){
                # Principal submatrices
                joint.density.R <- joint.density.R * det(matrix(R[-i,-i],nrow=K-1))^(-(K+1)/2)
            }
            log.prior.R[1,iteration]    <- log(joint.density.R)

        }else{
            # Heteroskedasticity
            for(regime in 1:M){
                R   <-  diag(K)
                R[upper.tri(R, diag=FALSE)] <- posteriors$R[,regime,iteration]
                R[lower.tri(R, diag=FALSE)] <- posteriors$R[,regime,iteration]

                joint.density.R <- det(R)^((K*(K-1)/2)-1)
                for(i in 1:K){
                    # Principal submatrices
                    joint.density.R <- joint.density.R * det(matrix(R[-i,-i],nrow=K-1))^(-(K+1)/2)
                }
                log.prior.R[regime,iteration] <- log(joint.density.R)
                }
            
        }
    }      
    # Sum for all iterations
    log.prior.R.joint   <- colSums(log.prior.R)


    # Intercepts and autoregressive coefficients
    #-----------------------------------------------------------------------------
    if(debug) cat("\nPriors, intercepts and autoregressive coefficients\n") 
        for(iteration in 1:N){
        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 

        for(regime in 1:M){
            par.x       <- posteriors$Beta[,regime,iteration]
            par.mean    <- priors$Beta.bar
            par.sigma   <- solve(priors$V.Beta.bar.m1)
            ldensity    <- dmvnorm(x=par.x, mean=par.mean, sigma=par.sigma, log=TRUE)

            log.prior.Beta[regime,iteration]   <- ldensity
        }
    }      
    # Sum for all iterations
    log.prior.Beta.joint    <-  colSums(log.prior.Beta)


    #-----------------------------------------------------------------------------------------------
    # 3/4: Truncated multivariate normal
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n3/4: Truncated multivariate normal, loop over all iterations of the Gibbs sampler")
    
    for(iteration in 1:N){
        if((debug) && ((iteration %% print.iterations)==0))  cat(" ",iteration) 
        
        # Density from truncated multivariate normal
        #-------------------------------------------------------------------------
        par.x       <- as.vector(posteriors.vector[,iteration])
        par.mean    <- as.vector(posterior.mean)
        par.sigma   <- posterior.var

        par.lower   <- truncation.interval[,1]
        par.upper   <- truncation.interval[,2]
        log.truncated.normal[iteration] <- dtmvnorm(x=par.x, mean=par.mean, sigma=par.sigma, lower=par.lower, upper=par.upper,log=TRUE)
    }




    #-----------------------------------------------------------------------------------------------
    # 4/4: Modied Harmonic Mean
    #-----------------------------------------------------------------------------------------------
    if(debug) cat("\n4/4: Computation of the MHM")
    # Equation (15), Sims, Wagoner, Zha (2008)
    iteration.product   <- array(NA, N)
    for(iteration in 1:N){
        iteration.truncated             <- log.truncated.normal[iteration]
        iteration.likelihood            <- log.likelihood[iteration]
        iteration.priors                <- log.prior.PR_TR.joint[iteration] + log.prior.S.joint[iteration] + log.prior.R.joint[iteration] + log.prior.Beta.joint[iteration]
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




MSVAR.posteriors.matrix <- function(Gibbs.output, restrictions)
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
    N           <- dim(posteriors$Beta)[3]
    Q           <- length(restrictions$dj)    
 
    # Vectorize all the posteriors for all Gibbs iterations: number of parameters * N
    # PR_TR
    number.parameters   <- sum(restrictions$dj -1)
    cat("\nPR_TR:", number.parameters, "\n")
    # S
    number.parameters   <- number.parameters + sum(diag(restrictions$Sigma.0))
    for(regime in 1:M){
        number.parameters   <- number.parameters + sum(diag(restrictions$Sigma.St[,,regime]))
    }
    cat("S:", number.parameters, "\n")
    # R
    number.parameters   <- number.parameters + sum(restrictions$Sigma.0[upper.tri(restrictions$Sigma.0, diag=FALSE)])
    for(regime in 1:M){
        R   <- restrictions$Sigma.St[,,regime]
        number.parameters   <- number.parameters + sum(R[upper.tri(R, diag=FALSE)])
    } 
    cat("R:", number.parameters, "\n")
    # Beta
    number.parameters   <- number.parameters + sum(restrictions$Beta.0) + sum(restrictions$Beta.St) 
    cat("Beta:", number.parameters, "\n")
    # Posterior vectors
    posteriors.vector   <- array(NA, c(number.parameters,N))

    for(iteration in 1:N){
        # PR_TR
        ### PR_TR.vec           <- matrix(matrix(posteriors$PR_TR[,iteration], nrow=M,byrow=TRUE),ncol=1)
		### w                   <- solve(restrictions$M) %*% PR_TR.vec
        w                   <- posteriors$w[,iteration]
        w.minus.last        <- vector()
        for(j in 1:Q){
            index.begin <- sum(restrictions$dj[1:j-1]) + 1
            dj          <- restrictions$dj[j]
            index.end   <- index.begin + dj - 1

            w.minus.last<- append(w.minus.last,w[index.begin:(index.end-1)])
        }  
        iteration.vector    <- as.vector(w.minus.last)

        # S
        s.iteration     <- vector()
        s0.restriction  <- diag(restrictions$Sigma.0)
        s0              <- posteriors$S[,1,iteration]
        s.iteration     <- append(s.iteration, s0[which(s0.restriction != 0)])
        for(regime in 1:M){
            sm.restriction  <- diag(restrictions$Sigma.St[,,regime])
            sm              <- posteriors$S[,regime,iteration]
            s.iteration     <- append(s.iteration, sm[which(sm.restriction != 0)])  
        }
        iteration.vector    <- append(iteration.vector, s.iteration)

        # R
        r.iteration     <- vector()
        r0.restriction  <- restrictions$Sigma.0[upper.tri(restrictions$Sigma.0, diag=FALSE)]
        r0              <- posteriors$R[,1,iteration]
        r.iteration 	<- append(r.iteration, r0[which(r0.restriction != 0)])
        for(regime in 1:M){
            R               <- restrictions$Sigma.St[,,regime]   
            rm.restriction  <- R[upper.tri(R, diag=FALSE)]
            rm              <- posteriors$R[,regime,iteration]
            r.iteration <-append(r.iteration, rm[which(rm.restriction != 0)]) 
        }
        iteration.vector    <- append(iteration.vector, r.iteration)

        # Beta
        b.iteration <- vector()
        beta.0 <- matrix(posteriors$Beta[,1,iteration],ncol=K)
        b.iteration <- append(b.iteration, beta.0[which(restrictions$Beta.0 != 0)])
        for(regime in 1:M){
            beta.St <- matrix(posteriors$Beta[,regime,iteration],ncol=K)
            b.iteration <- append(b.iteration, beta.St[which( restrictions$Beta.St[,,regime] != 0)]) 
        }
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
    



#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------




MSVAR.remove.posteriors<- function(Gibbs.output, keep.every)
{
    # parameters
    #-----------------------------------------------------------------------------
    aux <- Gibbs.output$last.draws
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]
    d   <- sum(restrictions$dj)   
    N   <- dim(Gibbs.output$posteriors$Beta)[3] 

    # List of vectorized posteriors for each iteration
    #----------------------------------------------------------------------------- 
    new.N   <- trunc(N/keep.every)
    posteriors <- list(
        Beta    = array(NA, c(K*(1+K*p),M,new.N)),
        Sigma   = array(NA, c(K*K,M,new.N)),
        S       = array(NA, c(K,M,new.N)),
        R       = array(NA, c(K*(K-1)/2,M,new.N)),
        PR_TR   = array(NA, c(M*M,new.N)),
        w   	= array(NA, c(d,new.N)),
        S.t     = array(NA, c(TT,new.N)),
        F.S.t   = array(NA, c(h,new.N)),
        F.Y     = array(NA, c(h,K,new.N))
        )

    # Build new posterior structure
    #-----------------------------------------------------------------------------
    for(iteration in 1:new.N){
        ### print(iteration*keep.every)
        posteriors$Beta[,,iteration]    <- Gibbs.output$posteriors$Beta[,,iteration*keep.every]
        posteriors$Sigma[,,iteration]   <- Gibbs.output$posteriors$Sigma[,,iteration*keep.every]
        posteriors$S[,,iteration]       <- Gibbs.output$posteriors$S[,,iteration*keep.every]
        posteriors$R[,,iteration]       <- Gibbs.output$posteriors$R[,,iteration*keep.every]
        posteriors$PR_TR[,iteration]    <- Gibbs.output$posteriors$PR_TR[,iteration*keep.every]
        posteriors$w[,iteration]        <- Gibbs.output$posteriors$w[,iteration*keep.every]
        posteriors$S.t[,iteration]      <- Gibbs.output$posteriors$S.t[,iteration*keep.every]
        posteriors$F.S.t[,iteration]    <- Gibbs.output$posteriors$F.S.t[,iteration*keep.every]
        posteriors$F.Y[,,iteration]     <- Gibbs.output$posteriors$F.Y[,,iteration*keep.every]
    
    }

    # Output
    #-----------------------------------------------------------------------------
    output  <- list(
                last.draws  = aux,
                posteriors  = posteriors
               )

    return(output)  
}
