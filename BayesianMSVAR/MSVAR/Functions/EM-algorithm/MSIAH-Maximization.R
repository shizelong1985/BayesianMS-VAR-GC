
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  MSIAH: M-step
#---------------------------------------------------------------------------------------------------
  
M.step.MSIAH <- function(y, p, M, fp, sp, PR_TR)
{
	T   <- dim(y)[1]	# number of obs
    TT  <- T - p
	K   <- dim(y)[2]	# number of variables in VAR
    
    # Storage  	
    #------------------------------------------------------------------------------   
    df      <- rowSums(sp)   
    
    Beta    <- array(NA, c(1+K*p, K, M))
    Sigma   <- array(NA, c(K,K,M))

    u       <- array(NA, c(TT, K, M))

    eta.t 	<- matrix(NA, M, TT)


    # X.bar, Y, X 
    #------------------------------------------------------------------------------  
    X   <- matrix(data=1, nrow=TT, ncol=1)	# A column of ones for the intercept
    if(p > 0){
        for(lag in 1:p){
            Y.m.j   <- y[(p+1-lag):(T-lag),]
            X       <- cbind(X,Y.m.j)
        }
    }
    Y <- matrix(y[(p+1):T,],ncol=K)   
     
    
    # Beta, residuals, Sigma, densities 
    #------------------------------------------------------------------------------ 
    for(regime in 1:M){
        cross.LHS     <- crossprod(X, diag(sp[regime,]))%*%X
        cross.RHS     <- crossprod(X, diag(sp[regime,]))%*%Y

        # Beta
        Beta[,,regime]  <- solve(cross.LHS, cross.RHS)

        # Residuals and Sigma (based on Krolzig)
        u[,,regime]     <- matrix(Y - X%*%Beta[,,regime], ncol=K)
        res.tmp         <- matrix(u[,,regime], ncol=K)
        Sigma[,,regime] <- (crossprod(res.tmp,diag(sp[regime,]))%*%res.tmp)/df[regime]


        # Densities:
        eta.tmp         <- dmvnorm(res.tmp, sigma=matrix(Sigma[,,regime],ncol=K)) 
        eta.t[regime,]  <- matrix(eta.tmp, nrow=1)   
    }
    

    # Likelihood 
    #------------------------------------------------------------------------------
    # xi.tp1.t and log.lik
    xi.tp1.t    <- Ergodic.PR_TR(PR_TR)
    log.lik     <- log(t(eta.t[,1]) %*% xi.tp1.t)
    for(t in 1:(TT-1)){
        xi.tp1.t    <- t(PR_TR) %*% fp[,t]
        log.lik     <- log.lik +  log(t(eta.t[,t+1]) %*% xi.tp1.t)
    }  

    
    # Output
    #------------------------------------------------------------------------------
	output <- list( Beta        = Beta, 
                    Sigma       = Sigma, 
                    u           = u, 
                    likelihood  = log.lik)
    return(output)
}
