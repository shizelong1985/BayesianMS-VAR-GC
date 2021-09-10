
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  MSIH: M-step
#---------------------------------------------------------------------------------------------------
  
M.step.MSIH <- function(y, p, M, fp, sp, PR_TR, Sigma.old)
{
	TT  <- dim(y)[1]	# number of obs
    T   <- TT - p
	K   <- dim(y)[2]	# number of variables in VAR
    

    # Sigma from last iteration is used to calculate Beta, see Krolzig (1997), chapter 9

    # Storage 
    #_____________________________________________________________________________ 
    
    df          <- rowSums(sp)

    Sigma.new   <- array(NA, c(K, K, M))

    eta.t 	<- matrix(NA, M, T) 


    # X.bar, Xm.bar, Y 
    #_____________________________________________________________________________

    X.bar   <- matrix(data=1, nrow=T, ncol=1)	# A vector of ones for the intercept
    for(i in 1:p){
		Y.m.j   <- matrix(y[(p+1-i):(TT-i),], ncol=K)
		X.bar   <- cbind(X.bar,Y.m.j)
	}
    X.bar   <- X.bar[,-1] # Get rid of the intercepts
    
    Xm.bar  <- array(NA, c(T, M+K*p,M))
    for(regime in 1:M){
        iota.m          <- matrix(diag(M)[,regime], nrow=1)
        Xm.bar[,,regime]<- cbind(kronecker(matrix(1,nrow=T,ncol=1),iota.m), X.bar)
    }

    Y       <- matrix(y[(p+1):TT,], ncol=K)

    # y.reg
    y.reg   <- matrix(t(as.matrix(Y,ncol=K)),ncol=1,byrow=FALSE)



    # Beta
    #_____________________________________________________________________________
    LHS <- matrix(0, K*(K*p+M),  K*(K*p+M))
    MHS <- matrix(0, K*(K*p+M),  K*T)
    for(regime in 1:M){
        cross.Xm.Xi <- crossprod(Xm.bar[,,regime], diag(sp[regime,]))
        sigma.inv   <- solve(Sigma.old[,,regime])
        LHS <- LHS + kronecker(cross.Xm.Xi %*% Xm.bar[,,regime] , sigma.inv)
        MHS <- MHS + kronecker(cross.Xm.Xi, sigma.inv)
    }
    Beta        <- solve(LHS) %*% MHS %*% y.reg
    Beta.out    <- matrix(Beta, ncol=K, byrow=TRUE) 
    ### cat("Post Beta\n")


    # Residuals, Variance-covariance matrices + densities 
    #_____________________________________________________________________________
    u.tmp   <- array(NA, c(T*K,1,M))
    u       <- array(NA, c(T,K,M))

    for(regime in 1:M){

        # Residuals
        #_____________________________________________________________________________
        u.tmp[,,regime] <- y.reg - kronecker(Xm.bar[,,regime], diag(K)) %*% Beta
        u[,,regime]     <- matrix(u.tmp[,,regime], ncol=K,byrow = TRUE)

        # Variance-covariance matrices
        #_____________________________________________________________________________
        Sigma.new[,,regime] <- (crossprod(u[,,regime], diag(sp[regime,])) %*% u[,,regime]) / df[regime]  

        # Densities
        #_____________________________________________________________________________
        # Densities: matrix algebra >> slower
        ### sigma.inv		<- solve(Sigma.new[,,regime])
        ### det.sigma.m12	<- sqrt(det(as.matrix(sigma.inv))) 
        ### u.likel.tmp     <- matrix(u.tmp[,,regime], ncol=K,byrow = TRUE)
        ### eta.t[,regime]   <- apply(u.likel.tmp, 1, cond.dens, det.sigma.m12=det.sigma.m12, sigma.inv=sigma.inv, K=K)

        # Densities: Using R's normal density function >> faster
        res.tmp         <- matrix(u.tmp[,,regime], ncol=K,byrow = TRUE)
        eta.tmp         <- dmvnorm(res.tmp, sigma=matrix(Sigma.new[,,regime],ncol=K)) 
        eta.t[regime,]  <- matrix(eta.tmp, nrow=1)
    }

    # Likelihood
    #_____________________________________________________________________________    
    # xi.tp1.t and log.lik
    xi.tp1.t    <- Ergodic.PR_TR(PR_TR)
    log.lik     <- log(t(eta.t[,1]) %*% xi.tp1.t)
    for(t in 1:(T-1)){
        xi.tp1.t    <- t(PR_TR) %*% fp[,t]
        log.lik     <- log.lik +  log(t(eta.t[,t+1]) %*% xi.tp1.t)
    }   


    # Output
    #_____________________________________________________________________________

	output <- list( Beta        = Beta.out, 
                    Sigma       = Sigma.new, 
                    u           = u, 
                    likelihood  = log.lik)
    return(output)
}
