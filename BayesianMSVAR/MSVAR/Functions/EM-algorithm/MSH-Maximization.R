
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  MSH: M-step
#---------------------------------------------------------------------------------------------------

M.step.MSH <- function(y, p, M, Sigma.old, fp, sp, PR_TR)
{
	TT  <- dim(y)[1]	# number of obs
    T   <- TT - p
	K   <- dim(y)[2]	# number of variables in VAR
    

    # Sigma from last iteration is used to calculate Beta, see Krolzig (1997), chapter 9

    # Storage  
    #_____________________________________________________________________________ 
    
    df          <- rowSums(sp)

    Beta        <- matrix(NA, K*p+1, K)
    Sigma.new   <- array(0, c(K, K, M))   

    eta.t 	    <- matrix(NA, M, T) 

    # X.bar, Y 
    #_____________________________________________________________________________

    X.bar   <- matrix(data=1, nrow=T, ncol=1)	# A row of ones for the intercept

	for(i in 1:p){
		Y.m.j   <- matrix(y[(p+1-i):(TT-i),], ncol=K)
		X.bar   <- cbind(X.bar,Y.m.j)
	}
    Y <- matrix(y[(p+1):TT,], ncol=K)

    # y.reg
    y.reg   <- matrix(t(as.matrix(Y,ncol=K)),ncol=1)

    # Beta 
    #_____________________________________________________________________________ 
    LHS <- matrix(0, K*(K*p+1),  K*(K*p+1))
    MHS <- matrix(0, K*(K*p+1),  K*T)
    for(regime in 1:M){

        cross.X.Xi  <- crossprod(X.bar, diag(sp[regime,]))
        sigma.inv   <- solve(Sigma.old[,,regime])
        
        LHS <- LHS + kronecker(cross.X.Xi %*% X.bar , sigma.inv)
        MHS <- MHS + kronecker(cross.X.Xi, sigma.inv)
    }
    Beta    <- solve(LHS) %*% MHS %*% y.reg
    Beta.out    <- matrix(Beta, ncol=K, byrow=TRUE)  
     
    # Residuals 
    #_____________________________________________________________________________ 
    u   <- y.reg - kronecker(X.bar, diag(K)) %*% Beta
    u   <- matrix(u, ncol=K,byrow = TRUE)

    # Variance-covariance matrices + densities 
    #_____________________________________________________________________________ 
    for(regime in 1:M){
        Sigma.new[,,regime] <- (crossprod(u, diag(sp[regime,])) %*% u) / df[regime]

        # Densities: matrix algebra >> slower 
        ### sigma.inv		<- solve(Sigma.new[,,regime])
        ### det.sigma.m12	<- sqrt(det(as.matrix(sigma.inv)))   
        ### pYSt[,regime]   <- apply(u, 1, cond.dens, det.sigma.m12=det.sigma.m12, sigma.inv=sigma.inv, K=K)   

        # Densities: Using R's normal density function >> faster
        eta.tmp         <- dmvnorm(u, sigma=matrix(Sigma.new[,,regime],ncol=K)) 
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
    # Because of common functions for E step, return residuals in (T x K x M) format
    u.out   <- array(NA, c(T,K,M))
    for(regime in 1:M){
        u.out[,,regime] <- u
    }

	output <- list( Beta        = Beta.out, 
                    Sigma       = Sigma.new, 
                    u           = u.out,
                    likelihood  = log.lik)
    return(output)
}
