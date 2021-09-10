
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  MSIA: M-step
#---------------------------------------------------------------------------------------------------
 
M.step.MSIA <- function(y, p, M, fp, sp, PR_TR) 
{
	TT  <- dim(y)[1]	# number of obs
    T   <- TT - p
	K   <- dim(y)[2]	# number of variables in VAR


    # Storage  
    #_____________________________________________________________________________ 
    
    Beta        <- array(NA, c(K*(1+K*p),1, M))
    Beta.out    <- array(NA, c(1+K*p, K, M))

    eta.t 	<- matrix(NA, M, T)


    # X.bar, Y
    #_____________________________________________________________________________

    X.bar   <- matrix(data=1, nrow=T, ncol=1)	# A row of ones for the intercept
    for(i in 1:p){
		Y.m.j   <- matrix(y[(p+1-i):(TT-i),], ncol=K)
		X.bar   <- cbind(X.bar,Y.m.j)
	}
    Y       <- matrix(y[(p+1):TT,], ncol=K)
    # y.reg
    y.reg   <- matrix(t(as.matrix(Y,ncol=K)),ncol=1)


    # Beta, regime-dependent  
    #_____________________________________________________________________________
    #

    for(regime in 1:M){
        Xi.m    <- diag(sp[regime,])
        crossp  <- crossprod(X.bar,Xi.m)
        Beta[,,regime]      <- kronecker(solve(crossp %*% X.bar) %*% crossp, diag(K)) %*% y.reg
        Beta.out[,,regime]  <- matrix(Beta[,,regime],ncol=K,byrow=TRUE)
    }   


    # Residuals 
    #_____________________________________________________________________________
    #
    
    beta.cols   <- matrix(Beta[,,1], nrow=K)
    for(regime in 2:M){
        beta.cols   <-  cbind(beta.cols, matrix(Beta[,,regime],nrow=K))
    }
    LHS <- kronecker(array(1,c(M,1)),Y)
    RHS <- kronecker(diag(M), X.bar) %*% t(beta.cols)
    U   <- LHS - RHS
    # format residuals from (MT x K) into (T x K x M)
    u   <- array(NA, c(T, K, M))
    for(regime in 1:M){
        row.begin   <- 1 + (regime-1) * T
        row.end     <- row.begin + T - 1
        u[,,regime] <- U[row.begin:row.end,]
    }


    # Sigma  
    #_____________________________________________________________________________
    #

    # Xi.s
    xi      <- matrix(t(sp),ncol=1)
    df      <- as.vector(matrix(1,1,M*T) %*% xi)
    Xi      <- diag(as.vector(xi))

    Sigma   <- ( t(U) %*% Xi %*% U ) / df
    
    Sigma.out   <- array(NA,c(K,K,M))
    for(regime in 1:M){
        Sigma.out[,,regime] <- Sigma
    } 


    # Densities  
    #_____________________________________________________________________________
    #

    # Densities: matrix algebra >> slower  
    ### sigma.inv		<- solve(Sigma)
    ### det.sigma.m12	<- sqrt(det(as.matrix(sigma.inv)))  
    ### for(regime in 1:M){
        ### res.tmp         <- matrix(u[,,regime], T, K)  
        ### pYSt[,regime]   <- apply(res.tmp, 1, cond.dens, det.sigma.m12=det.sigma.m12, sigma.inv=sigma.inv, K=K)    
    ### } 
 
    # Densities: Using R's normal density function >> faster
    for(regime in 1:M){
        res.tmp         <- matrix(u[,,regime], T, K)
        eta.tmp         <- dmvnorm(res.tmp, sigma=matrix(Sigma,ncol=K)) 
        eta.t[regime,]  <- matrix(eta.tmp, nrow=1)   
    } 


    # Likelihood  
    #_____________________________________________________________________________
    #  
    # xi.tp1.t and log.lik
    xi.tp1.t    <- Ergodic.PR_TR(PR_TR)
    log.lik     <- log(t(eta.t[,1]) %*% xi.tp1.t)
    for(t in 1:(T-1)){
        xi.tp1.t    <- t(PR_TR) %*% fp[,t]
        log.lik     <- log.lik +  log(t(eta.t[,t+1]) %*% xi.tp1.t)
    }  

    #_____________________________________________________________________________

	output <- list( Beta        = Beta.out, 
                    Sigma       = Sigma.out, 
                    u           = u,
                    likelihood  = log.lik)
    return(output)
}
