
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  MSI: M-step
#---------------------------------------------------------------------------------------------------
 
M.step.MSI <- function(y, p, M, fp, sp, PR_TR)
{
	TT  <- dim(y)[1]	# number of obs
    T   <- TT - p
	K   <- dim(y)[2]	# number of variables in VAR


    # Storage 
    #_____________________________________________________________________________ 
    
    Beta    <- matrix(NA, K*p+1, K)
    Sigma   <- array(0, c(K, K))

    eta.t 	<- matrix(NA, M, T)    

    
    # X.bar, Y: cast everything as matrices for univariate series
    #_____________________________________________________________________________

    X.bar   <- matrix(data=1, nrow=T, ncol=1)	# A column of ones for the intercept
	for(i in 1:p){
		Y.m.j   <- matrix(y[(p+1-i):(TT-i),], ncol=K)
		X.bar   <- cbind(X.bar,Y.m.j)
	}
    X.bar   <- matrix(X.bar[,-1], ncol=K*p)   # Get rid of the first column of 1..
   
    Y       <- matrix(y[(p+1):TT,], ncol=K) 

    # Xi.big, Xi
    Xi.big  <- t(sp)
    tmp     <- as.vector(matrix(1,1,T) %*% Xi.big)
    Xi      <- diag(tmp)

    # Beta 
    #_____________________________________________________________________________
    #

    # LHS
    tmp.up  <- cbind( Xi, crossprod(Xi.big, X.bar) )
    tmp.down<- cbind( crossprod(X.bar, Xi.big), crossprod(X.bar, X.bar) )
    LHS     <- rbind(tmp.up, tmp.down)

    # RHS
    tmp.up  <- crossprod(Xi.big, Y)
    tmp.down<- crossprod(X.bar, Y)
    RHS     <- rbind(tmp.up, tmp.down)      
   
    # Beta
    Beta    <- solve(LHS) %*% RHS



    
    # Residuals  
    #_____________________________________________________________________________
    #

    # Z.bar
    Z.bar   <- cbind( kronecker( diag(M), matrix(1,T,1) ) , kronecker( matrix(1,M,1), X.bar ) )
  
    U       <- kronecker( matrix(1,M,1) , Y ) - Z.bar %*% Beta

    # u: (T x K x M) format
    u       <- array(NA, c(T,K,M))
    for(regime in 1:M){
        row.begin   <- (regime - 1)* T + 1
        row.end     <- row.begin + T - 1
        u[,,regime] <- U[row.begin:row.end,]
    }

    
    # Sigma  
    #_____________________________________________________________________________
    #
     
    # Xi.s
    xi.s    <- matrix(Xi.big,ncol=1)
    Xi.s    <- diag(as.vector(xi.s))
    # DF
    df      <- T

    Sigma   <- ( t(U) %*% Xi.s %*% U ) / df
    
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

    # likelihood  
    #_____________________________________________________________________________
    #  
    # xi.tp1.t and log.lik
    xi.tp1.t    <- Ergodic.PR_TR(PR_TR)
    log.lik     <- log(t(eta.t[,1]) %*% xi.tp1.t)
    for(t in 1:(T-1)){
        xi.tp1.t    <- t(PR_TR) %*% fp[,t]
        log.lik     <- log.lik +  log(t(eta.t[,t+1]) %*% xi.tp1.t)
    }
 
    
    # Output
    #_____________________________________________________________________________
    #
	output <- list( Beta        = Beta, 
                    Sigma       = Sigma.out, 
                    u           = u, 
                    likelihood  = log.lik)
    return(output)
}
