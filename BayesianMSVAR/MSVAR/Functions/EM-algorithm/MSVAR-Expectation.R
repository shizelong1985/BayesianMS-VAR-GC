
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# BHLK filter function
#---------------------------------------------------------------------------------------------------
#
# ARGUMENTS:    u       = residuals for the regression model in each regime (T x K x M)
#                         (models with regime invariant variances must call function with M #                          identical variances-covariances matrices) 
#               PR_TR   = transition probability matrix (M x M)
#               Sigma   = covariance of the residuals (K x K x M)
#                         (models with regime invariant variances must call function with M #                          identical variances-covariances matrices)
#               Xi_0_0  = initial state

# OUTPUT:       filtered probabilities

BHLK.filter <- function(u, PR_TR, Sigma, xi.bar){ 
    
    # Setup constants
    TT 	<- dim(u)[1]
    M 	<- dim(PR_TR)[1]
    K 	<- dim(u)[2]

	# Data recipients
    eta.t 	<- matrix(NA, M, TT)
    xi.t.t 	<- matrix(NA, M, TT)

 
    # Conditional densities for each regime
    for(regime in 1:M){
        res.tmp         <- matrix(u[,,regime], TT, K)     
        eta.tmp         <- dmvnorm(res.tmp, sigma=matrix(Sigma[,,regime],ncol=K))
        eta.t[regime,]  <- matrix(eta.tmp, nrow=1)        
    }

	# Smoothed prob or ergodic prob (after initialization only) for starting the algorithm
    xi.tm1.tm1   <- xi.bar    
    
    for(t in 1:TT){
		num 		<- eta.t[,t] * (t(PR_TR) %*% xi.tm1.tm1)
        den 		<- as.numeric(matrix(1,1,M) %*% num)
        xi.t.t[,t] 	<- num / den
        xi.tm1.tm1  <- xi.t.t[,t]
    }

	output <- list(fp = xi.t.t)
	return(output)
}    



#---------------------------------------------------------------------------------------------------
# BHLK smoother function
#---------------------------------------------------------------------------------------------------
BHLK.smoother   <- function(fp, PR_TR){

    TT 	<- ncol(fp)
    M 	<- nrow(fp)
    sp  <- matrix(0, M, TT)

    sp[,TT]	<- fp[,TT]

    for(t in (TT-1):1){
        sp[,t]	<- (PR_TR %*% (sp[,t+1]/(t(PR_TR)%*%fp[,t])))  * fp[,t]
    }

    return(sp)
}                             





  
#---------------------------------------------------------------------------------------------------
# Transition probabilities
#---------------------------------------------------------------------------------------------------
sp.to.PR_TR <- function(fp, sp, PR_TR)  {
   
    M   <- dim(PR_TR)[1]
    TT  <- dim(sp)[2]

    vec.P   <- matrix(PR_TR, ncol=1)
    
    xi.2.t.T    <- matrix(NA, M*M, TT-1)    
    xi.2        <- matrix(0, M*M, 1)


    # Xi.2 collection of smoothed regime probabilities
    for(t in 1:(TT-1)){
        xi.t.t      <- fp[,t]
        xi.tp1.t    <- t(PR_TR) %*% xi.t.t
        xi.tp1.T    <- sp[,t+1]

        xi.2.t.T[,t]<- vec.P * ( kronecker( (xi.tp1.T / xi.tp1.t) , xi.t.t ) )

        xi.2        <- xi.2 + xi.2.t.T[,t]
    }
    
    xi.1    <- kronecker(matrix(1,1,M), diag(M)) %*% xi.2
    
    # Formula 6.14 of Krolzig (1997), p.100
    rho     <- xi.2 / kronecker(matrix(1,M,1) , xi.1)

    PR_TR.return   <- matrix(rho, ncol=M, byrow=FALSE)
    
    output  <- list(    PR_TR       = PR_TR.return)
    return(output)
}    





#---------------------------------------------------------------------------------------------------
# Conditional density: (8) in Bellone
#---------------------------------------------------------------------------------------------------
### cond.dens   <-function(residuals, det.sigma.m12, sigma.inv, K) {

	### f <- (2*pi)^(-K/2) * det.sigma.m12 * exp(-1/2*(t(residuals)%*%sigma.inv%*%residuals))

	### return(f)
###}





#---------------------------------------------------------------------------------------------------
# Ergodic probabilities
#---------------------------------------------------------------------------------------------------
Ergodic.PR_TR <- function(PR_TR){

	# computes equilibrium stationary ergodic probabilities P * (M,1) 
	# See Hamilton(1994) , chap 22 
	M 			<- dim(PR_TR)[1]
	A 			<- rbind((diag(M)-PR_TR),matrix(1, nrow=1, ncol=M))
	EN 			<- rbind(matrix(0,nrow=M,ncol=1),1) 
	PR_ergodic 	<- solve(crossprod(A)) %*% t(A) %*% EN
    
	return(PR_ergodic)
}
