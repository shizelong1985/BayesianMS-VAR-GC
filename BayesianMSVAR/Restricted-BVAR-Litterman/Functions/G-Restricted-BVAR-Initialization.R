
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  VAR: Initial value for starting the Gibbs sampler
#---------------------------------------------------------------------------------------------------
  
G.Initialization.BVAR <- function(Y, p, G)
{   
    Y   <- as.matrix(Y)
	K	<- dim(Y)[2] 
	
    
	# OLS
    #----------------------------------------------------------------------------- 
    TMP <- G.olsvarc.VAR(Y, p) 

    A   <- matrix(TMP$A, ncol=(1+K*p))
    
    OLS_coeff <-  matrix(A[1:K,], ncol=K, byrow=TRUE) 

    OLS_residuals <- matrix(TMP$U[1:K,], ncol=K, byrow=TRUE)

    OLS_sigma <- TMP$SIGMA[1:K,1:K]


	# Decomposition of covariance into standard deviation and correlation
    #-----------------------------------------------------------------------------
    S   <- apply(OLS_residuals, 2, sd)
    R   <- cor(OLS_residuals)
	
	# Output 
    #-----------------------------------------------------------------------------     
    Gibbs.input <- list(
        Y       = Y,
        U       = OLS_residuals,
        Beta    = (OLS_coeff),
        Sigma   = OLS_sigma,
        S       = S,
        R       = R,
        G       = G
        ) 

    return(Gibbs.input)
}







#---------------------------------------------------------------------------------------------------
# Estimate a level VAR with intercept in companion format by OLS

G.olsvarc.VAR <- function(y, p)
{
	y	<- as.matrix(y)
	
	TT	<- dim(y)[1] 
    K	<- dim(y)[2]
	y	<- t(y)

	
	if(p<=1){
    	Y	<- as.matrix(y)
	}else{
		Y	<- y[,p:TT]
		for(i in 1:(p-1)){
	    	Y	 <- rbind(Y, y[,(p-i):(TT-i)])
		}
	}

	if(p>=1){
        X	<- rbind(matrix(1,nrow=1,ncol=TT-p), Y[,1:(TT-p)])
        Y	<- Y[,2:(TT-p+1)] 
    }else{
        X	<- matrix(1,nrow=1,ncol=TT)      
    }
	
 
	A 		<- (Y%*%t(X)) %*% solve(X%*%t(X))
    A_orig	<- A

	U		<- Y - A%*%X
	SIGMA	<- U%*%t(U) / (TT-p-(p*K)-1) 
	V		<- A[,1]
	A		<- A[,1:(K*p+1)]

	output 	<- list(A=A, SIGMA=SIGMA, U=U, V=V, X=X, A_orig=A_orig)
	return(output)
}      


