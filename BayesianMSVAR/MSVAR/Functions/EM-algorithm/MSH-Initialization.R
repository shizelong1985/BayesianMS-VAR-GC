
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#
#  MSH: Initial value for starting the EM algorithm  
#---------------------------------------------------------------------------------------------------
  

Initialization.MSH <- function(Y, M, p, diagonal.PR_TR)
{
    TT	<- dim(Y)[1] 
	K	<- dim(Y)[2] 
	
	OLS.coeff 		<- array(data = NA, c(K		, K*(p+1)))	 
	OLS.residuals 	<- array(data = NA, c(TT-p	, K))
	OLS.sigma 		<- array(data = NA, c(K		, K		, M))  



    # Intercepts and AR coefficients
    TMP             <- olsvarc.MSH(Y, p=p)
    ### print(TMP)
    OLS.coeff	    <- matrix(TMP$A_orig[1:K,], nrow=K)
    OLS.residuals   <- matrix(TMP$U[1:K,], ncol=K)
    ### print(dim(OLS.residuals))
    ### print(OLS.residuals)
    ### print(OLS.coeff)


    # Variance / Covariances matrices: split residuals into M parts, calculate variance of each
    #______________________________________________________________________________________________
	#
    # resize sample into a multiple of M, to have equally sized matrices.
    T.res   <- dim(OLS.residuals)[1]
	if( (T.res%%M) != 0){
	    res.mod <- matrix(OLS.residuals[((T.res%%M)+1):T.res,] , ncol=K)
	} else{
		res.mod <- matrix(OLS.residuals, ncol=K)
	}
    ### print(dim(res.mod))
    ### print(res.mod)      

    # Sort column of data
	### column.sort <- 1
    ### res.mod <- matrix(res.mod[order(res.mod[,column.sort]),], ncol=K)      
    
    # Segment sample into M pieces
	segment.size <- dim(res.mod)[1] / M
    for(regime in 1:M){
		row.begin 	<- segment.size*(regime-1) + 1
		row.end     <- row.begin + segment.size - 1
		### cat("begin", row_begin, "	end", row_end, "\n")
		residuals   <- res.mod[row.begin:row.end,]
        sigma	    <- t(residuals)%*%residuals / (segment.size-p-(p*K)-1)

        OLS.sigma[,,regime]    <- sigma
	}
    ### print(dim(OLS.sigma))
    ### print(OLS.sigma) 


	
	# Initialization of the transition probabilities.
	#_______________________________________________________________________________________________
	#
	# Symmetric automatic initialization with no limitation of number of states
    PR_TR   <- get.PR_TR(M,diagonal.PR_TR) 
      

	
	# Output    
	#_______________________________________________________________________________________________

    u.out   <- array(NA, c(TT-p,K,M))
    for(regime in 1:M){
        u.out[,,regime] <- OLS.residuals
    }  

	output <- list(Beta=OLS.coeff, Sigma=OLS.sigma, u=u.out, PR_TR=PR_TR)
	return(output)  
}







#_________________________________________________________________________________________________
#
# This program estimates a level VAR with intercept in companion format by OLS

olsvarc.MSH <- function(y, p)
{
	y	<- as.matrix(y)
	
	TT	<- dim(y)[1] 
    K	<- dim(y)[2]
	y	<- t(y)

	
	if(p==1){
    	Y	<- as.matrix(y)
	}else{
		Y	<- y[,p:TT]
		for(i in 1:(p-1)){
	    	Y	 <- rbind(Y, y[,(p-i):(TT-i)])
		}
	}

	X	<- rbind(matrix(1,nrow=1,ncol=TT-p), Y[,1:(TT-p)])
	Y	<- Y[,2:(TT-p+1)]
 
	A 		<- (Y%*%t(X)) %*% solve(X%*%t(X))
    A_orig	<- A

	U		<- Y - A%*%X
	SIGMA	<- U%*%t(U) / (TT-p-(p*K)-1) 
	V		<- A[,1]
	A		<- A[,2:(K*p+1)]

	output 	<- list(A=A, SIGMA=SIGMA, U=U, V=V, X=X, A_orig=A_orig)
	return(output)
}      





#_________________________________________________________________________________________________
#
# Initialization of the transition probabilities.
# Symmetric automatic initialization with no limitation of number of states
#                                    
get.PR_TR   <- function(M, diagonal.PR_TR){
    
    ksi <- diagonal.PR_TR

	if(M==2){
		tmp		<- log(ksi/(1-ksi))
        tmp2    <- log((1-ksi)/ksi)
 		PR_TR	<- cbind(tmp, tmp2)  
	} else{
		delta	<- (1-ksi)/(M-1)
 		chi		<- 1+(M-2)*(delta/ksi)
  
		tmp		<- log(ksi/(1-ksi*chi))
					# gives pii=f(thetaii), invariant probabilities : exp(p)/(1+exp(p)+(M-2)*exp(q))=ksi
		q 		<- tmp+log(delta/ksi)
					# pij=f(thetaij), for j=1..M-1 : exp(q)/(1+exp(p)+(M-2)*exp(q))=delta
		z 		<- log(delta/(1-(M-1)*delta))
					# z gives piMf(thetaiM), for i=1..M-1 : exp(z)/(1+(M-1)*exp(z))=delta
		A 		<- (tmp-q)*diag(M-1)+q*matrix(1, nrow=M-1, ncol=M-1)
					# First block of thetaij, for j=1..M-1
		B 		<- z*matrix(1, nrow=M-1, ncol=1)
					# Second block of thetaij, for j=M
		PR_TR 	<- cbind(A, B)
	}
   
	# Constraining n_states-1 transition  probabilities, the n_states-th to sum to one
	PR_TR 	<- exp(PR_TR)
	Sum 	<- matrix(1, nrow=1, ncol=M) + rowSums(t(PR_TR))
	for(i in 1:(M-1)){
		PR_TR[i,] <- PR_TR[i,] / Sum
	}
  
	# Add final row
	last_row	<- matrix(1,nrow=1,ncol=M) - colSums(PR_TR)
    PR_TR		<- rbind(PR_TR,last_row) 

    return(PR_TR)
}
