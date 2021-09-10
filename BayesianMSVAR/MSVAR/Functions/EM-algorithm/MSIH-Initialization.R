
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  MSIH: Initial value for starting the EM algorithm  
#---------------------------------------------------------------------------------------------------
  

Initialization.MSIH <- function(Y, M, p, diagonal.PR_TR)
{
    TT	<- dim(Y)[1] 
	K	<- dim(Y)[2] 
	
	# First, initial values of the vector of parameters are calculated. 
	#_______________________________________________________________________________________________
	#
	# The series are sorted, split in M parts on which initial conditional regressions are computed to launch the Maximum likelihood descent.  


    # resize sample into a multiple of n_states, to have equal matrices.
	# may remove some observations but better to make loops on identically sized matrices, for convenience
	if( (TT%%M) != 0){
		y_mat_mod <- Y[((TT%%M)+1):TT,] 
	} else{
		y_mat_mod <- Y
	}

    # Split the sample in M periods
	y_mat_ordered <- as.matrix(y_mat_mod)

	# Sort n-th column of data
	### n	<- 1
	### column_sort <- n
    ### y_mat_ordered <-  as.matrix(y_mat_mod[order(y_mat_mod[,column_sort]),],ncol=K)
 
	segment_size <- dim(y_mat_ordered)[1] / M
	y_mat_segmented <- array(data = NA, c(segment_size, K, M))
    for(regime in 1:M){
		row_begin 	<- segment_size*(regime-1)+1
		row_end     <- row_begin+segment_size-1
		### cat("begin", row_begin, "	end", row_end, "\n")
		y_mat_segmented[,,regime] <- y_mat_ordered[row_begin:row_end,]
	}

	# OLS on each sub-sample
	OLS_coeff 		<- array(data = NA, c(K		, K*p	, M))	 
	OLS_residuals 	<- array(data = NA, c(TT-p	, K		, M))
	OLS_sigma 		<- array(data = NA, c(K		, K		, M))
	for(regime in 1:M){
		### cat("y_mat_segmented[,,regime]", regime, "\n")
    	### print(y_mat_segmented[,,regime])
 
 		TMP <- olsvarc.MSIH(y_mat_segmented[,,regime], p=p) 
		
		A	<- as.matrix(TMP$A)
		OLS_coeff[,,regime] <-  A[1:K,]

		# ols_coeffs_resids() returns residuals given some OLS coefficients.
		TMP <-  ols_coeffs_resids.MSIH(y=Y, p=p, A=TMP$A_orig)
		if(K==1){
			OLS_residuals[,,regime] <- t(TMP$U)
		}else{
			OLS_residuals[,,regime] <- TMP$U
		}
		### cat("OLS_residuals[,,regime]", regime, "\n")
		### print(OLS_residuals[,,regime])

		OLS_sigma[,,regime] <- TMP$SIGMA
	    ### cat("OLS_sigma[,,regime]", regime, "\n")
    	### print(OLS_sigma[,,regime])
	}

	
	# Initialization of the transition probabilities.
	#_______________________________________________________________________________________________
	#
	# Symmetric automatic initialization with no limitation of number of states
    PR_TR   <- get.PR_TR(M,diagonal.PR_TR) 
      
	
	# Output    
	#_______________________________________________________________________________________________
	
	output <- list(A=OLS_coeff, Sigma=OLS_sigma, u=OLS_residuals, PR_TR=PR_TR)
	return(output)
}







#_________________________________________________________________________________________________
#
# This program estimates a level VAR with intercept in companion format by OLS

olsvarc.MSIH <- function(y, p)
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
# This program returns residuals from known OLS coefficients

ols_coeffs_resids.MSIH <- function(y, p, A)
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

	U		<- Y - A%*%X
    SIGMA	<- U%*%t(U) / (TT-p-(p*K)-1) 
	
	U		<- t(as.matrix(U[1:K,]))
    SIGMA 	<- SIGMA[1:K,1:K]

	output 	<- list(U=U, SIGMA=SIGMA)
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
