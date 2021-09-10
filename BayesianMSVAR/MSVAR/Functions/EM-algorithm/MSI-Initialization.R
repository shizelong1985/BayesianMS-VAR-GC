
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  MSI: Initial value for starting the EM algorithm  
#---------------------------------------------------------------------------------------------------
  

Initialization.MSI <- function(Y, M, p, diagonal.PR_TR)
{
    TT	<- dim(Y)[1] 
	K	<- dim(Y)[2] 
	

	OLS.AR.coeff 	<- array(data = NA, c(K*p	, K*p))
    OLS.coeff	    <- array(data = NA, c(K*p	, K*p+1 , M)) 
    OLS.coeff.out	<- array(data = NA, c(K 	, K*p+1 , M)) 
	OLS.residuals 	<- array(data = NA, c(TT-p	, K     , M))
	OLS.sigma 		<- array(data = NA, c(K		, K))  


    # OLS for AR coefficients and Sigma
    TMP             <- olsvarc.MSI(Y, p=p)
    OLS.AR.coef     <- matrix(TMP$A_orig[,-1],ncol=K*p)     # get rid of intercepts
    OLS.sigma       <- TMP$Sigma[1:K,1:K]
    ### cat("OLS.coeff\n")
    ### print(OLS.AR.coef)  
    ### cat("OLS.sigma\n")
    ### print(OLS.sigma)  




    # Intercepts: split residuals into M parts, calculate intercepts of each
    #______________________________________________________________________________________________


    ### # 1
    ### # resize sample into a multiple of M, to have equal matrices.
	### # remove some observations but better to make loops on identically sized matrices, for convenience
	### if( (TT%%M) != 0){
		### y.mat.mod <- matrix(Y[((TT%%M)+1):TT,], ncol=K)
	### } else{
		### y.mat.mod <- matrix(Y, ncol=K)
	### }     
    ### y.mat.ordered <- as.matrix(y.mat.mod)  

    ### # Sort column of data
	### column.sort <- 1
    ### y.mat.ordered <- y.mat.mod[order(y.mat.mod[,column.sort]),]

	### # Segment sample into n_states pieces
	### segment.size    <- dim(y.mat.ordered)[1] / M
	### y.mat.segmented <- array(data = NA, c(segment.size, K, M))

    ### for(regime in 1:M){
		### row.begin                   <- segment.size*(regime-1) + 1
		### row.end                     <- row.begin + segment.size - 1
		### # cat("begin", row_begin, "	end", row_end, "\n")
		### y.mat.segmented[,,regime]   <- y.mat.ordered[row.begin:row.end,]

        ### TMP                         <- olsvarc.MSI(y.mat.segmented[,,regime], p=p)
        ### OLS.inter	                <- matrix(TMP$A_orig[,1],ncol=1)
        ### cat("Regime: ",regime,"\n")
        ### print(TMP$A_orig)
        ### OLS.coeff[,,regime]         <- cbind(OLS.inter,OLS.AR.coef)
        ### OLS.coeff.out[,,regime]     <- OLS.coeff[1:K,,regime]
	### }


	# 2
    # 
    intercepts <- init.intercepts.MSI(Y, p, OLS.AR.coef, M)   

    for(regime in 1:M){
        OLS.coeff[,,regime]         <- cbind(intercepts[,,regime],OLS.AR.coef)
        OLS.coeff.out[,,regime]     <- OLS.coeff[1:K,,regime]
	}


    ### cat("OLS.coeff.out\n") 
    ### print(dim(OLS.coeff.out))
    ### print(OLS.coeff.out)



    # Residuals
    #______________________________________________________________________________________________
    # 
    for(regime in 1:M){
     	# ols_coeffs_resids() returns residuals given some OLS coefficients.
		TMP <-  ols.coeffs.resids.MSI(y=Y, p=p, A=OLS.coeff[,,regime])
		if(K==1){
			OLS.residuals[,,regime] <- t(TMP$U)
		}else{
			OLS.residuals[,,regime] <- TMP$U
		}
    }
    ### print(dim(OLS.residuals))
    ### print(OLS.residuals)  


    # Sigma : (K x K x M) format for common expectation function
    #______________________________________________________________________________________________
    #
    sigma.out   <- array(NA, c(K,K,M))
    for(regime in 1:M){
     	sigma.out[,,regime] <- OLS.sigma
    }


	# Initialization of the transition probabilities.
	#______________________________________________________________________________________________
	#
	# Symmetric automatic initialization with no limitation of number of states
    PR_TR   <- get.PR_TR(M,diagonal.PR_TR) 
      

	
	# Output    
	#______________________________________________________________________________________________
	
	output <- list( Beta    = OLS.coeff.out, 
                    Sigma   = sigma.out, 
                    u       = OLS.residuals, 
                    PR_TR   = PR_TR)
	return(output)  
}







#_________________________________________________________________________________________________
#
# This program estimates a level VAR with intercept in companion format by OLS

olsvarc.MSI <- function(y, p)
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

	output 	<- list(A=A, Sigma=SIGMA, U=U, V=V, X=X, A_orig=A_orig)
	return(output)
}      



#_________________________________________________________________________________________________
#
# This program returns the intercepts once the AR coefficients are known

init.intercepts.MSI <- function(y, p, A, M)
{

	# Phase 1
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

	X	<- Y[,1:(TT-p)]
	Y	<- Y[,2:(TT-p+1)]
	U		<- Y - A%*%X

    U.T    <- t(U) # Transpose
    ### print(U.T)


    # Phase 2
    TT  <- dim(U.T)[1] 
    K	<- dim(U.T)[2]
    
    if( (TT%%M) != 0){
		y.mat.mod <- matrix(U.T[((TT%%M)+1):TT,], ncol=K)
	} else{
		y.mat.mod <- matrix(U.T, ncol=K)
	}     
    y.mat.ordered <- matrix(y.mat.mod, ncol=K)  
	
    # Sort column of data
	column.sort <- 1
    y.mat.ordered <- matrix(y.mat.mod[order(y.mat.mod[,column.sort]),], ncol=K)
 

	# Segment sample into n_states pieces
	segment.size    <- dim(y.mat.ordered)[1] / M
	y.mat.segmented <- array(data = NA, c(segment.size, K, M))

    for(regime in 1:M){
		row.begin                   <- segment.size*(regime-1) + 1
		row.end                     <- row.begin + segment.size - 1
		### cat("begin: ", row.begin, "	end: ", row.end, "\n")
		y.mat.segmented[,,regime]   <- y.mat.ordered[row.begin:row.end,]
    }

    out.intercepts  <- array(NA, c(K, 1, M))

  
    for(regime in 1:M){
        ### out.intercepts[,1,regime]    <- array((colMeans(y.mat.segmented[,,regime])), c(K,1))
        out.intercepts[,1,regime]    <- array((apply(matrix(y.mat.segmented[,,regime],ncol=K), 2, mean)), c(K,1))
    }
    
    ### print(out.intercepts)

	output 	<- out.intercepts
	return(output)
}  


#_________________________________________________________________________________________________
#
# This program returns residuals from known OLS coefficients

ols.coeffs.resids.MSI <- function(y, p, A)
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
