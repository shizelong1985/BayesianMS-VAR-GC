
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  MSIAH: Initial value for starting the EM algorithm  
#---------------------------------------------------------------------------------------------------
  
  
Initialization.MSIAH <- function(Y, M, p, diagonal.PR_TR)
{
    T	<- dim(Y)[1]
    TT  <- T - p
	K	<- dim(Y)[2] 
	
	# First, initial values of the vector of parameters are calculated. 
	#-----------------------------------------------------------------------------------------------
	# The series are sorted, split in M parts on which initial conditional regressions are computed to launch the Maximum likelihood descent.  


    # resize sample into a multiple of M, to have equal matrices, for convenience
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
    ### y_mat_ordered <- as.matrix(y_mat_mod[order(y_mat_mod[,column_sort]),])
 
	segment_size <- dim(y_mat_ordered)[1] / M
	y_mat_segmented <- array(data = NA, c(segment_size, K, M))
    for(regime in 1:M){
		row_begin 	<- segment_size*(regime-1)+1
		row_end     <- row_begin+segment_size-1
		y_mat_segmented[,,regime] <- y_mat_ordered[row_begin:row_end,]
	}

	# OLS on each sub-sample
	OLS_residuals 	<- array(data = NA, c(TT, K, M))
	OLS_sigma 		<- array(data = NA, c(K, K, M))
	for(regime in 1:M){
 		TMP.seg <- olsvarc.MSIAH(y=y_mat_segmented[,,regime], p=p) 
		
		# ols_coeffs_resids() returns residuals given some OLS coefficients.
		TMP <-  ols_coeffs_resids.MSIAH(y=Y, p=p, A=TMP.seg$A)
		OLS_residuals[,,regime] <- matrix(TMP$U, ncol=K) 

		OLS_sigma[,,regime] <- TMP.seg$SIGMA
	}

	
	# Initialization of the transition probabilities.
	#-----------------------------------------------------------------------------------------------
	# Symmetric automatic initialization with no limitation of number of states
    PR_TR   <- get.PR_TR(M,diagonal.PR_TR) 
      
	
	# Output    
	#-----------------------------------------------------------------------------------------------
	output <- list(Sigma=OLS_sigma, u=OLS_residuals, PR_TR=PR_TR)
	return(output)
}







#---------------------------------------------------------------------------------------------------
# Estimate a level VAR with intercept in companion format by OLS

olsvarc.MSIAH <- function(y, p)
{
	y	<- as.matrix(y)
	T   <- dim(y)[1]	# number of obs
    TT  <- T - p
	K   <- dim(y)[2]	# number of variables in VAR  
   	
    if(p == 0){
        X   <- matrix(data=1, nrow=T, ncol=1)	# A column of ones for the intercept    
        Y <- matrix(y,ncol=K)     
    }else{
        X   <- matrix(data=1, nrow=TT, ncol=1)	# A column of ones for the intercept    
        for(lag in 1:p){
            Y.m.j   <- y[(p+1-lag):(T-lag),]
            X       <- cbind(X,Y.m.j)
        }
        Y <- matrix(y[(p+1):T,],ncol=K)
    } 
	
 
	A 		<-  solve(crossprod(X)) %*% crossprod(X,Y)

	U		<- Y - X%*%A
	SIGMA	<- crossprod(U) / (TT-p-(p*K)-1) 

	output 	<- list(A       = A,
                    SIGMA   = SIGMA,
                    U       =U
                   )

	return(output)
}      

    

#---------------------------------------------------------------------------------------------------
# This program returns residuals from known OLS coefficients

ols_coeffs_resids.MSIAH <- function(y, p, A)
{
    y	<- as.matrix(y)
	T   <- dim(y)[1]	# number of obs
    TT  <- T - p
	K   <- dim(y)[2]	# number of variables in VAR  
   	
    if(p == 0){
        X   <- matrix(data=1, nrow=T, ncol=1)	# A column of ones for the intercept    
        Y <- matrix(y,ncol=K)     
    }else{
        X   <- matrix(data=1, nrow=TT, ncol=1)	# A column of ones for the intercept    
        for(lag in 1:p){
            Y.m.j   <- y[(p+1-lag):(T-lag),]
            X       <- cbind(X,Y.m.j)
        }
        Y <- matrix(y[(p+1):T,],ncol=K)
    } 

	U		<- Y - X%*%A
    SIGMA	<- crossprod(U) / (TT-p-(p*K)-1) 
	

	output 	<- list(U       = U,
                    SIGMA   = SIGMA
                   )

	return(output)
}



#---------------------------------------------------------------------------------------------------
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
