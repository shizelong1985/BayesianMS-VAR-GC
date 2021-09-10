
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



# Initialization


Initialization.MSM <- function(Y, M, p, diagonal.PR_TR){

    TT	<- dim(Y)[1] 
	K	<- dim(Y)[2]

    N   <- M^(p+1)
	
	# First, initial values of the vector of parameters are calculated. 
	#_________________________________________________________________________________________________
	#
	# The series are sorted, split in N parts on which initial conditional regressions are computed to launch the Maximum likelihood descent.  


    # resize sample into a multiple of n_states, to have equal matrices.
	# may remove some observations but better to make loops on identically sized matrices, for convenience
	if( (TT%%M) != 0){
		y_mat_mod <- matrix(Y[((TT%%M)+1):TT,],ncol=K)
	} else{
		y_mat_mod <- matrix(Y,ncol=K)
	}
    # Divide the sample in M periods by chronology
	y_mat_ordered <- y_mat_mod

	# Sort n-th column of data
	n	            <- 1
	column_sort     <- n
    y_mat_ordered   <- matrix(y_mat_mod[order(matrix(y_mat_mod[,column_sort],ncol=K),decreasing=FALSE),],ncol=K)

	# Segment sample into M pieces
	segment_size <- dim(y_mat_ordered)[1] / M
	y_mat_segmented <- array(data = NA, c(segment_size, K, M))
    for(regime in 1:M){
		row_begin 	<- segment_size*(regime-1)+1
		row_end     <- row_begin+segment_size-1
		### cat("begin", row_begin, "	end", row_end, "\n")
		y_mat_segmented[,,regime] <- y_mat_ordered[row_begin:row_end,]
	}

	# OLS on each piece
    OLS_sigma 		<- array(NA, c(K   	, K     ))
	OLS_AR   		<- array(NA, c(K   	, K     , p ))	 
 	OLS_AR.comp   	<- array(NA, c(K   	, K*p   ))	 
  	OLS_nu   		<- array(NA, c(K	, 1     , M ))	 
	OLS_mu   		<- array(NA, c(K	, 1     , M ))	 
    OLS_residuals 	<- array(NA, c(TT-p , K		, N ))
 


    # AR and Sigma are regime invariant
    #_________________________________________________________________________________________________
	out.ols     <- olsvarc(Y,p)
    OLS_AR      <- out.ols$A
    OLS_AR.comp <- out.ols$A.comp
    OLS_sigma   <- out.ols$Sigma



    # Mean is regime variant
    #_________________________________________________________________________________________________
	for(regime in 1:M){
		### cat("y_mat_segmented[,,regime]", regime, "\n")
    	### print(y_mat_segmented[,,regime])
 		out.ols             <- olsvarc(y_mat_segmented[,,regime], p) 
        ### print(out.ols)
		OLS_nu[,,regime]    <- out.ols$nu
        OLS_mu[,,regime]    <- out.ols$mu
	}
    ### cat("OLS_mu\n")
    ### print(OLS_mu)
    
    # Generate M^(1+p) series of residuals, for all possible histories of regimes
    #_________________________________________________________________________________________________
	states  <- states(M,p)
    ### print(states)
    for(n in 1:N){
        # ols_coeffs_resids() returns residuals given some OLS coefficients and past means
        out.ols             <- ols_coeffs_resids(y=Y, p=p, A=OLS_AR.comp, mu=OLS_mu, states=states[n,])
        OLS_residuals[,,n]  <- out.ols$U
        ### cat("OLS_residuals[,,regime]", regime, "\n")
        ### print(OLS_residuals[,,regime])
    }
    ### print(OLS_residuals)
 

    ### cat("OLS mu \n")
    ### print(OLS_mu) 
    ### cat("OLS A \n")
    ### print(OLS_AR)
    ### cat("OLS sigma \n")
    ### print(OLS_sigma)  


	# Initialization of the transition probabilities.
	F   <- get.F(M, p, diagonal.PR_TR)$F
                                
    
    # Sigma : (K x K x M) format for common expectation function
    #______________________________________________________________________________________________
    #
    sigma.out   <- array(NA, c(K,K,N))
    for(regime in 1:N){
     	sigma.out[,,regime] <- OLS_sigma
    }         
	
	# Output    
	#_________________________________________________________________________________________________
	#
	output <- list( A           = OLS_AR, 
                    mu          = OLS_mu, 
                    Sigma       = sigma.out, 
                    residuals   = OLS_residuals, 
                    F           = F
                  )

	return(output)
}













#____________________________________________________________________________________________________________
#
# Estimates a level VAR with intercept in companion format by OLS

olsvarc <- function(y, p){

	y	<- as.matrix(y)
	
	TT	<- dim(y)[1] 
    K	<- dim(y)[2]
	y	<- t(y)

    # A as a 3D array, or A.comp as the companion form matrix
    A   <- array(NA, c(K, K, p))
	
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
 
	A.comp 	<- (Y%*%t(X)) %*% solve(X%*%t(X))
    A_orig	<- A.comp

	U		<- Y - A.comp%*%X
	sig.tmp	<- U%*%t(U) / (TT-p-(p*K)-1)
    Sigma   <- sig.tmp[1:K,1:K]
	nu		<- A.comp[1:K,1]
	A.comp  <- matrix(A.comp[1:K,2:(K*p+1)],nrow=K)

    
    A       <- Comp.A.to.A(A.comp)
    # means, p.11 on Krolzig 1997 
    LHS <- diag(K)
    for(j in 1:p){
        LHS         <- LHS - A[,,j]
    }
    mu  <- solve(LHS) %*% nu

    ### cat("OLS nu \n")
    ### print(nu)
    ### cat("OLS mu \n")
    ### print(mu) 
    ### cat("OLS A \n")
    ### print(A)
    ### cat("OLS sigma \n")
    ### print(Sigma)   
	
    output 	<- list(    A.comp  = A.comp, 
                        A       = A,
                        Sigma   = Sigma,
                        U       = U,
                        X       = X,
                        A_orig  = A_orig,
                        nu      = nu,
                        mu      = mu    )
	return(output)
}      



#____________________________________________________________________________________________________________
# 
# returns residuals from known OLS coefficients

ols_coeffs_resids <- function(y, p, A, mu, states){

    # companion form? NO, loop
    y	<- as.matrix(t(y))
	TT	<- dim(y)[2]
    K	<- dim(y)[1]
    
    U   <- matrix(NA, K, TT-p)

	for(t in (p+1):TT){
          tmp<- y[,t]- mu[, , states[1]]
          
          for(lag in 1:p){
              col.begin <- (lag-1)*K + 1
              col.end   <- col.begin + K - 1
              tmp <- tmp - A[, col.begin:col.end] %*% (y[,t-p]- mu[, , states[lag+1]])
          }
          ### cat("t",t,"tmp",tmp,"\n")
          U[,t-p] <- tmp 
    }
	
    ### cat("states",states,"\n")
    ### print(U)

    output 	<- list(U=t(U))
	return(output)
}



# Initialization of the transition probabilities.
#_________________________________________________________________________________________________
#
# Symmetric automatic initialization with no limitation of number of states
get.F   <- function(M, p, diagonal.PR_TR){

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

    # Constraining M-1 transition probabilities, the M-th to sum to one
    PR_TR 	<- exp(PR_TR)
    Sum 	<- matrix(1, nrow=1, ncol=M) + rowSums(t(PR_TR))
    for(i in 1:(M-1)){
        PR_TR[i,] <- PR_TR[i,] / Sum
    }

    # Add final row
    last_row	<- matrix(1,nrow=1,ncol=M) - colSums(PR_TR)
    PR_TR		<- rbind(PR_TR,last_row)

    #________________________________
    #
    F   <- stacked.PR_TR(PR_TR,M,p)

    output  <- list(    PR_TR = PR_TR,
                        F     = F)
    return(output)
}

