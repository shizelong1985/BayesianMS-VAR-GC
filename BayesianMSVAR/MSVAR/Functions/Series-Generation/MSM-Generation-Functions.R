
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Generate MSM process
#---------------------------------------------------------------------------------------------------

MSM.Series.Generation <- function(mu, A, U, Xi, M, burn)   {

	K 	<- dim(A)[1]
	p 	<- dim(A)[3]
	TT 	<- dim(Xi)[1]

    states      <- states(M,p)
    residuals   <- t(matrix(U$U,ncol=K))

    Y   <- matrix(0, K, TT)
    for(t in 1:p){
        regime.t.base.n <- Xi[t] 
        regime.t.base.m <- states[regime.t.base.n,1]  
        
        # set first p values to the mean of their regime, plus add some residuals
        # for initializing the process
        # values will be dropped later
        Y[,t]   <- mu[,,regime.t.base.m] + residuals[,t] 
    }

    for(t in (p+1):TT){
        regime.t.base.n <- Xi[t]
        regime.t.base.m <- states[regime.t.base.n,1]
        Y.t             <- mu[,,regime.t.base.m] + residuals[,t] # mean plus residuals
        
        for(lag in 1:p){
            regime.lag.base.n   <- states[regime.t.base.n,lag+1]
            regime.lag.base.m   <- states[regime.lag.base.n,1]

            ### cat("\nt",t,"regime.t.base.n",regime.t.base.n,"\n")
            ### cat("states",matrix(states[regime.t.base.n,],nrow=1),"\n")
            ### cat("   lag",lag,"regime.lag.base.n",regime.lag.base.n,"A lag",A[,,lag],"mu lag",mu[,,regime.lag.base.n],"\n")

            A.lag       <- A[,,lag]
            Y.lag       <- Y[,t-lag]
            mu.lag      <- mu[,,regime.lag.base.m]
            Y.t         <- Y.t  + A.lag %*% (Y.lag - mu.lag)
        }
        Y[,t]   <- Y.t
    }

    # remove first burn values
    Y   <- Y[,-(1:burn)]

    output  <- t(matrix(Y,nrow=K))
    return(output)
}        




#_________________________________________________________________________________________________
#
MSM.gen.A   <- function(K, p, max.iter=1000)    {

	A   <- matrix(NA, K, K*p)

    iteration <- 1 
    condition <- TRUE
    while(condition == TRUE){
        
        # Each coefficient will be decreasing in absolute value compared to the one from the previous lag... why not! may decrease explosiveness of processes with more lags
        A[,1:K] <- matrix(runif(K*K, min=-2, max=2), ncol=K)
        ### A[,1:(K+1),regime] <- matrix(rnorm(K*(K+1), mean = 0, sd = 0.5), ncol=K+1)

        if(p>1){
            for(lag in 1:(p-1)){
                ### cat("lag", lag, "\n")
                column.begin    <- (K*lag)+1
                column.end      <- column.begin + K - 1
                for(column in column.begin:column.end){
                    ### cat("column", column, "\n")
                    for(row in 1:K){
                        ### cat("row", row, "\n")
                        # 3 ideas

                        # 1: normal distr
                        ### A[row,column,regime] <- rnorm(1, mean = 0, sd = 0.5) 
                        
                        # 2: Force the coefficients to change signs and decrease at each lag
                        ### pm1.value <-  A[row,(column-K),regime]
                        ### if(pm1.value <= 0){
                          ### min.val = pm1.value
                          ### max.val = 0
                        ### }else{
                          ### min.val = 0
                          ### max.val = pm1.value
                        ### }
                        ### A[row,column,regime] <- runif(1, min=min.val, max=max.val)
                        
                        
                        # 3: Force the coefficients to be decreasing in the lags
                        pm1.value <-  abs(A[row,(column-K)])
                        ### cat("pm1.value", pm1.value, "\n")
                        A[row,column] <- runif(1, min=-pm1.value, max=pm1.value)
                    }
                }
            }
        }  

        # Transform into VAR(1) process -> B, ((K*p)*(K*p)) is the outcome
        # Remove constant
        B <- A[,1:(K*p)]
        # Identity matrix (K*(p-1),K*(p-1))
        tmp <- diag(K*(p-1))
        # Append to a matrix of zeros (K*(p-1),K)
        tmp <- cbind(tmp, matrix(0,nrow=K*(p-1), ncol=K))
        # Append B and tmp
        B <- rbind(B,tmp)           

        ### print(B)

        # check for stability and reject explosive processes, for which some eigenvalues are outside the unit circle
        eigenvalues <- eigen(B)$values
        max.mod.eigen.regimes <- max(Mod(eigenvalues))
        if(max.mod.eigen.regimes < 0.95) {
            condition <- FALSE	# Get out of the while loop
            ### print(A)
            A.output    <- Comp.A.to.A(A) 
        }

        ### cat("\n-------------------------------------------\nRegime", regime, "	Iteration", iteration)
        ### cat("\nA[,,regime]", "\n") 
        ### print(A[,,regime])  
        ### cat("regime", regime, "max eigenvalue",  max.mod.eigen.regimes[,regime], "\n")
        
        iteration <- iteration + 1
        
        if(iteration > max.iter){
            stop("Too many iterations, algorithm not able to generate non explosive coefficients.")
        }  
    }    

    ### cat("max.mod.eigen.regimes",max.mod.eigen.regimes, "\n")        

    output  <- A.output
	return(output)
}                        


#_________________________________________________________________________________________________
# 
Comp.A.to.A <- function(Comp.A) {

	K	<- dim(Comp.A)[1]
	p	<- dim(Comp.A)[2] / K

    A   <- array(0, c(K, K, p))

	for(lag in 1:p) {
		col.begin   <- (lag-1)*K + 1
        col.end     <- col.begin + K - 1
        A[,,lag]    <- Comp.A[,col.begin:col.end]
	}

	output	<- A
	return(output)
}



   
#_________________________________________________________________________________________________
#
gen.mu <- function(K, M, mean.range)     {

    mu  <- array(NA,c(K,1,M))

    for(k in 1:K){
        mean.min    <- -mean.range
        mean.max    <- mean.range
        
        for(regime in 1:M){ 
            mu[k,,regime]  <- matrix(runif(1, min=mean.min, max=mean.max), ncol=1)

            if(mu[k,,regime] <= 0){
                mean.max <- mean.range
                mean.min <- 0
            }else{
                mean.max <- 0
                mean.min <- -mean.range
            }
        }
    }

	return(mu)
} 






#_________________________________________________________________________________________________
#
gen.PR_TR <- function(M)    {
	PR_TR <- matrix(NA, nrow=M, ncol=M)

	for(regime in 1:M){
		iteration <- 1 
		condition <- TRUE
		while(condition == TRUE){
			
			M_m_1 <- matrix(runif(M-1, min=0, max=1), ncol=M-1)
			
			# Random weight on diagonal term
			weight <- runif(1, min=0.5, max=1)
			if(rowSums(M_m_1) < (1-weight)){
				condition <- FALSE	# Get out of the while loop
				PR_TR[regime,1:(M-1)] <- M_m_1
                PR_TR[regime,M] <- (1-weight) - rowSums(M_m_1)
				Diag_weight <- diag(M) * weight
				PR_TR[regime,] <- PR_TR[regime,] + Diag_weight[regime,]
			}

			### cat("--------------------\nRegime", regime, "	Iteration", iteration, "\n")
			iteration <- iteration + 1
		}
	}

	return(PR_TR)
}    





#_________________________________________________________________________________________________
#
Xi.base.n.to.Xi.base.m  <- function(Xi.base.n,M,p)  {
    
    ### TT      <- dim(Xi.base.n)[1]
    ### N       <- M^(1+p)
    states  <- states(M,p)
    
    Xi.base.m   <- apply(Xi.base.n,2,which.Xi.base.n,states=states)
    return(Xi.base.m)
}


which.Xi.base.n <- function(Xi.base.m, states)  {
    
    Xi.base.n   <- states[Xi.base.m,1]
    return(Xi.base.n)
}




#_________________________________________________________________________________________________
#
gen.U <- function(Sigma, p, TT) {

	K <- dim(Sigma)[1]

	U 			<- array(NA, c(K, TT))
	U.return	<- array(NA, c(TT, K)) 
	Comp.U 		<- array(0, c(K*p,TT))
	Sigma.new 	<- array(NA, c(K, K))

    var.chol <- chol(Sigma)

    for(t in 1:TT){
        res 			<- matrix(rnorm(K, mean=0, sd=1), ncol=1)
        U[,t ] 			<- var.chol %*% res
        Comp.U[1:K,t] 	<- t(U[,t])
    }

    U.return[,] 	<- t(U[,])
    Sigma.new[,] 	<- var(U.return[,])

	output <- list(U=U.return, Comp.U=Comp.U, Sigma=Sigma.new)
	return(output)
}
     






#_________________________________________________________________________________________________
#
# returns residuals

residuals.max.step <- function(U, M, p, xi, alpha, mu)  {

	TT	<- dim(Y)[1]
    T   <- TT - p
    K	<- dim(Y)[2]
    N   <- M^(1+p)

    U   <- array(NA,c(T,K,N))
    
    states  <- states(M,p) 

    for(n in 1:N){
        regime.Y    <- states[n,1]
        mu.begin    <- (regime.Y-1)*K + 1
        mu.end      <- mu.begin -1 + K
        mu.Y        <- matrix(mu[mu.begin:mu.end,],ncol=1)

        ### print(states[n,-1])
        X.bar.star  <- X.bar.star(Y,p,mu,states[n,-1])

        ### cat("dimensions\n")
        ### print(dim(X.bar.star))
        ### print(dim(alpha))
        ### print(dim(t(alpha)))
        ### print(dim(kronecker(matrix(1,T,1), t(mu.Y))))

        U[,,n]  <- Y[-(1:p),] - X.bar.star %*% alpha - kronecker(matrix(1,T,1), t(mu.Y))
    }
          
    ### cat("U\n")
    ### print(U)

    output 	<- U
	return(output)
}




#_________________________________________________________________________________________________
#
# Estimate means of generated series

estimate.means <- function(mu, A, U, M, p, burn)  {

    
    TT	    <- dim(U$U)[1]
    T       <- TT - burn
    N       <- M^(1+p)
    seq.N   <- 1:N 
    states  <- states(M,p)
    
    OLS_mu  <- array(NA, c(K	, 1     , M )) 


    for(regime in 1:M){

        M.base.n        <- drop(return.M.base.n(states, regime))
        ### cat("\tregime", regime, "state", M.base.n, "\n")  
        Xi.base.n       <- matrix(M.base.n,TT,1)
        Y	            <- gen.MSM(mu, A, U, Xi.base.n, M, burn)
        
        plot(ts(Y))     
        out.ols         <- olsvarc(Y,p)
        OLS_mu[,,regime]<- out.ols$mu
        ### print(OLS_mu[,,regime])
        ### print("Press enter")
        ### scan()   
    }

    output 	<- OLS_mu
	return(output)
}      



# recursive function
return.M.base.n <- function(states.mat,regime)    {

    col.states  <- dim(states.mat)[2]
    ### cat("states.mat\n")
    ### print(states.mat)

    if(col.states > 1)  {
        candidates  <- which(states.mat[,1]==regime)
        states.mat[-candidates,]    <- 0
        ### cat("states.mat\n")
        ### print(states.mat)  
        states.rec  <- matrix(states.mat[,-1],ncol=col.states-1)
        output      <- return.M.base.n(states.rec, regime)
    }else   {
        ### cat("last observation\n")
        output      <-which(states.mat[,1]==regime)
        return(output)
    }
    return(output)
}

