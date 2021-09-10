
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Functions common to all MSVARs
#---------------------------------------------------------------------------------------------------
      
#---------------------------------------------------------------------------------------------------
gen.A <- function(K, p, M, max.iter=10000, low.max.eigen=0.5, up.max.eigen=0.99, ratio.zeros=0.5, max.values=1)     {

	A.output				<- array(NA, c(K*p+1, K, M))
    A.lag                   <- array(NA, c(K, K, p))
    nu                      <- matrix(NA, K, 1)
	max.mod.eigen.regimes	<- matrix(0, nrow=1, ncol=M)
    min.mod.eigen.regimes	<- matrix(0, nrow=1, ncol=M)

    for(regime in 1:M)  {

		iteration   <- 1 
		condition   <- TRUE
		while(condition == TRUE)    {
            nu  <- matrix(runif(K, min=-max.values, max=max.values), ncol=1)
            
            for(lag in 1:p)     {
                A       <- matrix(0,K,K)

                if(K > 2)   {
                    diag(A) <- runif(K, min=-max.values/lag, max=max.values/lag)
                
                    A.up.tri.ind    <- which(upper.tri(A))
                    A.lo.tri.ind    <- which(lower.tri(A)) 
                    
                    sample.size     <- round((1-ratio.zeros) * length(A.up.tri.ind))
                    A.up.sample     <- sample(A.up.tri.ind, sample.size)
                    A.lo.sample     <- sample(A.lo.tri.ind, sample.size)

                    A[A.up.sample]  <- runif(sample.size, min=-max.values/lag, max=max.values/lag)
                    A[A.lo.sample]  <- runif(sample.size, min=-max.values/lag, max=max.values/lag)

                    A.lag[,,lag]    <- A 
                }else   { # not useful for 1 or 2 variables
                    A               <- matrix(runif(K*K, min=-max.values/lag, max=max.values/lag),K,K)
                    A.lag[,,lag]    <- A
                }
                
            }

            # Transform into VAR(1) process -> B, ((K*p)*(K*p)) is the outcome
            B.tmp  <- matrix(NA, K, K*p)
            for(lag in 1:p) {
                col.b               <- (lag-1) * K + 1
                col.e               <- col.b + K - 1
                B.tmp[,col.b:col.e] <- A.lag[,,lag]
            }
            # Identity matrix (K*(p-1),K*(p-1))
			tmp <- diag(K*(p-1))
			# Append to a matrix of zeros (K*(p-1),K)
			tmp <- cbind(tmp, matrix(0,nrow=K*(p-1), ncol=K))
			# Append B and tmp
			B <- rbind(B.tmp,tmp)

            ### print(B)
    

			# check for stability and reject explosive processes, for which some eigenvalues are outside the unit circle
            eigenvalues <- eigen(B,only.values=TRUE)$values
			max.mod.eigen.regimes[,regime] <- max(Mod(eigenvalues))
			min.mod.eigen.regimes[,regime] <- min(Mod(eigenvalues))
			if((max.mod.eigen.regimes[,regime] <= up.max.eigen) && (max.mod.eigen.regimes[,regime] >= low.max.eigen)) {
				condition   <- FALSE	# Get out of the while loop
                return.code <- 0
                A.output[,,regime] <- t(cbind(nu,B.tmp))
			}

 			### cat("\n-------------------------------------------\nRegime", regime, "	Iteration", iteration)
			### cat("\nA.output[,,regime]", "\n") 
			### print(A.output[,,regime])  
			### cat("regime", regime, "max eigenvalue",  max.mod.eigen.regimes[,regime], "\n")
			
			iteration <- iteration + 1
            
            if(iteration > max.iter){
                cat("Too many iterations, the algorithm not able to generate coefficients for a process with its highest eigenvalue in the required range.\n")
                return.code <- 1
                A.output    <- NULL
                output  <- list(A.output    = A.output,
                                return.code = return.code)
	            return(output)
            }  
        }    

    }
    
    if(return.code==0){
        cat("\n_______________________\nmin.mod.eigen.regimes:",min.mod.eigen.regimes, "\n")
        cat("max.mod.eigen.regimes:",max.mod.eigen.regimes, "\n_______________________\n")
    }

    output  <- list(A.output    = A.output,
                    return.code = return.code)
	return(output)
}                         


   


#---------------------------------------------------------------------------------------------------
gen.PR_TR.random <- function(M, min.diag=0.5, max.diag=1)  {

	PR_TR <- matrix(0, nrow=M, ncol=M)

    # Diagonal elements first
    diag(PR_TR)  <- runif(M, min=min.diag, max=max.diag)   
	
    for(row in 1:M){

        for(column in sample(1:M))  {
			if(row != column)   {
                rowsum              <- sum(PR_TR[row,])
                PR_TR[row,column]   <- runif(1, min=0, max=1-rowsum)
            }
		}
        
        if(row != M){
            PR_TR[row,M]    <- PR_TR[row,M] + 1 - sum(PR_TR[row,])
        }else{
            PR_TR[row,M-1]    <- PR_TR[row,M-1] + 1 - sum(PR_TR[row,])
        }
	}

	return(PR_TR)
}    


gen.PR_TR.diagonals <- function(diagonal.values)  {
    # argument: vector of diagonals
    
    M       <- length(diagonal.values)
	PR_TR   <- matrix(NA, nrow=M, ncol=M)

    # Diagonal elements first
    diag(PR_TR)  <- diagonal.values

    # Off diagonals get equal weight
    
    off.diagonal.values <- (seq(1,1,M) - diagonal.values)/(M-1)
	
    for(row in 1:M){

        for(column in 1:M)  {
            if(row != column)PR_TR[row,column]   <- off.diagonal.values[row]
		}
	}

	return(PR_TR)
}
 
#---------------------------------------------------------------------------------------------------
# regime_multiplier: the variance-covariance of each regime is multiplied by a random number, sampled from 0 to regime_multiplier 

gen.Sigma <- function(K, M, regime.multiplier)
{
	Sigma <- array(NA, c(K, K, M))

	for(regime in 1:M){
		random.int <- runif(1, min=0, max=regime.multiplier)
		Sigma[,,regime] <- genPositiveDefMat("eigen",dim=K)$Sigma * random.int
	}

	return(Sigma)
}



#---------------------------------------------------------------------------------------------------
gen.Xi <- function(PR_TR, TT)
{   
	M 	<- dim(PR_TR)[1] 
	Xi 	<- matrix(NA, nrow=TT, ncol=1)
     
	Xi[1] <- sample(1:M, 1)
	
	for(t in 2:TT){
		state.tm1 	<- Xi[t-1]
		Xi[t] 		<- sample(1:M, 1, prob = PR_TR[state.tm1,])
	}

	return(Xi)
}


#---------------------------------------------------------------------------------------------------
# Transform coefficients into companion form

A.to.Comp.A <- function(A)
{
	K	<- dim(A)[2]
	p	<- (dim(A)[1] - 1) / K
    M	<- dim(A)[3]

    Comp.A		<- array(0, c(K*p, K*p, M))
	Comp.Int    <- array(0, c(K*p, 1, M))

	for(regime in 1:M)
	{
		Comp.Int[1:K,,regime]	<- t(A[1,,regime])
		Comp.A[1:K,,regime]  	<- t(A[2:(K*p+1),,regime])
		if(p > 1) Comp.A[(K+1):(K*p),1:(K*(p-1)),regime] <- diag(K*(p-1))
	}

	output	<- list(Comp.A=Comp.A, Comp.Int=Comp.Int)
	return(output)
}





#---------------------------------------------------------------------------------------------------
# Generate matrix of probabilities from Xi

sp.from.Xi <- function(Xi, M)
{
    TT <- dim(matrix(Xi,ncol=1))[1]
    
	sp <- matrix(0, nrow=TT, ncol=M)
    for (t in 1:TT){ 
		sp[t,Xi[t]] <- 1
    }
    
	return(sp)
}


# Counts the transitions in the discrete state space.
count.transitions.from.Xi <- function(Xi, M)
{
    TT <- length(Xi)
	PR_TR <- matrix(0, M, M)
    for (t in 2:TT){ 
		Xi_tm1 <- Xi[t-1]
        Xi_t <- Xi[t]
        PR_TR[Xi_tm1,Xi_t] <- PR_TR[Xi_tm1,Xi_t] + 1
    }
    
	PR_TR <- PR_TR/rowSums(PR_TR)
	return(PR_TR)
}                   
