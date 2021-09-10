
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Generate MSH process
#---------------------------------------------------------------------------------------------------
# K equations
# p lags
# M states
# TT: size of series



#---------------------------------------------------------------------------------------------------
# Wrapper: Generate parameters, generate series
#

MSH.Parameters.Generation <- function(K, p, M, TT, regime.multiplier=1,A.max.iter=10000, A.low.max.eigen=0.5, A.up.max.eigen=0.99, A.ratio.zeros=1/2, A.max.values=1, PR.min.diag=0.9, PR.max.diag=0.95)  {
    

    # Generate only one set of coefficients, MSH
    out.A <- gen.A(K=K, p=p, M=1, low.max.eigen=A.low.max.eigen, up.max.eigen=A.up.max.eigen, ratio.zeros=A.ratio.zeros, max.values=A.max.values, max.iter=A.max.iter)   # Sometimes problematic function, can not generate non explosive coefficients for high K, p... Need to work on that.
    if(out.A$return.code==1){
        return.code <- 1
        A           <- NULL
        PR_TR 	    <- NULL
        Sigma       <- NULL
    }else{
        return.code <- 0
        
        A           <- array(NA, c((1+K*p),K,M))
        for(regime in 1:M){
            A[,,regime] <- out.A$A.output
        }

        PR_TR 	    <- gen.PR_TR.random(M=M, min.diag=PR.min.diag, max.diag=PR.max.diag)
        Sigma       <- gen.Sigma(K=K,M=M,regime.multiplier=regime.multiplier)
    }

    output   <- list(   return.code = return.code,
                        A           = A,
                        PR_TR 	    = PR_TR,
                        Sigma       = Sigma
                    )
    return(output)
}




MSH.Series.Generation   <- function(K, p, M, TT, min.regime.occurrence=30, max.iter.hist=200, burnin=100, A, PR_TR, Sigma) { 

    condition       <- TRUE
    counter.hist    <- 1

    while(condition == TRUE)    {
        ### cat("Gen.Xi \n")
        Gen.Xi      <- gen.Xi(PR_TR,TT+burnin)
        # Avoid generating histories without any regime changes. Each regime has to be present at least a number of time (set in the "min.regime" argument)
        regime.min	<- min(colSums(sp.from.Xi(matrix(Gen.Xi[-(1:burnin)], ncol=1), M)))
        
        if(regime.min >= min.regime.occurrence)    {
            condition   <- FALSE
            return.code <- 0
            
            ### cat("Iteration:",counter.hist,"\n")
            ### cat("Gen.U \n")
            U.tmp	   	<- MSH.gen.U(Sigma, p, M, TT+burnin)
 
            ### cat("Gen.Y \n")
            Gen.Y	    <- matrix(MSH.gen.Y(A, U.tmp, Gen.Xi)$Y, ncol=K)

            # Get rid of burnin part from series, history of regimes, residuals, and calculate Sigma
	        ### cat("Gen.Sigma \n")
            Gen.Y   <- Gen.Y[-(1:burnin),]
            Gen.Xi  <- Gen.Xi[-(1:burnin)]
            
            Gen.U       <- array(NA,c(TT,K,M))
            Gen.Sigma   <- array(NA,c(K,K,M))
            for(regime in 1:M){
                Gen.U[,,regime]     <- U.tmp$U[-(1:burnin),,regime]
		        Gen.Sigma[,,regime] <- var(Gen.U[,,regime])  
            }

            ### cat("Gen.PR_TR \n")
            Gen.PR_TR   <- count.transitions.from.Xi(matrix(Gen.Xi), M)

        }else{
            counter.hist    <- counter.hist + 1
            if(counter.hist>max.iter.hist)  {
                # Too many unsuccessful iterations
                return.code <- 1
                
                output      <- list(return.code = return.code)
                return(output)
            }
        }
    }

    output  <- list(return.code = return.code,
                    Gen.Xi      = Gen.Xi,
                    Gen.Y       = Gen.Y,
                    Gen.Sigma   = Gen.Sigma,
                    Gen.PR_TR   = Gen.PR_TR,
                    Gen.U       = Gen.U)
    return(output)
}
 






#---------------------------------------------------------------------------------------------------
MSH.gen.U <- function(Sigma, p, M, TT)
{
	K <- dim(Sigma)[1]

	U 			    <- array(NA, c(K, TT, M))
	Comp.U.return   <- array(0, c(K*p,TT,M))
 	U.return	    <- array(NA, c(TT, K, M)) 

    for(regime in 1:M){
        var.chol <- chol(Sigma[,,regime])

		for(t in 1:TT){
			res 					    <- rnorm(K, mean=0, sd=1)
			U[,t , regime] 			    <- var.chol %*% t(t(res))
			Comp.U.return[1:K,t,regime] <- t(U[,t , regime])
		}

  		U.return[,,regime] 	    <- t(U[,,regime])

	}     



	output <- list(U=U.return, Comp.U=Comp.U.return)
	return(output)
}



#---------------------------------------------------------------------------------------------------
MSH.gen.Y <- function(A, U, Xi)
{
	M 	<- dim(A)[3]
	TT 	<- length(Xi) 
    K 	<- dim(A)[2]
	p 	<- (dim(A)[1]-1) / K

	Ym1 <- matrix(runif(K*p, min=-10, max=10), ncol=1)	# Initialize with random numbers. Draws from uniform distribution on (-10, 10)

    A.tmp 		<- A.to.Comp.A(A)
    Comp.A 		<- A.tmp$Comp.A
	Comp.Int 	<- A.tmp$Comp.Int
	Comp.Y		<- matrix(0, nrow=(K*p), ncol=TT)

	for(t in 1:TT){
		regime <- Xi[t]
		Comp.Y[,t] <- as.matrix(Comp.Int[,,1]) + as.matrix(Comp.A[,,1]) %*% as.matrix(Ym1) + as.matrix(U$Comp.U[,t,regime])
		Ym1 <- Comp.Y[,t]
	}
	Y <- t(Comp.Y)[,1:K]

	output <- list(Y=matrix(Y,ncol=K), Comp.Y=Comp.Y)
    return(output)
}
