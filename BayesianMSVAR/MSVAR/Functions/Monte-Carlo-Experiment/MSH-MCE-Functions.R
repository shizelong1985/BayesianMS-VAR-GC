
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#_________________________________________________________________________________________________
#
# 1: Monte Carlo: generation and estimation of N series
#_________________________________________________________________________________________________


MSH.MC.Experiment <- function(K, p, M, TT, N=1000, min.regime.occurrence, A, Sigma, PR_TR)   {


    # Store variables
    Gen.Xi	        <- array(NA, c(TT,1,N))
	Gen.SP          <- array(NA, c(TT-p,M,N))
    Gen.PR_TR       <- array(NA, c(M,M,N))
    Gen.Sigma	    <- array(NA, c(K,K,M,N))
    Gen.Y		    <- array(NA, c(TT,K,N))

	Est.Sigma	    <- array(NA, c(K,K,M,N))
    Est.A	        <- array(NA, c((1+K*p),K,N))
	Est.PR_TR       <- array(NA, c(M,M,N))
	Est.SP   	    <- array(NA, c(TT-p,M,N))
	Est.Xi		    <- array(NA, c(TT-p,1,N))
    
    mapping.mat	    <- array(NA, c(N,M))
    Est.likelihood  <- array(NA, c(1,N))

	cat("\n______________________________________________________________________________\nMonte Carlo experiment, MSH-VAR K:", K, " p:", p, "   M:", M, "   TT:", TT, "\n______________________________________________________________________________\n") 
    

    # Initialize:
    iterations      <- 1    # counter for  iterations
    iter.failures   <- 0    # counter for total failures of estimations


    while(iterations <= N)   { 
# BEGIN LOOP

        cat("Iteration",iterations,"/",N,"\t failed estimations so far",iter.failures,"\n")

        # Generation of series
        ### cat("  Generation of series\n")

        out.ser.gen     <- MSH.Series.Generation(K, p, M, TT, min.regime.occurrence, max.iter.hist, burnin=100, A, PR_TR, Sigma)
        if(out.ser.gen$return.code == 1)    {
            cat("Problem during generation of series, MC experiment aborted.\n")
            return.code <- 1
            output  <- list(return.code     = return.code)
	        return(output)
        }
        tmp.Gen.Xi      <- out.ser.gen$Gen.Xi
        tmp.Gen.Y       <- out.ser.gen$Gen.Y
        tmp.Gen.Sigma   <- out.ser.gen$Gen.Sigma
        tmp.Gen.PR_TR   <- out.ser.gen$Gen.PR_TR

        Gen.Xi[,,iterations]     <- tmp.Gen.Xi
        Gen.SP[,,iterations]     <- sp.from.Xi(tmp.Gen.Xi, M)[(p+1):TT,]
        Gen.Y[,,iterations]      <- tmp.Gen.Y

        for(regime in 1:M)  {
            Gen.Sigma[,,regime,iterations] <- tmp.Gen.Sigma[,,regime]
        }
        Gen.PR_TR[,,iterations]  <- tmp.Gen.PR_TR  

        # Estimation of series: wrapper for EM algorithm, returning 0=success, 1=failure
        ### cat("  Estimation of series\n")
        out.EM.wrapper  <- MSH.EM.wrapper(tmp.Gen.Y,M,p,tmp.Gen.Xi)
        
        # Successful estimation
        if(out.EM.wrapper$return.code == 0){
            
            Est.Sigma[,,,iterations]    <- out.EM.wrapper$Est.Sigma
            Est.A[,,iterations]	        <- out.EM.wrapper$Est.A
            Est.likelihood[,iterations] <- out.EM.wrapper$Est.likelihood
            Est.PR_TR[,,iterations]     <- out.EM.wrapper$Est.PR_TR 
            Est.SP[,,iterations]        <- out.EM.wrapper$Est.SP
            Est.Xi[,,iterations]        <- out.EM.wrapper$Est.Xi
            mapping.mat[iterations,]	<- out.EM.wrapper$mapping

        # Unsuccessful estimation 
        }else{
            # increment counter of overall unsuccessful iterations
            iter.failures   <- iter.failures + 1
        }
    
    # increment inner loop counter
    iterations      <- iterations + 1 

# END LOOP
    } 



    output <- list(	M			    = M, 
					p		        = p,
                    A		        = A,
                    Sigma		    = Sigma,
                    PR_TR	        = PR_TR,
                    Gen.Xi		    = Gen.Xi,
                    Gen.SP          = Gen.SP,
                    Gen.PR_TR       = Gen.PR_TR,
                    Gen.Sigma	    = Gen.Sigma,
                    ### Gen.Y		    = Gen.Y,
                    Est.Sigma	    = Est.Sigma,
                    Est.A		    = Est.A,
                    Est.PR_TR       = Est.PR_TR,
                    Est.SP   	    = Est.SP,
                    Est.Xi		    = Est.Xi,
                    Est.likelihood  = Est.likelihood,
					mapping		    = mapping.mat,
                    iter.failures   = iter.failures,
                    return.code     = 0
		   		)

	return(output)


}













#___________________________________________________________________________________________________
#
MSH.EM.wrapper  <- function(Y,M,p,Gen.Xi)   {

    Est.output  <- try(EM.MSH(Y, M, p, max.iter=1000, convergence=0.00001, diagonal.PR_TR=0.9, print=TRUE, plot=FALSE))
    if(inherits(Est.output, "try-error"))   {
        output  <- list(return.code = 1)
        return(output)  
    }

    ### if((Est.output$convergence.code != 0) && (Est.output$likelihood.climbing < 5))  {
    if(Est.output$convergence.code != 0)  {
        output  <- list(return.code = 1)
        return(output)
    }

    # Mapping regimes
    Gen.Xi.tmp	<- as.matrix(Gen.Xi, ncol=1)
    Gen.Xi.tmp	<- as.matrix(Gen.Xi.tmp[-(1:p),], ncol=1) # remove p rows from Gen.Sp
    Est.Xi.tmp 	<- Est.output$Xi

    out.mapping	<- mapping.Xi(Gen.Xi.tmp, Est.Xi.tmp,M) 
    if(out.mapping$return.code != 0) {
        output  <- list(return.code = 1)
        return(output)
    }

    # Good estimation, prepare variables to return
	mapping     <- out.mapping$mapping
    Est.Sigma   <- array(NA, c(K,K,M))
    Est.A	    <- array(NA, c((1+K*p),K))
	Est.PR_TR   <- array(NA, c(M,M))
    Est.SP      <- Est.output$sp


    Est.likelihood  <- Est.output$likelihood[length(Est.output$likelihood)]
    Est.A           <- Est.output$A   

    for(regime in 1:M)  {
        Est.Sigma[,,regime]    <- Est.output$Sigma[,,mapping[regime]]

        # Reorder transition probabilities 
        for(column in 1:M)  {
            row.est		<- mapping[regime]
            col.est		<- mapping[column]
            Est.PR_TR[regime,column]  <- Est.output$PR_TR[row.est,col.est] 
        }
        
        # Reorder smoothed probabilities
        Est.SP[,regime] <- Est.output$sp[,mapping[regime]]
    }
    # Calculate Xi again
    Est.Xi  <- t(t(max.col(Est.SP)))

    output  <- list(    return.code     = 0,
                        Est.Sigma       = Est.Sigma,
                        Est.A           = Est.A,
                        Est.PR_TR       = Est.PR_TR,
                        Est.SP          = Est.SP,
                        Est.Xi          = Est.Xi,
                        Est.likelihood  = Est.likelihood,
                        mapping         = mapping
                   )

    return(output)
}











#____________________________________________________________________________________________
#
# From true and estimated state vectors, match regimes by minimizing distance between the two.

mapping.Xi <- function(Gen.Xi, Est.Xi,M)
{
   # For all the observations of every state (in Gen.Xi), average estimated states (in Est.Xi)
	mapping <- matrix(0, nrow=1, ncol=M)
	
    ### print(dim(Gen.Xi))
    ### print(dim(Est.Xi))

	for(regime in 1:M){
		### cat("regime", regime, "\n")
		### print(round(mean(Est.Xi[which(Gen.Xi==regime)])))
		mapping[1,regime] <-abs(round(mean(Est.Xi[which(Gen.Xi==regime)])))
	}
	
	mapping.tmp <- as.vector(mapping)
	if(length(unique(mapping.tmp)) < M){
		cat("\t\tWarning, at least one regime is used twice for mapping!\n")
		return.code <- 1
	}else{
		return.code <- 0
	}
	
	output <- list(mapping=mapping, return.code=return.code)
	return(output)
}                         

