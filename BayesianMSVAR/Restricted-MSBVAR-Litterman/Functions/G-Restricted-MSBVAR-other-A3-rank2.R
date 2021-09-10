
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.





G.MSBVAR.griddy.A3.Rank2 <- function(N, h, priors, restrictions, starting.values, C=1, debug=FALSE, print.iterations=50) {

    # parameters
    #-----------------------------------------------------------------------------
    aux <- starting.values
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]
    d   <- sum(restrictions$dj) 

    # List of vectorized posteriors for each iteration
    #----------------------------------------------------------------------------- 
    posteriors <- list(
        Beta    = array(NA, c(K*(1+K*p),M,N)),
        Sigma   = array(NA, c(K*K,M,N)),
        S       = array(NA, c(K,M,N)),
        R       = array(NA, c(K*(K-1)/2,M,N)),  # only upper-triangular elements 
        PR_TR   = array(NA, c(M*M,N)),
        w   	= array(NA, c(d,N)),
        S.t     = array(NA, c(TT,N)),
        F.S.t   = array(NA, c(h,N)),
        F.Y     = array(NA, c(h,K,N)),
		c		= array(NA, N)
        )
        
    # N iterations of the Gibbs
    #----------------------------------------------------------------------------- 
    for(iteration in 1:N){
        
        if((iteration %% 100)==0) cat(" ",iteration)

        # Filtering and smoothing step
        if(debug) print("Filtering and smoothing step")
        aux <- G.filtering.smoothing(aux)

        # Hidden Markov Chain step
        if(debug) print("Hidden Markov Chain step") 
        aux <- G.hidden.markov.chain.Rank2(aux, priors)

        # Standard deviations step
        if(debug) print("Standard deviation step")
        tmp <- try(G.standard.deviation(aux, priors, restrictions))
		if(!inherits(tmp, "try-error")){
	    	aux <- tmp
    	}else{
			 cat(iteration,"Standard deviation error\n")
		}

        # Correlations step
        if(debug) print("Correlation step")
        aux <- G.correlation(aux, priors, restrictions) 

        # Regression step
        if(debug) print("Regression step")
        aux <- G.regression.Beta.A3.rank2(aux, priors, restrictions, C=C) 

        # Forecasting 
        if(debug) print("Forecasting step")
        posteriors  <- G.forecast.MSVAR(aux, posteriors, h, iteration)

        # Save posteriors as vectors
        for(regime in 1:M){
            posteriors$Beta[,regime,iteration]  <- matrix(aux$Beta[,,regime], ncol=1)
            posteriors$S[,regime,iteration]     <- matrix(aux$S[,regime], ncol=1)
            R   <- aux$R[,,regime]
            posteriors$R[,regime,iteration]     <- matrix(R[upper.tri(R,diag=FALSE)],ncol=1) 
            posteriors$Sigma[,regime,iteration] <- matrix(aux$Sigma[,,regime], ncol=1)    
        }
        posteriors$PR_TR[,iteration]    <- matrix(aux$PR_TR, ncol=1)
        posteriors$c[iteration] <- aux$c
        posteriors$w[,iteration]    	<- matrix(aux$w, ncol=1)
        posteriors$S.t[,iteration]      <- matrix(max.col(t(aux$xi)), ncol=1)

        # Print iteration results
        if((iteration %% print.iterations)==0){
            cat("\n---------------------------------------------------------- \nIteration:",iteration,"/", N,"\n") 
            cat("Count transitions\n")
            print(count.regime.transitions(aux$xi))
            cat("Transition probabilities\n")
            print(aux$PR_TR)

            ### cat("S\n") 
            ### print(aux$S)
            ### cat("R\n") 

            print(aux$R)
            cat("Sigma\n")
            print(aux$Sigma)
            cat("Intercepts + AR coefficients\n")
            print(aux$Beta)
        }
    }


    # Output
    #-----------------------------------------------------------------------------
    output  <- list(
                last.draws  = aux,
                posteriors  = posteriors
               )

    return(output)
                    
}





#---------------------------------------------------------------------------------------------------
# Gibbs: hidden Markov step, drawing from the Dirichlet distribution
# Restricted matrix of transition probabilities, P, Sims et al.(2008)
# 
# Metropolis-Hastings algorithm for stationary Markov Chains (Fruehwirth-Schnatter 2006, p.341)
# For stationary Markov chains, the initial distribution is equal to the ergodic distribution. The rows of the transition matrix are not independent a posteriori. One can use a Metropolis-Hastings algorithm.
#---------------------------------------------------------------------------------------------------
  
G.hidden.markov.chain.Rank2 <- function(aux, priors){
    
    # Setup constants 
    M <- dim(aux$PR_TR)[1]
    p <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T <- dim(aux$Y)[1]
    TT <- T - p
    K <- dim(aux$Y)[2]   

    PR_TR.old <- aux$PR_TR
       
    # count transitions from xi 
    transitions <- count.regime.transitions(aux$xi)
    # incorporate priors
    alpha <- transitions + matrix(priors$w,ncol=M)


    # 0) Swap regimes in the transition probability matrix
	swapseq <- sample(M, M, replace = FALSE)
	T.alpha <- swap.PR_TR(alpha, swapseq)
	T.PR_TR.old <- swap.PR_TR(PR_TR.old, swapseq) 
    
	# Draw S0 states from previous ergodic probabilities
	ergodic.old <- G.Ergodic.PR_TR(T.PR_TR.old) 
    S0 <- sample(1:M, 1, prob=ergodic.old) 
	
	# 1) Draw from Dirichlet distribution to all rows except the last
    T.PR_TR.partial <- t(apply(T.alpha[-M,], 1, rdirichlet, n=1))  

	# 2) Last row as the linear combination of the others, maximizing the kernel f(c) (except DEoptim minimizes..)
  	outDEoptim <- DEoptim(	fn=kernel.fc, 
							lower=0, 
							upper=1, 
							control=DEoptim.control(strategy=2, trace=FALSE), 
							PR_TR.partial=T.PR_TR.partial, 
							alpha=T.alpha, 
							S0=S0)

    c <- outDEoptim$optim$bestmem
	T.PR_TR.new <- rbind(T.PR_TR.partial , c*T.PR_TR.partial[1,] + (1-c)*T.PR_TR.partial[2,]) 
    aux$c <- c

    # Ergodic probabilities
    ergodic.new <- G.Ergodic.PR_TR(T.PR_TR.new)

    # MH acceptance ratio
    AA  <- ddirichlet(T.PR_TR.new[3,], T.alpha[3,]) * ergodic.new[S0] / (ddirichlet(T.PR_TR.old[3,], T.alpha[3,]) * ergodic.old[S0] )

     # Metropolis-Hastings acceptance or rejection
    if(runif(1) <= AA ){
		# Swap regimes back to original
        aux$PR_TR   <- reverse.swap.PR_TR(T.PR_TR.new, swapseq)
    }else{
        aux$PR_TR   <- PR_TR.old
    }
# 	print(aux$PR_TR)

    # Output 
    return(aux)  

}


#---------------------------------------------------------------------------------------------------
# evaluates the kernel
#---------------------------------------------------------------------------------------------------

kernel.fc <- function (c, PR_TR.partial, alpha, S0) {
	
	# Transition prob matrix associated with c
	PR_TR <- rbind(PR_TR.partial , c*PR_TR.partial[1,] + (1-c)*PR_TR.partial[2,])
	ergodic.new <- G.Ergodic.PR_TR(PR_TR) 
	### print(PR_TR)
	# Kernel (in logs)
	powermatrix <- PR_TR ^ alpha  
	kernel <- - ( log(ergodic.new[S0]) + sum(log(powermatrix))) # DEoptim minimizes

	return(kernel)
}


#---------------------------------------------------------------------------------------------------
# swap regimes in the transition probability matrix
#---------------------------------------------------------------------------------------------------

swap.PR_TR <- function (PR_TR, swapseq) {
	M <- dim(PR_TR)[1]
	out.PR_TR <- matrix(NA, M, M)

	for (row in 1:M) {
		i <- swapseq[row] 
		for (col in 1:M) {
			j <- swapseq[col]
			out.PR_TR[i,j] <- PR_TR[row,col]
		}
	}
	return(out.PR_TR)
}

reverse.swap.PR_TR <- function (PR_TR, swapseq) {
	M <- dim(PR_TR)[1]
	out.PR_TR <- matrix(NA, M, M)

	for (row in 1:M) {
		i <- swapseq[row] 
		for (col in 1:M) {
			j <- swapseq[col]
			out.PR_TR[row,col] <- PR_TR[i,j]
		}
	}
	return(out.PR_TR)
}  



#---------------------------------------------------------------------------------------------------
# restricted regression step: nonlinear restrictions A3 rank2
#---------------------------------------------------------------------------------------------------
G.regression.Beta.A3.rank2 <- function(aux, priors, restrictions, C=1){
   
   # Setup constants 
   #-----------------------------------------------------------------------------
   M   <- dim(aux$PR_TR)[1]
   p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
   T   <- dim(aux$Y)[1]
   TT  <- T - p
   K   <- dim(aux$Y)[2]   
   
   # X, Y
   #-----------------------------------------------------------------------------
   X   <- matrix(data=1, nrow=TT, ncol=1)   # A column of ones for the intercept
   if(p > 0){
      for(lag in 1:p){
         Y.m.j   <- aux$Y[(p+1-lag):(T-lag),]
         X       <- cbind(X,Y.m.j)
      }
   }
   Y <- matrix(aux$Y[(p+1):T,],ncol=K)   
   
   X.reg   <- X
   for(regime in 1:M){
      X.reg   <- cbind(X.reg,X)
   }
   
   # Metropolis Hastings
   #-----------------------------------------------------------------------------
   ### try         <- 1
   ### accepted    <- FALSE
   ### while((try <= 1) && (accepted==FALSE)){
   
   # Candidate Draw of Beta
   #-----------------------------------------------------------------------------
   # Inverse of variance
   Sm1 <- array(NA, c(K,K,M))
   for(regime in 1:M) Sm1[,,regime]   <- solve(aux$Sigma[,,regime])
   
   
   # Crossproducts
   Zp.Sm1.Z <- matrix(0, K*(1+K*p)*(1+M), K*(1+K*p)*(1+M))
   Zp.Sm1.y <- matrix(0, K*(1+K*p)*(1+M), 1)
   
   for(t in 1:TT){
      # Restrictions at time t
      restrictions.t  <- array(NA, c((1+K*p)*(1+M),K))
      # Regime-stable
      restrictions.t[1:(1+K*p),]    <- restrictions$Beta.0
      # Regime-1, regime-2, ..., regime-M
      for(regime in 1:M){
         row.b   <- 1+ (1+K*p) + (1+K*p) * (regime - 1)
         row.e   <- row.b + (1+K*p) - 1
         restrictions.t[row.b:row.e,]  <- aux$xi[regime,t] * restrictions$Beta.St[,,regime]
      }               
      
      # Z.t matrix, including restrictions
      Z.t <- matrix(0, K, K*(1+K*p)*(1+M))
      for(k in 1:K){
         col.b   <- (1+K*p)*(1+M) * (k-1) + 1  
         col.e   <- col.b + (1+K*p)*(1+M) - 1   
         
         Z.t[k,col.b:col.e]  <- matrix(X.reg[t,],nrow=1) * restrictions.t[,k]
      }
      
      Zp.Sm1      <- crossprod(Z.t,Sm1[,,which(aux$xi[,t]==1)])
      
      Zp.Sm1.Z    <- Zp.Sm1.Z + Zp.Sm1 %*% Z.t
      Zp.Sm1.y    <- Zp.Sm1.y + Zp.Sm1 %*% matrix(Y[t,],nrow=K)
   }
   
   # Variance of Beta.0
   prior.V.Beta.bar.m1 <- kronecker(priors$V.Beta.bar.m1, diag(1+M))
   V.Beta.bar          <- solve(prior.V.Beta.bar.m1 + Zp.Sm1.Z)
   
   # Mean of Beta.0
   prior.Beta.bar  <- kronecker(priors$Beta.bar, array(1,c(1+M,1)))
   Beta.bar        <- V.Beta.bar %*% ( prior.V.Beta.bar.m1 %*% prior.Beta.bar + Zp.Sm1.y )
   
   # Draw from multivariate normal (mvtnorm package)
   draw        <- matrix(rmvnorm(n=1, mean=Beta.bar, sigma=C*V.Beta.bar, method="chol"), ncol=K)
   
   
   # Store
   #----------------------------------------------------------------------------- 
   # Beta.0 (impose restrictions again, due to priors)
   Beta.0  <- matrix(draw[1:(1+K*p),],ncol=K) * restrictions$Beta.0
   
   # Beta.St (impose restrictions again, due to priors)
   Beta.St <- array(NA,dim(restrictions$Beta.St))
   for(regime in 1:M){
      row.b   <- 1+ (1+K*p)*regime 
      row.e   <- row.b + (1+K*p) - 1
      Beta.St[,,regime] <- matrix(draw[row.b:row.e,],ncol=K) * restrictions$Beta.St[,,regime]
   }
   
   # Beta as sum of regime invariant and regime dependent draws
   Beta    <- array(NA,dim(restrictions$Beta.St)) 
   for(regime in 1:M){
      Beta[,,regime]  <- Beta.0 + Beta.St[,,regime] 
   }
   
   # Restriction A3 rank2 (iii): 
   seq.AR   = seq(3,(1+K*p),K)
   for(i in seq.AR){
      Beta[i,1,1]	<- Beta[i,1,3]*((aux$PR_TR[1,3]*aux$PR_TR[2,2]-aux$PR_TR[2,3]*aux$PR_TR[1,2])/(aux$PR_TR[1,1]*aux$PR_TR[2,2]-aux$PR_TR[2,1]*aux$PR_TR[1,2]))
   }

   # Restriction A3 rank2 (i) and (ii): 
   seq.ARR  = (1:(1+K*p))[-seq.AR]
   for(i in seq.ARR){
      Beta[i,1,1]   <- - Beta[i,1,2]*((aux$PR_TR[1,2]-aux$PR_TR[2,2])/(aux$PR_TR[1,1]-aux$PR_TR[1,2])) - Beta[i,1,3]*((aux$PR_TR[1,3]-aux$PR_TR[2,3])/(aux$PR_TR[1,1]-aux$PR_TR[2,1]))
   }
   
   
   # MH: acceptance rate
   #-----------------------------------------------------------------------------
   
   # 1/4: Likelihood old
   #----------------------------------
   # Residuals
   U   <- array(NA, dim(aux$U))
   for(regime in 1:M){
      res.tmp             <- Y - X %*% aux$Beta[,,regime]
      U[,,regime]     <- matrix(res.tmp, ncol=K)
   }   
   # Densities for each regimes
   eta.t 	<- matrix(NA, M, TT)
   for(regime in 1:M){
      eta.tmp         <- dmvnorm(matrix(U[,,regime],ncol=K), sigma=matrix(aux$Sigma[,,regime],ncol=K)) 
      eta.t[regime,]  <- matrix(eta.tmp, nrow=1)
   }

   PR.forecasted   = G.BHLK.forecasting(aux)
   # Log likelihood 
   log.likelihood.old   <- sum( log(colSums(eta.t * PR.forecasted)))

   # 2/4: Likelihood new
   #----------------------------------
   # Residuals
   U   <- array(NA, dim(aux$U))
   for(regime in 1:M){
      res.tmp             <- Y - X %*% Beta[,,regime]
      U[,,regime]     <- matrix(res.tmp, ncol=K)
   }   
   # Densities for each regimes
   eta.t 	<- matrix(NA, M, TT)
   for(regime in 1:M){
      eta.tmp         <- dmvnorm(matrix(U[,,regime],ncol=K), sigma=matrix(aux$Sigma[,,regime],ncol=K)) 
      eta.t[regime,]  <- matrix(eta.tmp, nrow=1)
   }
   # Log likelihood 
   PR.forecasted   = G.BHLK.forecasting(aux)
   # Log likelihood 
   log.likelihood.new   <- sum( log(colSums(eta.t * PR.forecasted)))
   
   # 3/4 Priors old
   #----------------------------------
   log.prior.Beta <- array(NA, M)
   for(regime in 1:M){
      par.x       <- matrix(aux$Beta[,,regime], nrow=1)
      par.mean    <- priors$Beta.bar
      par.sigma   <- solve(priors$V.Beta.bar.m1)
      ldensity    <- dmvnorm(x=par.x, mean=par.mean, sigma=par.sigma, log=TRUE)
      
      log.prior.Beta[regime]   <- ldensity
   }
   log.prior.Beta.old    <-  sum(log.prior.Beta)
   
   # 4/4 Priors old
   #----------------------------------
   log.prior.Beta <- array(NA, M)
   for(regime in 1:M){
      par.x       <- matrix(Beta[,,regime], nrow=1)
      par.mean    <- priors$Beta.bar
      par.sigma   <- solve(priors$V.Beta.bar.m1)
      ldensity    <- dmvnorm(x=par.x, mean=par.mean, sigma=par.sigma, log=TRUE)
      
      log.prior.Beta[regime]   <- ldensity
   }
   log.prior.Beta.new    <-  sum(log.prior.Beta)
   
   if (is.na(log.likelihood.new)){ log.likelihood.new = -Inf}
   if (is.na(log.prior.Beta.new)){ log.prior.Beta.new = -Inf}
   
   
   # MH: acceptance rate 
   #----------------------------------   
   AA    = exp(log.likelihood.new + log.prior.Beta.new - log.likelihood.old - log.prior.Beta.old)
#    cat("\nlog.likelihood.new: ")
#    cat(log.likelihood.new)
#    cat("\nlog.prior.Beta.new: ")
#    cat(log.prior.Beta.new)
#    cat("\nAA: ")
#    cat(AA)
   
   if(!is.numeric(AA)){AA=2e-10}
   
   # Metropolis-Hastings acceptance or rejection
   if(runif(1) <= AA ){
      aux$Beta.0  <- Beta.0
      aux$Beta.St <- Beta.St
      aux$Beta    <- Beta
      aux$U       <- U
      aux$acceptance  <- aux$acceptance + 1
   }
   
   # Output 
   #-----------------------------------------------------------------------------
   return(aux)
}                
