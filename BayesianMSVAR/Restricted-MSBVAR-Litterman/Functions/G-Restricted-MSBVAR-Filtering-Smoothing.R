
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Gibbs: filtering and smoothing step
#---------------------------------------------------------------------------------------------------
 
G.filtering.smoothing   <- function(aux, iterations=501){

    # Setup constants 
    #----------------------------------------------------------------------------- 
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]  

    # BHLK filter 
    #----------------------------------------------------------------------------- 
    fp <- G.BHLK.filtering(aux)
    sp <- G.BHLK.smoothing(aux, fp)

    xi_T_T  <- sp[,TT]
    xi.tmp  = aux$xi
    # draw from xi(T|T)
    #----------------------------------------------------------------------------- 
    draw        <- sample(1:M, 1, prob=xi_T_T) 
    xi.tmp[,TT] <- as.matrix(diag(M)[,draw], ncol=1)

    # backwards loop 
    #----------------------------------------------------------------------------- 

    # Introduced Min. No. of appearances of staying in same state occurrences
	mnas = 3
	nas = 2
	iteration   <- 1
	while (nas<mnas & iteration<iterations) {
	   if((iteration %% 50)==0) cat(iteration, " ")
	    for(j in 1:(TT-1)){
	        xi_Tmj <- (aux$PR_TR %*% (xi.tmp[,TT-j+1]/(t(aux$PR_TR)%*%fp[,TT-j]))) * fp[,TT-j]
	        # draw from xi(T-j|T)
	        draw            <- sample(1:M, 1, prob=xi_Tmj) 
	        xi.tmp[,TT-j]   <- as.matrix(diag(M)[,draw], ncol=1)
     	}
    transitions.matrix  <- count.regime.transitions(xi.tmp) 
	nas <- min(diag(transitions.matrix))
	iteration <- iteration + 1
	}      
	if (iteration < iterations){
	   aux$xi   = xi.tmp
	} else {
	   cat("\n")
	}
	# Output 
    #----------------------------------------------------------------------------- 
    return(aux)
}




#---------------------------------------------------------------------------------------------------
# BHLK filter function
#---------------------------------------------------------------------------------------------------

G.BHLK.filtering <- function(aux){

    # Table 9.3 of Krolzig (1997)
    
    # Setup constants 
    #----------------------------------------------------------------------------- 
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   


	# Data recipients 
    #-----------------------------------------------------------------------------
    eta.t 	<- matrix(NA, M, TT)
    xi.t.t 	<- matrix(NA, M, TT)

 
    # Conditional densities for each regime 
    #-----------------------------------------------------------------------------
    for(regime in 1:M){
        # Sigma:
        sigma           <- matrix(aux$Sigma[,,regime], ncol=K)
        res.tmp         <- matrix(aux$U[,,regime], ncol=K)

        eta.t[regime,]  <- matrix(dmvnorm(res.tmp, sigma=sigma), nrow=1)
        # Copied from MSBVAR package: replace zeros by small numbers in order to avoid 
        # divisions by 0 in the forward iteration of the filtering, see below
        eta.t[regime,]  <- ifelse(eta.t[regime,]==0, 1e-300, eta.t[regime,])
    }

	# Smoothed prob or ergodic prob (after initialization only) for starting the algorithm 
    #-----------------------------------------------------------------------------
    xi.tm1.tm1   <- G.Ergodic.PR_TR(aux$PR_TR)   
    
    # Forward iteration 
    #-----------------------------------------------------------------------------
    for(t in 1:TT){
		num 		<- eta.t[,t] * (t(aux$PR_TR) %*% xi.tm1.tm1)
        den 		<- as.numeric(matrix(1,1,M) %*% num)
        xi.t.t[,t] 	<- num / den
        xi.tm1.tm1  <- xi.t.t[,t]
    }

    # Output 
    #-----------------------------------------------------------------------------
	return(xi.t.t)
}    


#---------------------------------------------------------------------------------------------------
# BHLK forecatsing function
#---------------------------------------------------------------------------------------------------
G.BHLK.forecasting <- function(aux){
   
   # equation 22.4.6 of Hamilton's book, page 692
   
   # Setup constants 
   #----------------------------------------------------------------------------- 
   M   <- dim(aux$PR_TR)[1]
   p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
   T   <- dim(aux$Y)[1]
   TT  <- T - p
   K   <- dim(aux$Y)[2]   
   
   
   # Data recipients 
   #-----------------------------------------------------------------------------
   eta.t    <- matrix(NA, M, TT)
   xi.t.t 	<- matrix(NA, M, TT)
   xi.t.tm1 <- matrix(NA, M, TT)
   
   
   # Conditional densities for each regime 
   #-----------------------------------------------------------------------------
   for(regime in 1:M){
      # Sigma:
      sigma           <- matrix(aux$Sigma[,,regime], ncol=K)
      res.tmp         <- matrix(aux$U[,,regime], ncol=K)
      
      eta.t[regime,]  <- matrix(dmvnorm(res.tmp, sigma=sigma), nrow=1)
      # Copied from MSBVAR package: replace zeros by small numbers in order to avoid 
      # divisions by 0 in the forward iteration of the filtering, see below
      eta.t[regime,]  <- ifelse(eta.t[regime,]==0, 1e-300, eta.t[regime,])
   }
   
   # Smoothed prob or ergodic prob (after initialization only) for starting the algorithm 
   #-----------------------------------------------------------------------------
   xi.tm1.tm1   <- G.Ergodic.PR_TR(aux$PR_TR)   
   
   # Forward iteration 
   #-----------------------------------------------------------------------------
   for(t in 1:TT){
      xi.t.tm1[,t] = t(aux$PR_TR) %*% xi.tm1.tm1
      num 		<- eta.t[,t] * (t(aux$PR_TR) %*% xi.tm1.tm1)
      den 		<- as.numeric(matrix(1,1,M) %*% num)
      xi.t.t[,t] 	<- num / den
      xi.tm1.tm1  <- xi.t.t[,t]
   }
   
   # Output 
   #-----------------------------------------------------------------------------
   return(xi.t.tm1)
}    




#---------------------------------------------------------------------------------------------------
# BHLK smoother function
#---------------------------------------------------------------------------------------------------
G.BHLK.smoothing   <- function(aux, fp){

    # Table 9.3 of Krolzig (1997) 
    
    # Setup constants 
    #----------------------------------------------------------------------------- 
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   

    # Data recipients 
    #-----------------------------------------------------------------------------
    sp  <- matrix(0, M, TT)

    # Last: smoothed = filtered 
    #-----------------------------------------------------------------------------
    sp[,TT]    <- fp[, TT]

    # Backwards iteration 
    #-----------------------------------------------------------------------------
    for(t in (TT-1):1){
        sp[,t]	<- (aux$PR_TR %*% (sp[,t+1]/(t(aux$PR_TR)%*%fp[,t]))) * fp[,t]
    }

    # Output 
    #-----------------------------------------------------------------------------
	return(sp) 
}




#---------------------------------------------------------------------------------------------------
# Ergodic probabilities
#---------------------------------------------------------------------------------------------------
G.Ergodic.PR_TR <- function(PR_TR){

	# computes equilibrium stationary ergodic probabilities P * (M,1) 
	# See Hamilton(1994) , section 22.2 
	M 			<- dim(PR_TR)[1]
	A 			<- rbind((diag(M)-t(PR_TR)),matrix(1, nrow=1, ncol=M))
	EN 			<- rbind(matrix(0, nrow=M, ncol=1), 1) 
	PR_ergodic 	<- solve(crossprod(A)) %*% t(A) %*% EN
    
	return(PR_ergodic)
}
