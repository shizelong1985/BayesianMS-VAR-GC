
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# MSM: EM-algorithm
#---------------------------------------------------------------------------------------------------
  
EM.MSM <- function(Y, M, p, max.iter=100, convergence=0.0001, diagonal.PR_TR=0.5, print=TRUE, plot=TRUE)
{


    if(M==1){ stop("\n\n\t -- Enter more than one regime.") }
	
	Y	<- as.matrix(Y)
	TT 	<- dim(Y)[1]
	K 	<- dim(Y)[2]
    N   <- M^(1+p)


# Initialization
#_________________________________________________________________________________________________
#
cat("\nINITIALIZATION\n______________________________________________\n")
    
    initialization  <- FALSE
	init.results 	<- try(Initialization.MSM(Y, M, p, diagonal.PR_TR))
    if(inherits(init.results, "try-error")){
	    initialization  <- TRUE
    }else{ 
        F    			<- init.results$F
        mu              <- init.results$mu
        A               <- init.results$A
        Sigma           <- init.results$Sigma
        U               <- init.results$residuals
        xi.bar          <- Ergodic.PR_TR(F)
    }
    ### cat("mu")
    ### print(mu)
    ### cat("A")
    ### print(A)
    ### cat("Sigma")
    ### print(Sigma)
    ### cat("F")
    ### print(F)



# Loop for the EM algorithm, starting with a M-step
#_________________________________________________________________________________________________
#
cat("\nEM LOOP\n____________________________________________\n")
	iter 					<- 1
    likelihood.val          <- 1 
	convergence.criterion 	<- 1
    max.change              <- 1
    expectation             <- FALSE
    maximization            <- FALSE
    parameters              <- rbind(matrix(F,ncol=1),matrix(A,ncol=1),matrix(mu,ncol=1),matrix(Sigma,ncol=1))	
	
	while( (initialization == FALSE) && ((abs(convergence.criterion) > convergence) && (max.change > convergence)) && (iter <= max.iter) && (expectation == FALSE) && (maximization == FALSE))    { 

# E-step
cat("E-step\n") 
#_________________________________________________________________________________________________
# 
        BHLK.results	<- try(BHLK.filter(u=U, PR_TR=F, Sigma=Sigma, xi.bar=xi.bar)) # Filtering
        if(inherits(BHLK.results, "try-error")){
            expectation <- TRUE
        }else{   
            fp 	    <- BHLK.results$fp
            sp	    <- BHLK.smoother(fp, F) # Smoothing


# M-step
cat("M-step\n")
#_________________________________________________________________________________________________
#           
            # Transition matrix from discrete space of smoothed probabilities
            out.PR_TR   <- sp.to.PR_TR(fp, sp, F)
            print(22)
            F           <- out.PR_TR$PR_TR
            print(F)    
            xi.bar      <- Ergodic.PR_TR(F)



            # Regression step to get the regime-specific estimates
            max.results <- try(M.step.MSM(Y=Y, p=p, M=M, sp=sp, fp=fp, A=A, mu=mu, Sigma=Sigma, F=F))
            if(inherits(max.results, "try-error")){
                    maximization <- TRUE
            }else{
                U       <- max.results$U
                Sigma   <- max.results$Sigma
                A       <- max.results$A
                mu      <- max.results$mu
            
				
# Convergence criteria
#_____________________________________________________________________________________________
#
				# Likelihood
                likelihood.val  <- cbind(likelihood.val, max.results$likelihood)
				if(likelihood.val[iter+1] <= 0){
                    convergence.criterion	<- -100 * (likelihood.val[iter+1] - likelihood.val[iter])/(likelihood.val[iter])
                }else{
                    convergence.criterion	<- 100 * (likelihood.val[iter+1] - likelihood.val[iter])/(likelihood.val[iter])
                } 

                # Change in parameters
                parametersm1<- parameters
                parameters  <- rbind(matrix(F,ncol=1),matrix(A,ncol=1),matrix(mu,ncol=1),matrix(Sigma,ncol=1))
                change      <- parameters - parametersm1
                max.change  <- abs(max(change))   
                ### cat("Iteration:", iter, "\tMax.change:", max.change, "\n")     


# Convergence Infos
#_____________________________________________________________________________________________
#
                
                ### cat("Iteration:", iter, "\tLikelihood:", likelihood.val[iter+1], "\t% change:", convergence.criterion, "\n")

                if(plot==TRUE){
                    op <- par(mfrow = c(N+1,1), mar=c(2,1,3,0.5), oma=c(1.5,2,1,1)) # N+1 pictures on one plot, and margins settings
                    # Plot of likelihood over iterations
                    plot(seq(1:(length(likelihood.val)-1)), likelihood.val[-1], type="l", main="Likelihood climbing")
                    # Plot of smoothed probabilities
                    for(regime in 1:N){
                        sp.plot <- ts(t(sp), start=c(1, 1+p), frequency=1)
                        title   <- paste("Smoothed probabilities, regime ", regime, sep="")
                        plot(sp.plot[,regime], main=title ,ylim=c(0, 1))
                    }
                    # par(ask=TRUE)
                    par(op) # At end of plotting, reset to previous settings
                }

                iter <- iter + 1

			}
		}
	}

    
# Return code for inform about the success of the EM convergence
#_________________________________________________________________________________________________
#
	if(!((abs(convergence.criterion) > convergence) && (max.change > convergence))){
		### cat("\n______________________________________\nConvergence achieved.")
		convergence.code <- 0
	}
	if(!(iter <= max.iter)){
		cat("Convergence NOT achieved within the maximum authorized number of iterations.\n")
		convergence.code <- 2 
	}
	if(!(expectation == FALSE)){
		cat("\nProblem occurred during the expectation step, EM stops.\n")
		convergence.code <- 3
	}
    if(!(maximization == FALSE)){
		cat("\nProblem occurred during the maximization step, EM stops.\n")
		convergence.code <- 4
	}   
	if(!(initialization == FALSE)){
		cat("Singularity problem during initialization, EM stops.\n")
		convergence.code <- 5
	} 	


	
# Output
#_________________________________________________________________________________________________
#
	if(convergence.code <= 1){
        # Criteria for determining whether the algorithm successfully estimated
        #_________________________________________________________________________________________________
        #
        # Idea: overall likelihood climbing result
        likelihood.climbing <- - 100 * (likelihood.val[length(likelihood.val)] - likelihood.val[2]) / likelihood.val[2]    
        ### cat("____________________________________________\nOverall % likelihood change between first and last iteration:", likelihood.climbing, "\n")
 

        # Regimes vector from smoothed probabilities 
		Xi	<- matrix(max.col(t(sp)),ncol=1)

		
		# Set up the output
		output <- list(	Sigma               = Sigma,
						A                   = A,
                        mu                  = mu,
						### PR_TR.base.m        = PR_TR,
                        PR_TR.base.n        = F,
						sp.base.n           = sp,
						sp.base.m           = Sp.base.n.to.base.m(t(sp), M, p),
						Xi.base.n           = Xi,
                        Xi.base.m           = Xi.base.n.to.base.m(Xi, M, p),
						Y                   = Y,
						U                   = U,
                        likelihood          = likelihood.val,
                        likelihood.climbing = likelihood.climbing,
						convergence.code    = convergence.code)
	}else{
		output <- list( Y                   = Y,
						convergence.code    = convergence.code)
	}
    return(output)
}
 



   











 


#___________________________________________________________________________________
#
# Counts the transitions in the discrete state space.
count.transitions.discrete <- function(sp){

	M   <- ncol(sp)
    TT  <- nrow(sp)
    Xi  <- max.col(sp)

    trans  <- matrix(0, M, M)
    for (t in 2:TT){ 
		stm1            <- Xi[t-1]
        st              <- Xi[t]
        trans[stm1,st]  <- trans[stm1, st] + 1
    }
    
	return(trans)
}   

# Counts the transitions in the continuous state space.
count.transitions.continuous <- function(sp){

	M   <- ncol(sp)
    TT  <- nrow(sp)

    trans  <- matrix(0, M, M)
    for (t in 2:TT){
        ### print(t)
        sptm1           <- matrix(sp[t-1,], nrow=1)
        spt             <- matrix(sp[t,], nrow=1)
        ### print(sptm1)
        ### print(spt)
        for(regime in 1:M){
            ### print(regime)
            trans[regime, ] <- trans[regime, ] + matrix(sptm1[,regime] * spt, nrow=1)
        }
    }
    output  <- round(trans)
	return(output)
}  


 
#___________________________________________________________________________________
# 
# Stacked version of the transition probability matrix
stacked.PR_TR   <- function(PR_TR,M,p){

    ### F   <- diag(as.vector(kronecker(vec(PR_TR), matrix(1,M^(p-1),1))))%*% kronecker(kronecker(matrix(1,M,1),diag(nrow=M^p)),matrix(1,1,M))   
    F   <- diag(as.vector(vec(kronecker(PR_TR, matrix(1,M^(p-1),1))))) %*% kronecker(kronecker(matrix(1,M,1),diag(nrow=M^p)),matrix(1,1,M))

    output  <- t(F)
    return(output)
}





#___________________________________________________________________________________
# 
# states: matrix of all possible past means 

states  <- function(M,p){
    nrows   <- M^(1+p)
    ncols   <- 1+p
    ### if(p==0) return(matrix(1:M,M,1))
    
    states  <- matrix(0, nrows, ncols)
    for(column in 1:ncols){
        divide.in   <- M^column
        block       <- nrows / divide.in
        for(i in 1:divide.in){
            row.begin   <- (i-1)*block + 1
            row.end     <- row.begin + block - 1
            states[row.begin:row.end,column]    <-i - ((i-1)%/%M) * M
        }
    }
    return(states)
}







#___________________________________________________________________________________
#
# vectorizes an (a x b) matrix y.  The resulting vector vecy has dimension (a*b x 1).

vec <- function(y){

	row <- dim(y)[1] 
	column <- dim(y)[2]
	dim(y) <- c(row*column,1)
    return(y)
}

# inverse operation...
vec.inv <- function(y,row){

    column <- dim(y)[1] / row
    dim(y)  <- c(row, column)
    return(y)

}


#___________________________________________________________________________________
# 
# From true and estimated state vectors, match regimes by minimizing distance between the two.
mapping.Xi <- function(Gen.Xi, Est.Xi,M)
{
   # For all the observations of every state (in Gen.Xi), average estimated states (in Est.Xi)
	mapping <- matrix(0, nrow=1, ncol=M)

	for(regime in 1:M){
		### cat("regime", regime, "\n")
        ### print(Est.Xi[which(Gen.Xi==regime)])
		### print(round(mean(Est.Xi[which(Gen.Xi==regime)])))
		mapping[1,regime] <-abs(round(mean(Est.Xi[which(Gen.Xi==regime)])))
	}

	mapping.tmp <- as.vector(mapping)
	if(length(unique(mapping.tmp)) < M){
		cat("Warning, at least one regime is used twice for mapping!\n")
		return.code <- 1
	}else{
		return.code <- 0
	}

	output <- list(mapping=mapping, return.code=return.code)
	return(output)
}                         

### mapping.Xi <- function(Gen.Xi, Est.Xi,M)
### {
   ### # For all the observations of every state (in Gen.Xi), average estimated states (in Est.Xi)
	### mapping <- matrix(0, nrow=1, ncol=M)
### 	
	### for(regime.est in 1:M){
		### cat("regime.est", regime.est, "\n")
		### for(regime.gen in 1:M){
            ### cat("\tregime.gen", regime.est, "\n")
            ### dist    <- Est.Xi[,regime.est]-Gen.Xi[,regime.gen]
            ### print(dist)
            ### print(sum(dist))
        ### }

	### }
### 	
    ### print("Press enter")
    ### scan()      

### 	
	### output <- list(mapping=mapping, return.code=return.code)
	### return(output)
### }        


#___________________________________________________________________________________
#
# Generate matrix of probabilities from Xi
sp.from.Xi <- function(Xi, M)
{
    TT <- dim(Xi)[1]
    
	sp <- matrix(0, nrow=TT, ncol=M)
    for (t in 1:TT){ 
		sp[t,Xi[t]] <- 1
    }
    
	return(sp)
}


#___________________________________________________________________________________
#
# Convert Xi from base N to base M
Xi.base.n.to.base.m <- function(Xi.n, M, p)
{
    TT      <- dim(Xi.n)[1]
    Xi.m    <- matrix(NA, TT, 1)

    states  <- states(M,p)

    for (t in 1:TT){
        regime.m    <- states[Xi.n[t],1]
		Xi.m[t]     <- regime.m 
    }
    
	return(Xi.m)
}


#___________________________________________________________________________________
#
# Convert Sp from base N to base M
Sp.base.n.to.base.m <- function(Sp.n, M, p)
{
    TT      <- dim(Sp.n)[1]
    Sp.m    <- matrix(NA, TT, M)

    N       <- M^(p+1)
    N.step  <- N / M  

    for(m in 1:M){
        n.begin <- (m-1) * N.step + 1
        n.end   <- n.begin + N.step - 1
        n.seq.m <- n.begin:n.end

        Sp.m[,m]    <- matrix(rowSums(Sp.n[,n.seq.m]), ncol=1)
    }
    
	return(Sp.m)
} 
