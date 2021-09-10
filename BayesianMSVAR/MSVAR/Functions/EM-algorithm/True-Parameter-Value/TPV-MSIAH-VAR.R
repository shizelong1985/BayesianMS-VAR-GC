
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  EM with true parameter values instead of initialization
#---------------------------------------------------------------------------------------------------
  
TPV.EM.MSIAH <- function(Y, M, p, TPV.PR_TR, TPV.Sigma, TPV.U, TPV.Xi, max.iter=1000, convergence=0.000001) 
{   
    
    if(M==1){ stop("\n\n\t -- Enter more than one regime.") }
	
	Y	<- as.matrix(Y)
	TT 	<- dim(Y)[1] 
	K 	<- dim(Y)[2]

# Initialization
#_________________________________________________________________________________________________
#
### cat("\nINITIALIZATION\n______________________________________________\n")
	

	    initialization  <- FALSE
        PR_TR 			<- TPV.PR_TR
        Sigma           <- TPV.Sigma
        u               <- TPV.U
        xi.bar          <- Ergodic.PR_TR(PR_TR) 

# Loop for the EM algorithm, starting with a M-step
#_________________________________________________________________________________________________
#
### cat("\nEM LOOP\n____________________________________________\n")
	iter 					<- 1
    likelihood.val          <- 1
    qps.val                 <- 1
	convergence.criterion 	<- 1
    max.change              <- 1
    expectation             <- FALSE
    maximization            <- FALSE 
    number.parameters       <- M^2 + M*K*(1+K*p) + M*K^2
	parameters              <- matrix(0, number.parameters, 1) 


	while( (initialization == FALSE) && ((abs(convergence.criterion) > convergence) && (max.change > convergence)) && (iter <= max.iter) && (expectation == FALSE) && (maximization == FALSE))    {   

# E-step
### cat("E-step\n") 
#_________________________________________________________________________________________________
#
        BHLK.results	<- try(BHLK.filter(u=u, PR_TR=PR_TR, Sigma=Sigma, xi.bar=xi.bar)) # Filtering
        if(inherits(BHLK.results, "try-error")){
            initialization <- TRUE
        }else{   
            fp 	    <- BHLK.results$fp
            sp	    <- BHLK.smoother(fp, PR_TR) # Smoothing 


# M-step
### cat("M-step\n")
#_________________________________________________________________________________________________
#
# Transition matrix from discrete space of smoothed probabilities
            out.PR_TR   <- sp.to.PR_TR(fp, sp, PR_TR)
            PR_TR       <- out.PR_TR$PR_TR
            xi.bar      <- Ergodic.PR_TR(PR_TR)


# Regression step to get the regime-specific estimates
            max.results <- try(M.step.MSIAH(y=Y, p=p, M=M, fp=fp, sp=sp, PR_TR=PR_TR))
            # Error handling
            if(inherits(max.results, "try-error")){
                    maximization <- TRUE
            }else{
                Beta    <- max.results$Beta
                Sigma   <- max.results$Sigma
                u       <- max.results$u


				
# Convergence criteria
#_____________________________________________________________________________________________
#
				likelihood.val  <- cbind(likelihood.val, max.results$likelihood)
                out.qps         <- QPS.mapping(Sp=sp, Xi=TPV.Xi) 
                qps.val         <- cbind(qps.val, out.qps$qps)
                regime1         <- out.qps$regime1


				if(likelihood.val[iter+1] <= 0){
                    convergence.criterion	<- -100 * (likelihood.val[iter+1] - likelihood.val[iter])/(likelihood.val[iter])
                }else{
                    convergence.criterion	<- 100 * (likelihood.val[iter+1] - likelihood.val[iter])/(likelihood.val[iter])
                }

				# Change in parameters
                parametersm1<- parameters
                parameters  <- rbind(matrix(PR_TR,ncol=1),matrix(Beta,ncol=1),matrix(Sigma,ncol=1))
                change      <- parameters - parametersm1
                max.change  <- abs(max(change))   

# Convergence Infos
#_____________________________________________________________________________________________
#               
                ### cat("Iteration:", iter, "\n\tLikelihood:", likelihood.val[iter+1], "\n\t% change:", convergence.criterion, "\n\tMax.change:", max.change, "\n")
              
                op <- par(mfrow = c(M+2,1), mar=c(2,1,3,0.5), oma=c(1.5,2,1,1)) # M+2 pictures on one plot, and margins settings
                # Plot of likelihood over iteration
                plot(seq(1:(length(likelihood.val)-1)), likelihood.val[-1], type="l", main="Likelihood climbing")
                # Plot of qps over iterations
                plot(seq(1:(length(qps.val)-1)), qps.val[-1], type="l", main="qps descent")
               

                # Plot of smoothed probabilities
                for(regime in 1:M){
                    sp.plot <- ts(t(sp), start=c(1, 1+p), frequency=1)
                    title   <- paste("Smoothed probabilities, regime ", regime, sep="")
                    plot(sp.plot[,regime], main=title ,ylim=c(0, 1))
                    lines(TPV.Xi - 1, type="l", col = "red", lty=3)
                }
                # par(ask=TRUE)
                par(op) # At end of plotting, reset to previous settings 
               
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
       	if(likelihood.val[length(likelihood.val)]  <= 0){
            likelihood.climbing	<- -100 * (likelihood.val[length(likelihood.val)] - likelihood.val[2])/likelihood.val[2] 
        }else{
            likelihood.climbing	<- 100 * (length(likelihood.val) - likelihood.val[2])/(likelihood.val[2])
        }  
        cat("Overall % likelihood change between first and last iteration:", likelihood.climbing, "\n")
        cat("\t\tLast likelihood value:", likelihood.val[length(likelihood.val)], "\n")
 

        # Information Criteria
        IC  <- get.IC(model.type="MSIAH",K=K,p=p,M=M,T=TT-p,max.log.likelihood=likelihood.val[length(likelihood.val)],print.results=FALSE)
 

		# Regimes vector from smoothed probabilities
		Xi	<- matrix(max.col(t(sp)), ncol=1)

		
		# Set up the output
		output <- list(	Sigma               = Sigma,
						A                   = Beta,
						PR_TR               = PR_TR,
						fp                  = t(fp),
						sp                  = t(sp),
						Xi                  = Xi,
						Y                   = Y,
						U                   = u,
                        likelihood          = likelihood.val,
                        likelihood.climbing = likelihood.climbing,
                        IC                  = IC,
						convergence.code    = convergence.code)
	}else{
		output <- list( Y                   = Y,
						convergence.code    = convergence.code,
                        likelihood.climbing = 0)
	}
    return(output)
}
