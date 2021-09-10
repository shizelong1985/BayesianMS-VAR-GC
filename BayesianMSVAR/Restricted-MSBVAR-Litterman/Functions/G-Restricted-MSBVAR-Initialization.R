
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Gibbs: Function converting EM output for being used in the Gibbs
#---------------------------------------------------------------------------------------------------

G.EM.to.Gibbs <- function(model.type, EM.output, priors, restrictions, G){

    M       <- dim(EM.output$PR_TR)[1]
    K       <- dim(EM.output$U)[2]
    T       <- dim(EM.output$Y)[1]
    TT      <- dim(EM.output$U)[1] 
    p       <- T - TT


    if(model.type == "MSIAH") {
        
        Beta.0  <- matrix(matrix(EM.output$Beta[,,1],ncol=K) * restrictions$Beta.0, ncol=K)
        Beta.St <- array(NA,c(1+K*p, K, M))
        for(regime in 1:M){
            Beta.St[,,regime] <- matrix(matrix(EM.output$Beta[,,regime],ncol=K) * restrictions$Beta.St[,,regime], ncol=K)
        }

        # Decomposition of covariance into standard deviation and correlation
        #-----------------------------------------------------------------------------
        S   <- array(NA, c(K,M))
        R   <- array(NA, c(K,K,M))

        for(regime in 1:M){
            xi.m    <- max.col(EM.output$sp)
            T.m     <- sum(xi.m)

            # Consider residuals when xi.m = 1
            U.m <- matrix(EM.output$U[which(xi.m==1),,regime],ncol=K)

            S[,regime]  <- apply(U.m, 2, sd)
            R[,,regime] <- cor(U.m)

        }
            

        Gibbs.input <- list(
            Y       = EM.output$Y,
            U       = EM.output$U,
            xi      = array(NA, c(M,TT)),
            Beta    = EM.output$Beta,
            Beta.0  = Beta.0,
            Beta.St = Beta.St,
            Sigma   = EM.output$Sigma,
            PR_TR   = EM.output$PR_TR,
			w		= priors$w / sum(priors$w),
            ### likelihood = EM.output$likelihood[length(EM.output$likelihood)] - 100000
            S       = S,
            R       = R,
            G       = G 
            ) 


    }else{
        stop("\n\n\t The model entered as parameter is not implemented.")
    }  


    return(Gibbs.input)
}




#---------------------------------------------------------------------------------------------------
# Gibbs: Function Generating restrictions
#---------------------------------------------------------------------------------------------------

G.set.restrictions <- function(model.type, EM.output){

    M       <- dim(EM.output$PR_TR)[1]
    K       <- dim(EM.output$U)[2]
    T       <- dim(EM.output$Y)[1]
    TT      <- dim(EM.output$U)[1] 
    p       <- T - TT

    Beta.0  <- array(0, c(1+K*p, K)) 
    Beta.St <- array(0, c(1+K*p, K, M))

    
    if(model.type == "MSI") {
        Beta.0[-1,]  <- 1
        for(regime in 1:M){
            Beta.St[1,,regime]   <- 1
        }
        heteroskedasticity <- FALSE
        Sigma.0  <- array(1, c(K, K))
        Sigma.St <- array(0, c(K, K, M))    

    }else if(model.type == "MSH") {
        Beta.0  <- array(1, c(1+K*p, K))
        heteroskedasticity <- TRUE
        Sigma.0  <- array(0, c(K, K))
        Sigma.St <- array(1, c(K, K, M))

    }else if(model.type == "MSA") {
        Beta.0[1,]  <- 1
        for(regime in 1:M){
            Beta.St[-1,,regime]   <- 1
        }
        heteroskedasticity <- FALSE
        Sigma.0  <- array(1, c(K, K))
        Sigma.St <- array(0, c(K, K, M))

    }else if(model.type == "MSIA") {
        for(regime in 1:M){
            Beta.St[,,regime]   <- 1
        }
        heteroskedasticity <- FALSE
        Sigma.0  <- array(1, c(K, K))
        Sigma.St <- array(0, c(K, K, M))

    }else if(model.type == "MSIH") {
        Beta.0[-1,]  <- 1
        for(regime in 1:M){
            Beta.St[1,,regime]   <- 1
        }
        heteroskedasticity <- TRUE
        Sigma.0  <- array(0, c(K, K))
        Sigma.St <- array(1, c(K, K, M))

    }else if(model.type == "MSAH") {
        Beta.0[1,]  <- 1
        for(regime in 1:M){
            Beta.St[-1,,regime]   <- 1
        }
        heteroskedasticity <- TRUE
        Sigma.0  <- array(0, c(K, K))
        Sigma.St <- array(1, c(K, K, M))

    }else if(model.type == "MSIAH") {
        for(regime in 1:M){
            Beta.St[,,regime]   <- 1
        }
        heteroskedasticity <- TRUE
        Sigma.0  <- array(0, c(K, K))
        Sigma.St <- array(1, c(K, K, M))

    }else{
        stop("\n\n\t The model entered as parameter is not implemented.")
    }  


    restrictions    <- list(Beta.0  = Beta.0,
                            Beta.St = Beta.St,
                            ### heteroskedasticity = heteroskedasticity,
                            Sigma.0 = Sigma.0,
                            Sigma.St= Sigma.St
                        )

    return(restrictions)
}


                                    
