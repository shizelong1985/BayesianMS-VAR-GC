
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# MSVAR Gibbs: forecast  
#---------------------------------------------------------------------------------------------------
    
G.forecast.MSVAR  <- function(aux, posteriors, h, n){

    # Setup constants 
    #-----------------------------------------------------------------------------
    Y   <- as.matrix(aux$Y)
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   

    # Matrix of last p values of Y and h forecasted values of Y
    YY      <- array(NA, c(p+h, K))
    if(p>0) YY[1:p,]<- Y[(T-p+1):T,]  
    
    # Transition probabilities
    F   <- t(aux$PR_TR)

    # Current probabilities of regime
    xi.T.T  <- matrix(aux$xi[,TT], ncol=1)

    for(horizon in 1:h){
        # Draw next period's regime 
        #------------------------------------------------------------------------- 
        xi.Tp1.T        <- F %*% xi.T.T
        xi.draw         <- sample(1:M, 1, prob=xi.Tp1.T)

        # Next period's mean and variance
        #-------------------------------------------------------------------------  
        Beta            <- aux$Beta[,,xi.draw]
        Y.index         <- p + horizon
        
        X.bar           <- array(1, c(1,1))
        if(p>0) for(lag in 1:p) X.bar <- cbind(X.bar, matrix(YY[Y.index-lag,],ncol=K))
        
        YY[Y.index,]    <- X.bar %*% Beta

        var.YY          <- aux$Sigma[,,xi.draw]
        
        # Draw from multivariate normal
        #-------------------------------------------------------------------------
        par.n       <- 1
        par.mean    <- as.vector(YY[Y.index,])
        par.sigma   <- var.YY

        par.lower   <- rep(-Inf, length = length(par.mean))
        par.upper   <- rep( Inf, length = length(par.mean)) 
        YY[Y.index,]<- rtmvnorm(n=par.n, mean=par.mean, sigma=par.sigma, lower=par.lower, upper=par.upper,algorithm="rejection")  

        ### YY[Y.index,]<- rmvnorm(n=par.n, mean=par.mean, sigma=par.sigma) # Use drawn data as basis for next horizon

        # Update posteriors
        #-------------------------------------------------------------------------  
        posteriors$F.S.t[horizon,n] <- xi.draw 
        
        # For next iteration
        #-------------------------------------------------------------------------
        xi.T.T  <- xi.Tp1.T
    }
    
    if(p==0){
        posteriors$F.Y[,,n] <- YY
    }else{
        posteriors$F.Y[,,n] <- YY[-(1:p),]
    }


    return(posteriors)
}





