
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# BVAR Gibbs: forecast  
#---------------------------------------------------------------------------------------------------
    
G.forecast.BVAR  <- function(aux, posteriors, h, n){

    # Setup constants 
    #-----------------------------------------------------------------------------
    Y   <- as.matrix(aux$Y)
    TT 	<- dim(aux$U)[1]
    K 	<- dim(aux$U)[2] 
    p   <- (dim(aux$Beta)[1] - 1) / K
    T   <- TT + p

    # Matrix of last p values of Y and h forecasted values of Y
    YY      <- array(NA, c(p+h, K))
    if(p>0) YY[1:p,]<- Y[(T-p+1):T,]
    
    for(horizon in 1:h){
        # Next period's mean and variance
        #-------------------------------------------------------------------------  
        Beta            <- aux$Beta
        Y.index         <- p + horizon
        
        X.bar           <- array(1, c(1,1))
        if(p>0) for(lag in 1:p) X.bar <- cbind(X.bar, matrix(YY[Y.index-lag,],ncol=K))
        
        YY[Y.index,]    <- X.bar %*% Beta

        var.YY          <- aux$Sigma
        
        # Draw from multivariate normal
        #-------------------------------------------------------------------------
        par.n       <- 1
        par.mean    <- as.vector(YY[Y.index,])
        par.sigma   <- var.YY

        par.lower   <- rep(-Inf, length = length(par.mean))
        par.upper   <- rep( Inf, length = length(par.mean)) 
        YY[Y.index,]<- rtmvnorm(n=par.n, mean=par.mean, sigma=par.sigma, lower=par.lower, upper=par.upper,algorithm="rejection")  

        ### YY[Y.index,]<- rmvnorm(n=par.n, mean=par.mean, sigma=par.sigma) # Use drawn data as basis for next horizon

    }
    if(p==0){
        posteriors$F.Y[,,n] <- YY
    }else{
        posteriors$F.Y[,,n] <- YY[-(1:p),]
    }


    return(posteriors)
}


