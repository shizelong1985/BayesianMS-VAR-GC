
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.




#---------------------------------------------------------------------------------------------------
# VAR Gibbs: regression step 
# Independent Normal-Wishart priors form Koop-Korobilis (2010)
#---------------------------------------------------------------------------------------------------
    
G.regression.BVAR.INWP <- function(aux, priors, restrictions){

    # Setup constants 
    #-----------------------------------------------------------------------------
    TT 	<- dim(aux$U)[1]
    K 	<- dim(aux$U)[2] 
    p   <- (dim(aux$Beta)[1] - 1) / K
    T   <- TT + p


    # X, Y
    #-----------------------------------------------------------------------------
 	if(p == 0){
        X   <- matrix(data=1, nrow=T, ncol=1)	# A column of ones for the intercept    
        Y <- matrix(aux$Y,ncol=K)     
    }else{
        X   <- matrix(data=1, nrow=TT, ncol=1)	# A column of ones for the intercept    
        for(lag in 1:p){
            Y.m.j   <- aux$Y[(p+1-lag):(T-lag),]
            X       <- cbind(X,Y.m.j)
        }
        Y <- matrix(aux$Y[(p+1):T,],ncol=K)
    }
     
    # Beta
    #----------------------------------------------------------------------------- 
    # Inverse of var
    Sm1 <- solve(aux$Sigma)

    # Crossproducts
    Zp.Sm1.Z <- matrix(0, K*(1+K*p), K*(1+K*p))
    Zp.Sm1.y <- matrix(0, K*(1+K*p), 1)
    for(t in 1:TT){
        # Z.t matrix, including restrictions
        Z.t <- matrix(0, K, K*(1+K*p))
        for(k in 1:K){
            col.b   <- (1+K*p) * (k-1) + 1  
            col.e   <- col.b + (1+K*p) - 1
            Z.t[k,col.b:col.e]  <- matrix(X[t,],nrow=1) * restrictions$Beta[,k]
        }  
        Zp.Sm1      <- crossprod(Z.t,Sm1)
        Zp.Sm1.Z    <- Zp.Sm1.Z + Zp.Sm1 %*% Z.t
        Zp.Sm1.y    <- Zp.Sm1.y + Zp.Sm1 %*% matrix(Y[t,],nrow=K)
    }

    # Variance of Beta
    V.Beta.bar  <- solve(priors$V.Beta.bar.m1 + Zp.Sm1.Z)

    # Mean of Beta
    Beta.bar    <- V.Beta.bar %*% ( priors$V.Beta.bar.m1 %*% priors$Beta.bar + Zp.Sm1.y )

    # Draw from multivariate normal (mvtnorm package)
    draw        <- rmvnorm(n=1, mean=Beta.bar, sigma=V.Beta.bar, method="chol")
    
    # Store
    #-----------------------------------------------------------------------------
    # Beta (impose restrictions again, due to priors)
    aux$Beta    <- matrix(draw, ncol=K) * restrictions$Beta

    # Residuals: from drawn Beta
    res.tmp <- Y - X %*% aux$Beta
    aux$U   <- matrix(res.tmp, ncol=K)

    
    # Output 
    #-----------------------------------------------------------------------------
    return(aux)
} 
