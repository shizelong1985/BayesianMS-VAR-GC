
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.




#---------------------------------------------------------------------------------------------------
# MSIAH Gibbs: regression step for Beta
# Independent Normal-Wishart priors form Koop-Korobilis (2010)
#---------------------------------------------------------------------------------------------------
    
G.regression.Beta.INWP <- function(aux, priors, restrictions){

    # Setup constants 
    #-----------------------------------------------------------------------------
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   

    # X, Y
    #-----------------------------------------------------------------------------
    X   <- matrix(data=1, nrow=TT, ncol=1)	# A column of ones for the intercept
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


    # Draw Beta
    #-----------------------------------------------------------------------------
    # Beta as stacked parameters for: regime-stable, regime-1, regime-2, ..., regime-M
    # See Fruhwirth-Schnatter (2006), p.257.

    ### Beta    <- array(NA, c((1+K*p)*(1+M),K))
    ### # Regime-stable
    ### Beta[1:(1+K*p),]    <- aux$Beta.0
    ### # Regime-1, regime-2, ..., regime-M
    ### for(regime in 1:M){
        ### row.b   <- 1+ (1+K*p) +  (1+K*p) * (regime - 1)
        ### row.e   <- row.b + (1+K*p) - 1
        ### Beta[row.b:row.e,]  <- aux$Beta.St[,,regime]
    ### }


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
            row.b   <- 1+ (1+K*p) +  (1+K*p) * (regime - 1)
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

    # Variance of Beta
    prior.V.Beta.bar.m1 <- kronecker(priors$V.Beta.bar.m1, diag(1+M))
    V.Beta.bar          <- solve(prior.V.Beta.bar.m1 + Zp.Sm1.Z)
    
    # Mean of Beta
    prior.Beta.bar  <- kronecker(priors$Beta.bar, array(1,c(1+M,1)))
    Beta.bar        <- V.Beta.bar %*% ( prior.V.Beta.bar.m1 %*% prior.Beta.bar + Zp.Sm1.y )

    # Draw from multivariate normal (mvtnorm package)
    draw        <- matrix(rmvnorm(n=1, mean=Beta.bar, sigma=V.Beta.bar, method="chol"), ncol=K)


    # Store
    #----------------------------------------------------------------------------- 
    # Beta.0 (impose restrictions again, due to priors)
    aux$Beta.0  <- matrix(draw[1:(1+K*p),],ncol=K) * restrictions$Beta.0

    # Beta.St (impose restrictions again, due to priors)
    for(regime in 1:M){
        row.b   <- 1+ (1+K*p)*regime 
        row.e   <- row.b + (1+K*p) - 1
        aux$Beta.St[,,regime] <- matrix(draw[row.b:row.e,],ncol=K) * restrictions$Beta.St[,,regime]
    }

    # Beta as sum of regime invariant and regime dependent draws
    for(regime in 1:M){
        aux$Beta[,,regime]  <- aux$Beta.0 + aux$Beta.St[,,regime] 
    }

    # Residuals
    #-----------------------------------------------------------------------------
    for(regime in 1:M){
        aux$U[,,regime]  <- matrix(Y - X %*% aux$Beta[,,regime], ncol=K)
    }

    
    # Output 
    #-----------------------------------------------------------------------------
    return(aux)
} 
