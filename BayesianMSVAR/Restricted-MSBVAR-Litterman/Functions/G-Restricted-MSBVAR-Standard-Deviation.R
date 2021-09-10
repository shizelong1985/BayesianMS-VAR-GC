
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Gibbs: standard deviation step
#---------------------------------------------------------------------------------------------------

G.standard.deviation    <- function(aux,priors,restrictions){
    
    # Setup constants 
    #-----------------------------------------------------------------------------
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   

    
    
    # Regime invariant parameters
    #-----------------------------------------------------------------------------
    regime.stable.SD    <- diag(restrictions$Sigma.0)  
    
    # 0: Residuals which occur
    #----------------------------------------
    xi.m    <- aux$xi[1,] 
    U   <-matrix(aux$U[which(xi.m==1),,1],ncol=K)
    for(regime in 2:M){
        xi.m    <- aux$xi[regime,] 
        U   <- rbind(U, matrix(aux$U[which(xi.m==1),,regime],ncol=K))
    }
    U   <- matrix(U, ncol=K)    
    
    # 1: Sample moments
    #----------------------------------------
    S   <- apply(U, 2, sd)  # Standard deviations
    V   <- 1/TT * ( TT - 1 - 2 * exp(2*(lgamma(TT/2)-lgamma((TT-1)/2)))) * S    # Variances of S (http://mathworld.wolfram.com/StandardDeviationDistribution.html)
    
    G.Si<- array(NA, aux$G)

    for(i in 1:length(regime.stable.SD)){
        # Check for regime independence of the parameter restriction
        if(regime.stable.SD[i] != 0){

            # 2: Initial grid for Si
            #----------------------------------------
            lower.bound <- S[i] - 3 * sqrt(V[i])
            upper.bound <- S[i] + 3 * sqrt(V[i])
            S.grid      <- seq(from=lower.bound, to=upper.bound, length=G)

            # 3: Evaluate kernel at the grid points
            #----------------------------------------
            for(j in 1:aux$G){
                # Build covariance matrices at grid point
                Sigma.grid  <- aux$Sigma
                # Likelihood
                eta.t 	<- matrix(NA, M, TT)
                for(regime in 1:M){
                    Si.grid     <- aux$S[,regime]
                    Si.grid[i]  <- S.grid[j]
                    Sigma.grid[,,regime]    <- diag(Si.grid) %*% aux$R[,,regime]%*% diag(Si.grid)  
                    # Densities:
                    eta.tmp         <- dmvnorm(U, sigma=matrix(Sigma.grid[,,regime],ncol=K)) 
                    eta.t[regime,]  <- matrix(eta.tmp, nrow=1)
                }
                
                PR.forecasted    = G.BHLK.forecasting(aux)
                l.likelihood     = colSums(eta.t * PR.forecasted)
                l.likelihood[l.likelihood==-Inf]    = log(10e-30)
                log.likelihood <- sum(l.likelihood)

                # Prior for Si, log normal
                log.Si.prior    <- dlnorm(S.grid[j], meanlog=priors$S.mean, sdlog=priors$S.sd, log=TRUE)
                if (log.Si.prior==-Inf) log.Si.prior = log(10e-40)
                # Store in G.Si
                G.Si[j] <- log.likelihood + log.Si.prior
            }

            # Rescale with constant (for exponential to compute) and transform with exponential
            G.Si    <- G.Si + abs(max(G.Si))
            G.Si    <- exp(G.Si)

            # 4: Deterministic integration and normalization
            #----------------------------------------
            G.Phi <- array(0, aux$G)
            # Integration
            for(j in 2:G){
                G.Phi[j]    <- G.Phi[j-1] + (S.grid[j]-S.grid[j-1]) * ( abs(G.Si[j-1]) + abs(G.Si[j]-G.Si[j-1])/2)
            }
            # Normalization to 1
            G.Phi   <- G.Phi / G.Phi[aux$G]
            
            # 5: Random draw of Si by numerical interpolation
            #----------------------------------------
            draw    <- runif(n=1, min=0, max=1)

            # Interpolation
            index.min   <- max(which(G.Phi<=draw))
            index.max   <- min(which(G.Phi>=draw))
            Si.new      <- S.grid[index.min] + ( S.grid[index.max]-S.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

            # 6: Store as an auxiliary variable, for all regimes
            #----------------------------------------   
            for(regime in 1:M){
                aux$S[i,regime]    <- Si.new
            }
         
        }

    }


    # Regime dependent parameters
    #-----------------------------------------------------------------------------      
    for(regime in 1:M){  

        regime.dependent.SD <- diag(restrictions$Sigma.St[,,regime])

        # 0: Residuals
        #----------------------------------------
        xi.m<- aux$xi[regime,]
        T.m <- sum(xi.m)
        # Consider residuals when xi.m = 1
        U.m <- matrix(aux$U[which(xi.m==1),,regime],ncol=K)
        # 1: Sample moments
        #---------------------------------------- 
        S   <- apply(U.m, 2, sd)    # Standard deviations
        V   <- 1/T.m * ( T.m - 1 - 2 * exp(2*(lgamma(T.m/2)-lgamma((T.m-1)/2)))) * S    # Variances of S (http://mathworld.wolfram.com/StandardDeviationDistribution.html)
        

        for(i in 1:length(regime.dependent.SD)){

            # Check for regime dependence of the parameter restriction
            if(regime.dependent.SD[i] != 0){    
            
                G.Si<- array(NA, aux$G) 

                # 2: Initial grid for Si
                #----------------------------------------
                lower.bound <- S[i] - 3 * sqrt(V[i])
                upper.bound <- S[i] + 3 * sqrt(V[i])
                S.grid      <- seq(from=lower.bound, to=upper.bound, length=G)

                # 3: Evaluate kernel at the grid points
                #----------------------------------------
                for(j in 1:aux$G){
                    # Build covariance matrices at grid point
                    Sigma.grid  <- aux$Sigma
                    Si.grid     <- aux$S[,regime]
                    Si.grid[i]  <- S.grid[j]
                    Sigma.grid[,,regime]    <- diag(Si.grid) %*% aux$R[,,regime]%*% diag(Si.grid)

                    # Likelihood
                    # Evaluate only in actual regime, since other regimes can be treated as constant and here we target the evaluation of a kernel
                    eta.tmp         <- dmvnorm(matrix(U.m,ncol=K), sigma=matrix(Sigma.grid[,,regime],ncol=K)) 
                    eta.tmp[eta.tmp=-Inf]
                    l.likelihood    = log(eta.tmp)
                    l.likelihood[l.likelihood==-Inf] = log(10e-30)
                    log.likelihood <- sum(l.likelihood) 

                    # Prior for Si, log normal
                    log.Si.prior    <- dlnorm(S.grid[j], meanlog=priors$S.mean, sdlog=priors$S.sd, log=TRUE)
                    if (log.Si.prior==-Inf) log.Si.prior = log(10e-40)
                    # Store in G.Si
                    G.Si[j] <- log.likelihood + log.Si.prior
                }

                # Rescale with constant (for exponential to compute) and transform with exponential
                G.Si    <- G.Si + abs(max(G.Si))
                G.Si    <- exp(G.Si)

                # 4: Deterministic integration and normalization
                #----------------------------------------
                G.Phi <- array(0, aux$G)
                # Integration
                for(j in 2:G){
                    G.Phi[j]    <- G.Phi[j-1] + (S.grid[j]-S.grid[j-1]) * ( abs(G.Si[j-1]) + abs(G.Si[j]-G.Si[j-1])/2)
                }
                # Normalization to 1
                G.Phi   <- G.Phi / G.Phi[aux$G]
                
                # 5: Random draw of Si by numerical interpolation
                #----------------------------------------
                draw    <- runif(n=1, min=0, max=1)

                # Interpolation
                index.min   <- max(which(G.Phi<=draw))
                index.max   <- min(which(G.Phi>=draw))
                Si.new      <- S.grid[index.min] + ( S.grid[index.max]-S.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

                # Store as an auxiliary variable
                aux$S[i,regime]    <- Si.new
            }
        }
    }

    # Output 
    #-----------------------------------------------------------------------------
    return(aux)
}








#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------







G.standard.deviation.adaptive.grid  <- function(aux,priors,restrictions,G.grid){
    
    # Setup constants 
    #-----------------------------------------------------------------------------
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   

    
    
    
    # Regime invariant parameters
    #-----------------------------------------------------------------------------
    regime.stable.SD    <- diag(restrictions$Sigma.0)
    
    # 0: Residuals which occur
    #----------------------------------------
    xi.m    <- aux$xi[1,] 
    U   <-matrix(aux$U[which(xi.m==1),,1],ncol=K)
    for(regime in 2:M){
        xi.m    <- aux$xi[regime,] 
        U   <- rbind(U, matrix(aux$U[which(xi.m==1),,regime],ncol=K))
    }
    U   <- matrix(U, ncol=K)    
    
    # 1: Sample moments
    #----------------------------------------
    S   <- apply(U, 2, sd)  # Standard deviations
    V   <- 1/TT * ( TT - 1 - 2 * exp(2*(lgamma(TT/2)-lgamma((TT-1)/2)))) * S    # Variances of S (http://mathworld.wolfram.com/StandardDeviationDistribution.html)
    
    G.Si<- array(NA, aux$G)

    for(i in 1:length(regime.stable.SD)){
        # Check for regime independence of the parameter restriction
        if(regime.stable.SD[i] != 0){

            # 2: Initial grid for Si
            #----------------------------------------
            S.grid  <- G.grid$S.grid[i,,regime]

            # 3: Evaluate kernel at the grid points
            #----------------------------------------
            for(j in 1:aux$G){
                # Build covariance matrices at grid point
                Sigma.grid  <- aux$Sigma
                # Likelihood
                eta.t 	<- matrix(NA, M, TT)
                for(regime in 1:M){
                    Si.grid     <- aux$S[,regime]
                    Si.grid[i]  <- S.grid[j]
                    Sigma.grid[,,regime]    <- diag(Si.grid) %*% aux$R[,,regime]%*% diag(Si.grid)  
                    # Densities:
                    eta.tmp         <- dmvnorm(U, sigma=matrix(Sigma.grid[,,regime],ncol=K)) 
                    eta.t[regime,]  <- matrix(eta.tmp, nrow=1)
                }
                log.likelihood <- sum( log(colSums(eta.t * aux$xi))) 

                # Prior for Si, log normal
                log.Si.prior    <- dlnorm(S.grid[j], meanlog=priors$S.mean, sdlog=priors$S.sd, log=TRUE)
                # Store in G.Si
                G.Si[j] <- log.likelihood + log.Si.prior
            }

            # Rescale with constant (for exponential to compute) and transform with exponential
            G.Si    <- G.Si + abs(max(G.Si))
            G.Si    <- exp(G.Si)

            # 4: Deterministic integration and normalization
            #----------------------------------------
            G.Phi <- array(0, aux$G)
            # Integration
            for(j in 2:G){
                G.Phi[j]    <- G.Phi[j-1] + (S.grid[j]-S.grid[j-1]) * ( abs(G.Si[j-1]) + abs(G.Si[j]-G.Si[j-1])/2)
            }
            # Normalization to 1
            G.Phi   <- G.Phi / G.Phi[aux$G]
            
            # 5: Random draw of Si by numerical interpolation
            #----------------------------------------
            draw    <- runif(n=1, min=0, max=1)

            # Interpolation
            index.min   <- max(which(G.Phi<=draw))
            index.max   <- min(which(G.Phi>=draw))
            Si.new      <- S.grid[index.min] + ( S.grid[index.max]-S.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

            # 6: Store as an auxiliary variable, for all regimes
            #----------------------------------------   
            for(regime in 1:M){
                aux$S[i,regime]    <- Si.new
            }
         
        }

    }


    # Regime dependent parameters
    #-----------------------------------------------------------------------------      
    for(regime in 1:M){  

        regime.dependent.SD <- diag(restrictions$Sigma.St[,,regime])

        # 0: Residuals
        #----------------------------------------
        xi.m<- aux$xi[regime,]
        T.m <- sum(xi.m)
        # Consider residuals when xi.m = 1
        U.m <- matrix(aux$U[which(xi.m==1),,regime],ncol=K)
        # 1: Sample moments
        #---------------------------------------- 
        S   <- apply(U.m, 2, sd)    # Standard deviations
        V   <- 1/T.m * ( T.m - 1 - 2 * exp(2*(lgamma(T.m/2)-lgamma((T.m-1)/2)))) * S    # Variances of S (http://mathworld.wolfram.com/StandardDeviationDistribution.html)
        

        for(i in 1:length(regime.dependent.SD)){

            # Check for regime dependence of the parameter restriction
            if(regime.dependent.SD[i] != 0){    
            
                G.Si<- array(NA, aux$G) 

                # 2: Initial grid for Si
                #----------------------------------------
                S.grid      <- G.grid$S.grid[i,,regime]

                # 3: Evaluate kernel at the grid points
                #----------------------------------------
                for(j in 1:aux$G){
                    # Build covariance matrices at grid point
                    Sigma.grid  <- aux$Sigma
                    Si.grid     <- aux$S[,regime]
                    Si.grid[i]  <- S.grid[j]
                    Sigma.grid[,,regime]    <- diag(Si.grid) %*% aux$R[,,regime]%*% diag(Si.grid)

                    # Likelihood
                    # Evaluate only in actual regime, since other regimes can be treated as constant and here we target the evaluation of a kernel
                    eta.tmp         <- dmvnorm(matrix(U.m,ncol=K), sigma=matrix(Sigma.grid[,,regime],ncol=K)) 
                    log.likelihood <- sum( log(eta.tmp)) 

                    # Prior for Si, log normal
                    log.Si.prior    <- dlnorm(S.grid[j], meanlog=priors$S.mean, sdlog=priors$S.sd, log=TRUE)
                    # Store in G.Si
                    G.Si[j] <- log.likelihood + log.Si.prior
                }

                # Rescale with constant (for exponential to compute) and transform with exponential
                G.Si    <- G.Si + abs(max(G.Si))
                G.Si    <- exp(G.Si)

                # 4: Deterministic integration and normalization
                #----------------------------------------
                G.Phi <- array(0, aux$G)
                # Integration
                for(j in 2:G){
                    G.Phi[j]    <- G.Phi[j-1] + (S.grid[j]-S.grid[j-1]) * ( abs(G.Si[j-1]) + abs(G.Si[j]-G.Si[j-1])/2)
                }
                # Normalization to 1
                G.Phi   <- G.Phi / G.Phi[aux$G]
                
                # 5: Random draw of Si by numerical interpolation
                #----------------------------------------
                draw    <- runif(n=1, min=0, max=1)

                # Interpolation
                index.min   <- max(which(G.Phi<=draw))
                index.max   <- min(which(G.Phi>=draw))
                Si.new      <- S.grid[index.min] + ( S.grid[index.max]-S.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

                # Store as an auxiliary variable
                aux$S[i,regime]    <- Si.new
            }
        }
    }

    # Output 
    #-----------------------------------------------------------------------------
    return(aux) 
}   
