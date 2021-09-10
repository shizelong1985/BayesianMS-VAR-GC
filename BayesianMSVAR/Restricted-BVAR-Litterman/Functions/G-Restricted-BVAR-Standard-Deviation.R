
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
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   

    
    
    # Restrictions
    #-----------------------------------------------------------------------------
    restrictions.SD    <- diag(restrictions$Sigma)  
    
    # 0: Residuals
    #----------------------------------------
    U   <- matrix(aux$U, ncol=K)    
    
    # 1: Sample moments
    #----------------------------------------
    S   <- apply(U, 2, sd)    # Standard deviations
    V   <- 1/TT * ( TT - 1 - 2 * exp(2*(lgamma(TT/2)-lgamma((TT-1)/2)))) * S    # Variances of S (http://mathworld.wolfram.com/StandardDeviationDistribution.html)
    
    G.Si<- array(NA, aux$G)

    for(i in 1:length(restrictions.SD)){
        # Check for regime independence of the parameter restriction
        if(restrictions.SD[i] != 0){

            # 2: Initial grid for Si
            #----------------------------------------
            lower.bound <- S[i] - 2 * sqrt(V[i])
            upper.bound <- S[i] + 2 * sqrt(V[i])
            S.grid      <- seq(from=lower.bound, to=upper.bound, length=G)

            # 3: Evaluate kernel at the grid points
            #----------------------------------------
            for(j in 1:aux$G){
                # Build covariance matrices at grid point
                Sigma.grid  <- aux$Sigma
                # Likelihood
                Si.grid     <- aux$S
                Si.grid[i]  <- S.grid[j]
                Sigma.grid  <- diag(Si.grid) %*% aux$R %*% diag(Si.grid)  
                # Densities:
                eta.tmp <- dmvnorm(U, sigma=matrix(Sigma.grid,ncol=K)) 
                eta.t   <- matrix(eta.tmp, nrow=1)
                log.likelihood <- sum(log(eta.t)) 

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

            # 6: Store as an auxiliary variable
            #----------------------------------------   
            aux$S[i]   <- Si.new
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
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   

    
    
    # Restritions
    #-----------------------------------------------------------------------------
    restrictions.SD    <- diag(restrictions$Sigma)  
    
    # 0: Residuals
    #----------------------------------------
    U   <- matrix(aux$U, ncol=K)    
    
    # 1: Sample moments
    #----------------------------------------
    S   <- sd(U)    # Standard deviations
    V   <- 1/TT * ( TT - 1 - 2 * exp(2*(lgamma(TT/2)-lgamma((TT-1)/2)))) * S    # Variances of S (http://mathworld.wolfram.com/StandardDeviationDistribution.html)
    
    G.Si<- array(NA, aux$G)

    for(i in 1:length(restrictions.SD)){
        # Check for regime independence of the parameter restriction
        if(restrictions.SD[i] != 0){

            # 2: Initial grid for Si
            #----------------------------------------
            S.grid  <- G.grid$S.grid[i,]       

            # 3: Evaluate kernel at the grid points
            #----------------------------------------
            for(j in 1:aux$G){
                # Build covariance matrices at grid point
                Sigma.grid  <- aux$Sigma
                # Likelihood
                Si.grid     <- aux$S
                Si.grid[i]  <- S.grid[j]
                Sigma.grid  <- diag(Si.grid) %*% aux$R %*% diag(Si.grid)  
                # Densities:
                eta.tmp <- dmvnorm(U, sigma=matrix(Sigma.grid,ncol=K)) 
                eta.t   <- matrix(eta.tmp, nrow=1)
                log.likelihood <- sum(log(eta.t)) 

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

            # 6: Store as an auxiliary variable
            #----------------------------------------   
            aux$S[i]    <- Si.new
        }
    }                                

    # Output 
    #-----------------------------------------------------------------------------
    return(aux) 
}   
