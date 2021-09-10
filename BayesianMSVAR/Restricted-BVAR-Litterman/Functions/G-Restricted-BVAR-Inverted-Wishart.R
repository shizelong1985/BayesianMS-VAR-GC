
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Gibbs: inverted Wishart step for BVAR models
# - Decomposition into standard deviation + covariance
# - Griddy Gibbs
#---------------------------------------------------------------------------------------------------

G.inverted.wishart.BVAR.Griddy.Gibbs  <- function(aux, priors){
    
    # Setup constants 
    #-----------------------------------------------------------------------------
    TT 	<- dim(aux$U)[1]
    K 	<- dim(aux$U)[2] 
    T   <- TT + p

# Standard deviations vector 
#----------------------------------------------------------------------------- 
    
    # 1: Sample moments
    #----------------------------------------
    S   <- sd(aux$U)    # Standard deviations
    V   <- 1/TT * ( TT - 1 - 2 * exp(2*(lgamma(TT/2)-lgamma((TT-1)/2)))) * S    # Variances of S (http://mathworld.wolfram.com/StandardDeviationDistribution.html)
    

    for(i in 1:K){
        
        G.Si<- array(NA, aux$G) 

        # 2: Initial grid for Si
        #----------------------------------------
        lower.bound <- S[i] - 3 * sqrt(V[i])
        upper.bound <- S[i] + 3 * sqrt(V[i])
        S.grid      <- seq(from=lower.bound, to=upper.bound, length=G)

        # 3: Evaluate kernel at the grid points
        #----------------------------------------
        for(j in 1:aux$G){
            # Build covariance matrix at grid point
            Si.grid     <- aux$S
            Si.grid[i]  <- S.grid[j]
            Sigma.grid  <- diag(Si.grid) %*% aux$R %*% diag(Si.grid)

            # Likelihood
            eta.tmp <- dmvnorm(matrix(aux$U,ncol=K), sigma=matrix(Sigma.grid,ncol=K), log=TRUE) 
            eta.t   <- matrix(eta.tmp, nrow=1)
            log.likelihood  <- sum(eta.t)

            # Prior for Si, log normal
            log.Si.prior    <- dlnorm(S.grid[j], meanlog=priors$S.mean, sdlog=priors$S.sd, log=TRUE)

            # Store in G.Si
            G.Si[j] <- log.likelihood + log.Si.prior
        }
        ### plot(x=S.grid, y=G.Si)

        # Rescale with constant (for exponential to compute) and transform with exponential
        G.Si    <- G.Si + abs(max(G.Si))
        G.Si    <- exp(G.Si)
        ### plot(x=S.grid, y=G.Si)

        # 4: Deterministic integration and normalization
        #----------------------------------------
        G.Phi <- array(0, aux$G)
        # Integration
        for(j in 2:G){
            G.Phi[j]    <- G.Phi[j-1] + (S.grid[j]-S.grid[j-1]) * ( abs(G.Si[j-1]) + abs(G.Si[j]-G.Si[j-1])/2)
        }
        ### plot(x=S.grid, y=G.Phi)
        # Normalization to 1
        G.Phi   <- G.Phi / G.Phi[aux$G]
        ### plot(x=S.grid, y=G.Phi)
        
        # 5: Random draw of Si by numerical interpolation
        #----------------------------------------
        draw    <- runif(n=1, min=0, max=1)

        # Interpolation
        index.min   <- max(which(G.Phi<=draw))
        index.max   <- min(which(G.Phi>=draw))
        Si.new      <- S.grid[index.min] + ( S.grid[index.max]-S.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

        # Store as an auxiliary variable
        aux$S[i]    <- Si.new
    }


# Covariance matrix
#-----------------------------------------------------------------------------
    R.length    <- length(aux$R[upper.tri(aux$R, diag=FALSE)])

    for(i in 1:R.length){

        G.Ri        <- array(NA, aux$G)
        R.elements  <- aux$R[upper.tri(aux$R, diag=FALSE)]

        # 6: Calculate support for correlation
        #----------------------------------------
        # f(0)
        R.0 <- aux$R
        R.elements.0    <- R.elements
        R.elements.0[i] <- 0
        R.0[upper.tri(R.0, diag=FALSE)] <- R.elements.0
        R.0[lower.tri(R.0, diag=FALSE)] <- R.elements.0
        f.0 <- det(R.0)

        # f(-1)
        R.m1    <- aux$R
        R.elements.m1   <- R.elements
        R.elements.m1[i]<- -1
        R.m1[upper.tri(R.m1, diag=FALSE)] <- R.elements.m1
        R.m1[lower.tri(R.m1, diag=FALSE)] <- R.elements.m1
        f.m1    <- det(R.m1)

        # f(1)
        R.1     <- aux$R
        R.elements.1    <- R.elements
        R.elements.1[i] <- 1
        R.1[upper.tri(R.1, diag=FALSE)] <- R.elements.1
        R.1[lower.tri(R.1, diag=FALSE)] <- R.elements.1
        f.1    <- det(R.1)

        # k, l, m
        k   <- (f.1 + f.m1 - 2*f.0) / 2
        l   <- (f.1 - f.m1) / 2
        m   <- f.0

        # Roots of f(r)=kr^2 + lr + m
        a   <- (-l - sqrt(l^2 - 4*k*m)) / (2*k)
        b   <- (-l + sqrt(l^2 - 4*k*m)) / (2*k)

        # 7: Grid for Ri
        #---------------------------------------- 
        lower.bound <- min(a,b)
        upper.bound <- max(a,b)
        R.grid      <- seq(from=lower.bound+0.000001, to=upper.bound-0.000001, length=G)

        # 8: Evaluate kernel at the grid points
        #----------------------------------------
        for(j in 1:aux$G){
            # Build covariance matrix at grid point
            Ri.grid <- aux$R
            R.elements.i    <- R.elements
            R.elements.i[i] <- R.grid[j] 
            Ri.grid[upper.tri(Ri.grid, diag=FALSE)] <- R.elements.i
            Ri.grid[lower.tri(Ri.grid, diag=FALSE)] <- R.elements.i 
            Sigma.grid  <- diag(aux$S) %*% Ri.grid %*% diag(aux$S)

            # Likelihood
            eta.tmp <- dmvnorm(matrix(aux$U,ncol=K), sigma=matrix(Sigma.grid,ncol=K), log=TRUE) 
            eta.t   <- matrix(eta.tmp, nrow=1)
            log.likelihood  <- sum(eta.t)

            # Prior for Si, log normal
            log.Ri.prior    <- dunif(R.grid[j], min=lower.bound, max=upper.bound, log=TRUE)

            # Store in G.Ri
            G.Ri[j] <- log.likelihood + log.Ri.prior
        }
        ### plot(x=R.grid, y=G.Ri)    
                
        # Rescale with constant (for exponential to compute) and transform with exponential
        G.Ri    <- G.Ri + abs(max(G.Ri))
        G.Ri    <- exp(G.Ri)
        ### plot(x=R.grid, y=G.Ri)


        # 9: Deterministic integration and normalization
        #----------------------------------------
        G.Phi <- array(0, aux$G)
        # Integration
        for(j in 2:G){
            G.Phi[j]    <- G.Phi[j-1] + (R.grid[j]-R.grid[j-1]) * ( abs(G.Ri[j-1]) + abs(G.Ri[j]-G.Ri[j-1])/2)
        }
        ### plot(x=R.grid, y=G.Phi)
        # Normalization to 1
        G.Phi   <- G.Phi / G.Phi[aux$G]
        ### plot(x=R.grid, y=G.Phi)
        
        # 10: Random draw of Si by numerical interpolation
        #----------------------------------------
        draw    <- runif(n=1, min=0, max=1)

        # Interpolation
        index.min   <- max(which(G.Phi<=draw))
        index.max   <- min(which(G.Phi>=draw))
        Ri.new      <- R.grid[index.min] + ( R.grid[index.max]-R.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

        # Store as an auxiliary variable
        R.new           <- aux$R
        R.elements.i    <- R.elements
        R.elements.i[i] <- Ri.new
        R.new[upper.tri(R.new, diag=FALSE)] <- R.elements.i
        R.new[lower.tri(R.new, diag=FALSE)] <- R.elements.i 
        aux$R    <- R.new
    }



    # 10: Update auxiliary Sigma matrix
    #----------------------------------------   
    aux$Sigma  <- diag(aux$S) %*% aux$R %*% diag(aux$S)


  

    # Output 
    #-----------------------------------------------------------------------------
    return(aux)
}











#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
 









#---------------------------------------------------------------------------------------------------
# Gibbs: inverted Wishart step for BVAR models
# - Decomposition into standard deviation + covariance
# - Griddy Gibbs
#---------------------------------------------------------------------------------------------------

G.inverted.wishart.BVAR.Adaptive.Griddy.Gibbs  <- function(aux, priors, G.grid){
    
    # Setup constants 
    #-----------------------------------------------------------------------------
    TT 	<- dim(aux$U)[1]
    K 	<- dim(aux$U)[2] 
    T   <- TT + p

# Standard deviations vector 
#----------------------------------------------------------------------------- 
    
    for(i in 1:K){
        
        G.Si<- array(NA, aux$G) 

        # 2: Grid for Si
        #----------------------------------------
        S.grid      <- G.grid$S[i,]

        # 3: Evaluate kernel at the grid points
        #----------------------------------------
        for(j in 1:aux$G){
            # Build covariance matrix at grid point
            Si.grid     <- aux$S
            Si.grid[i]  <- S.grid[j]
            Sigma.grid  <- diag(Si.grid) %*% aux$R %*% diag(Si.grid)

            # Likelihood
            eta.tmp <- dmvnorm(matrix(aux$U,ncol=K), sigma=matrix(Sigma.grid,ncol=K), log=TRUE) 
            eta.t   <- matrix(eta.tmp, nrow=1)
            log.likelihood  <- sum(eta.t)

            # Prior for Si, log normal
            log.Si.prior    <- dlnorm(S.grid[j], meanlog=priors$S.mean, sdlog=priors$S.sd, log=TRUE)

            # Store in G.Si
            G.Si[j] <- log.likelihood + log.Si.prior
        }
        ### plot(x=S.grid, y=G.Si)

        # Rescale with constant (for exponential to compute) and transform with exponential
        G.Si    <- G.Si + abs(max(G.Si))
        G.Si    <- exp(G.Si)
        ### plot(x=S.grid, y=G.Si)

        # 4: Deterministic integration and normalization
        #----------------------------------------
        G.Phi <- array(0, aux$G)
        # Integration
        for(j in 2:G){
            G.Phi[j]    <- G.Phi[j-1] + (S.grid[j]-S.grid[j-1]) * ( abs(G.Si[j-1]) + abs(G.Si[j]-G.Si[j-1])/2)
        }
        ### plot(x=S.grid, y=G.Phi)
        # Normalization to 1
        G.Phi   <- G.Phi / G.Phi[aux$G]
        ### plot(x=S.grid, y=G.Phi)
        
        # 5: Random draw of Si by numerical interpolation
        #----------------------------------------
        draw    <- runif(n=1, min=0, max=1)

        # Interpolation
        index.min   <- max(which(G.Phi<=draw))
        index.max   <- min(which(G.Phi>=draw))
        Si.new      <- S.grid[index.min] + ( S.grid[index.max]-S.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

        # Store as an auxiliary variable
        aux$S[i]    <- Si.new
    }


# Covariance matrix
#-----------------------------------------------------------------------------
    
    R.length    <- length(aux$R[upper.tri(aux$R, diag=FALSE)])

    for(i in 1:R.length){

        G.Ri        <- array(NA, aux$G)    
        R.elements  <- aux$R[upper.tri(aux$R, diag=FALSE)] 

        # 6: Calculate support for correlation
        #----------------------------------------
        # f(0)
        R.0 <- aux$R
        R.elements.0    <- R.elements
        R.elements.0[i] <- 0
        R.0[upper.tri(R.0, diag=FALSE)] <- R.elements.0
        R.0[lower.tri(R.0, diag=FALSE)] <- R.elements.0
        f.0 <- det(R.0)

        # f(-1)
        R.m1    <- aux$R
        R.elements.m1   <- R.elements
        R.elements.m1[i]<- -1
        R.m1[upper.tri(R.m1, diag=FALSE)] <- R.elements.m1
        R.m1[lower.tri(R.m1, diag=FALSE)] <- R.elements.m1
        f.m1    <- det(R.m1)

        # f(1)
        R.1     <- aux$R
        R.elements.1    <- R.elements
        R.elements.1[i] <- 1
        R.1[upper.tri(R.1, diag=FALSE)] <- R.elements.1
        R.1[lower.tri(R.1, diag=FALSE)] <- R.elements.1
        f.1    <- det(R.1)

        # k, l, m
        k   <- (f.1 + f.m1 - 2*f.0) / 2
        l   <- (f.1 - f.m1) / 2
        m   <- f.0

        # Roots of f(r)=kr^2 + lr + m
        a   <- (-l - sqrt(l^2 - 4*k*m)) / (2*k)
        b   <- (-l + sqrt(l^2 - 4*k*m)) / (2*k)

        # 7: Grid for Ri
        #---------------------------------------- 
        lower.bound <- min(a,b)
        upper.bound <- max(a,b)
        R.grid      <- G.grid$R[i,]

        # 8: Evaluate kernel at the grid points
        #----------------------------------------
        for(j in 1:aux$G){
            # Build covariance matrix at grid point
            Ri.grid <- aux$R
            R.elements.i    <- R.elements
            R.elements.i[i] <- R.grid[j] 
            Ri.grid[upper.tri(Ri.grid, diag=FALSE)] <- R.elements.i
            Ri.grid[lower.tri(Ri.grid, diag=FALSE)] <- R.elements.i 
            Sigma.grid  <- diag(aux$S) %*% Ri.grid %*% diag(aux$S)

            # Likelihood
            eta.tmp <- dmvnorm(matrix(aux$U,ncol=K), sigma=matrix(Sigma.grid,ncol=K), log=TRUE) 
            eta.t   <- matrix(eta.tmp, nrow=1)
            log.likelihood  <- sum(eta.t)

            # Prior for Si, log normal
            log.Ri.prior    <- dunif(R.grid[j], min=lower.bound, max=upper.bound, log=TRUE)

            # Store in G.Ri
            G.Ri[j] <- log.likelihood + log.Ri.prior
        }
        ### plot(x=R.grid, y=G.Ri)    
                
        # Rescale with constant (for exponential to compute) and transform with exponential
        G.Ri    <- G.Ri + abs(max(G.Ri))
        G.Ri    <- exp(G.Ri)
        ### plot(x=R.grid, y=G.Ri)


        # 9: Deterministic integration and normalization
        #----------------------------------------
        G.Phi <- array(0, aux$G)
        # Integration
        for(j in 2:G){
            G.Phi[j]    <- G.Phi[j-1] + (R.grid[j]-R.grid[j-1]) * ( abs(G.Ri[j-1]) + abs(G.Ri[j]-G.Ri[j-1])/2)
        }
        ### plot(x=R.grid, y=G.Phi)
        # Normalization to 1
        G.Phi   <- G.Phi / G.Phi[aux$G]
        ### plot(x=R.grid, y=G.Phi)
        
        # 10: Random draw of Si by numerical interpolation
        #----------------------------------------
        draw    <- runif(n=1, min=0, max=1)

        # Interpolation
        index.min   <- max(which(G.Phi<=draw))
        index.max   <- min(which(G.Phi>=draw))
        Ri.new      <- R.grid[index.min] + ( R.grid[index.max]-R.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

        # Store as an auxiliary variable
        R.new           <- aux$R
        R.elements.i    <- R.elements
        R.elements.i[i] <- Ri.new
        R.new[upper.tri(R.new, diag=FALSE)] <- R.elements.i
        R.new[lower.tri(R.new, diag=FALSE)] <- R.elements.i 
        aux$R    <- R.new
    }



    # 10: Update auxiliary Sigma matrix
    #----------------------------------------   
    aux$Sigma  <- diag(aux$S) %*% aux$R %*% diag(aux$S)


  

    # Output 
    #-----------------------------------------------------------------------------
    return(aux)
}  
