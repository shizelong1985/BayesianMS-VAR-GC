
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Gibbs: restricted correlation step
# G.correlation(): for griddy Gibbs
# G.correlation.adaptive.grid(): for adaptive griddy Gibbs
#---------------------------------------------------------------------------------------------------

G.correlation   <- function(aux,priors,restrictions){
    
    # Setup constants 
    #-----------------------------------------------------------------------------
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   


    # 1/2 Parameters restricted to 0 
    #-----------------------------------------------------------------------------------------------
    restrictions.R  <- restrictions$Sigma[upper.tri(restrictions$Sigma, diag=FALSE)]
    R.length        <- length(restrictions.R)
    
    R.regime    <- aux$R
    R.elements  <- R.regime[upper.tri(R.regime, diag=FALSE)]
    R           <- aux$R

    for(i in 1:R.length){
        if(restrictions.R[i] == 0){
         R.elements[i]  <- 0   
        }
    }
    
    R       = form.correlation.matrix(R.elements)
    aux$R[,] <- R 


    # 2/2 Regime invariant parameters
    #-----------------------------------------------------------------------------------------------
    restrictions.R  <- restrictions$Sigma[upper.tri(restrictions$Sigma, diag=FALSE)]
    R.length        <- length(restrictions.R) 

    # 0: Residuals
    #---------------------------------------- 
    U   <- matrix(aux$U, ncol=K) 


    # Loop on all correlation parameters
    for(i in 1:R.length){
        # Check for regime independence of the parameter restriction
        if(restrictions.R[i] != 0){

            G.Ri        <- array(NA, aux$G)

            # 1: Calculate support for correlation
            #----------------------------------------
            R.regime    <- aux$R
            R.elements  <- R.regime[upper.tri(R.regime, diag=FALSE)] 
            
            # f(0)
            R.0 <- R.regime
            R.elements.0    <- R.elements
            R.elements.0[i] <- 0
            R.0     = form.correlation.matrix(R.elements.0)
            f.0 <- det(R.0)

            # f(-1)
            R.m1    <- R.regime
            R.elements.m1   <- R.elements
            R.elements.m1[i]<- -1
            R.m1    = form.correlation.matrix(R.elements.m1 )  
            f.m1    <- det(R.m1)

            # f(1)
            R.1     <- R.regime
            R.elements.1    <- R.elements
            R.elements.1[i] <- 1
            R.1     = form.correlation.matrix(R.elements.1)
            f.1    <- det(R.1)

            # k, l, m
            k   <- (f.1 + f.m1 - 2*f.0) / 2
            l   <- (f.1 - f.m1) / 2
            m   <- f.0

            # Roots of f(r)=kr^2 + lr + m
            a   <- (-l - sqrt(l^2 - 4*k*m)) / (2*k)
            b   <- (-l + sqrt(l^2 - 4*k*m)) / (2*k)

            # 2: Grid for Ri: feasible set
            #---------------------------------------- 
            lower.bound <- min(a,b)
            upper.bound <- max(a,b)
            
            R.grid      <- seq(from=lower.bound+0.000001, to=upper.bound-0.000001, length=G)

            # 3: Evaluate kernel at the grid points
            #----------------------------------------
            for(j in 1:aux$G){
                # Build covariance matrix at grid point
                Sigma.grid  <- aux$Sigma
                # Likelihood
                R.regime    <- aux$R
                R.elements  <- R.regime[upper.tri(R.regime, diag=FALSE)]

                Ri.grid         <- aux$R
                R.elements.i    <- R.elements
                R.elements.i[i] <- R.grid[j] 

                Ri.grid     = form.correlation.matrix(R.elements.i)

                Sigma.grid  <- diag(aux$S) %*% Ri.grid %*% diag(aux$S)

                # Densities:
                eta.tmp <- dmvnorm(U, sigma=matrix(Sigma.grid,ncol=K)) 
                eta.t   <- matrix(eta.tmp, nrow=1)   

                log.likelihood <- sum( log(eta.t))

                # Prior for Ri, uniform on the support (feasible set)
                log.Ri.prior    <- dunif(R.grid[j], min=lower.bound, max=upper.bound, log=TRUE) 

                # Store in G.Ri
                G.Ri[j] <- log.likelihood + log.Ri.prior
            }
                    
            # Rescale with constant (for exponential to compute) and transform with exponential
            G.Ri    <- G.Ri + abs(max(G.Ri))
            G.Ri    <- exp(G.Ri)

            # 4: Deterministic integration and normalization
            #----------------------------------------
            G.Phi <- array(0, aux$G)
            # Integration
            for(j in 2:G){
                G.Phi[j]    <- G.Phi[j-1] + (R.grid[j]-R.grid[j-1]) * ( abs(G.Ri[j-1]) + abs(G.Ri[j]-G.Ri[j-1])/2)
            }
            # Normalization to 1
            G.Phi   <- G.Phi / G.Phi[aux$G]
            
            # 5: Random draw of Si by numerical interpolation
            #----------------------------------------
            draw    <- runif(n=1, min=0, max=1)

            # Interpolation
            index.min   <- max(which(G.Phi<=draw))
            index.max   <- min(which(G.Phi>=draw))
            Ri.new      <- R.grid[index.min] + ( R.grid[index.max]-R.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

            # Store as an auxiliary variable
            R.regime    <- aux$R
            R.elements  <- R.regime[upper.tri(R.regime, diag=FALSE)]

            R.elements.i    <- R.elements
            R.elements.i[i] <- Ri.new
            R.regime        = form.correlation.matrix(R.elements.i)
            aux$R       <- R.regime
        }
    }

    # Update auxiliary Sigma matrix
    #----------------------------------------   
    aux$Sigma   <- diag(aux$S) %*% aux$R %*% diag(aux$S)
                                 
    # Output 
    #-----------------------------------------------------------------------------
    return(aux)
}


form.correlation.matrix = function(elements){
  # Creates N x N correlation matrix from unique elements of the correlation matrix 
  # elements    - correlation coefficients obtained by:   corr.matrix[upper.tri(corr.matrix, diag=FALSE)]
  # output      - correlation matrix:                     corr.matrix
  elements  = as.vector(elements)
  N         = .5*(sqrt(8*length(elements) + 1 ) + 1)
  output    = matrix(0, N, N)
  output[upper.tri(output, diag=FALSE)] = elements
  output      = output + t(output) + diag(N)
  return(output)  
}






#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------


G.correlation.adaptive.grid   <- function(aux,priors,restrictions,G.grid){ 
    # Setup constants 
    #-----------------------------------------------------------------------------
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]   


    # 1/2 Parameters restricted to 0 
    #-----------------------------------------------------------------------------------------------
    restrictions.R  <- restrictions$Sigma[upper.tri(restrictions$Sigma, diag=FALSE)]
    R.length        <- length(restrictions.R)
    
    R.regime    <- aux$R
    R.elements  <- R.regime[upper.tri(R.regime, diag=FALSE)]
    R           <- aux$R

    for(i in 1:R.length){
        if(restrictions.R[i] == 0){
         R.elements[i]  <- 0   
        }
    }
    
    R[lower.tri(R, diag=FALSE)]   <- R.elements
    R[upper.tri(R, diag=FALSE)]   <- R.elements 

    aux$R[,] <- R 


    # 2/2 Regime invariant parameters
    #-----------------------------------------------------------------------------------------------
    restrictions.R  <- restrictions$Sigma[upper.tri(restrictions$Sigma, diag=FALSE)]
    R.length        <- length(restrictions.R) 

    # 0: Residuals
    #---------------------------------------- 
    U   <- matrix(aux$U, ncol=K) 


    # Loop on all correlation parameters
    for(i in 1:R.length){
        # Check for regime independence of the parameter restriction
        if(restrictions.R[i] != 0){

            G.Ri        <- array(NA, aux$G)

            # 1: Calculate support for correlation
            #----------------------------------------
            R.regime    <- aux$R
            R.elements  <- R.regime[upper.tri(R.regime, diag=FALSE)] 
            
            # f(0)
            R.0 <- R.regime
            R.elements.0    <- R.elements
            R.elements.0[i] <- 0
            R.0[lower.tri(R.0, diag=FALSE)] <- R.elements.0
            R.0[upper.tri(R.0, diag=FALSE)] <- R.elements.0
            f.0 <- det(R.0)

            # f(-1)
            R.m1    <- R.regime
            R.elements.m1   <- R.elements
            R.elements.m1[i]<- -1
            R.m1[lower.tri(R.m1, diag=FALSE)] <- R.elements.m1 
            R.m1[upper.tri(R.m1, diag=FALSE)] <- R.elements.m1  
            f.m1    <- det(R.m1)

            # f(1)
            R.1     <- R.regime
            R.elements.1    <- R.elements
            R.elements.1[i] <- 1
            R.1[lower.tri(R.1, diag=FALSE)] <- R.elements.1
            R.1[upper.tri(R.1, diag=FALSE)] <- R.elements.1
            f.1    <- det(R.1)

            # k, l, m
            k   <- (f.1 + f.m1 - 2*f.0) / 2
            l   <- (f.1 - f.m1) / 2
            m   <- f.0

            # Roots of f(r)=kr^2 + lr + m
            a   <- (-l - sqrt(l^2 - 4*k*m)) / (2*k)
            b   <- (-l + sqrt(l^2 - 4*k*m)) / (2*k)

            # 2: Grid for Ri: feasible set
            #---------------------------------------- 
            lower.bound <- min(a,b)
            upper.bound <- max(a,b)
            
            R.grid      <- G.grid$R.grid[i,]

            # 3: Evaluate kernel at the grid points
            #----------------------------------------
            for(j in 1:aux$G){
                # Build covariance matrix at grid point
                Sigma.grid  <- aux$Sigma
                # Likelihood
                R.regime    <- aux$R
                R.elements  <- R.regime[upper.tri(R.regime, diag=FALSE)]

                Ri.grid         <- aux$R
                R.elements.i    <- R.elements
                R.elements.i[i] <- R.grid[j] 

                Ri.grid[lower.tri(Ri.grid, diag=FALSE)] <- R.elements.i
                Ri.grid[upper.tri(Ri.grid, diag=FALSE)] <- R.elements.i

                Sigma.grid  <- diag(aux$S) %*% Ri.grid %*% diag(aux$S)

                # Densities:
                eta.tmp <- dmvnorm(U, sigma=matrix(Sigma.grid,ncol=K)) 
                eta.t   <- matrix(eta.tmp, nrow=1)   

                log.likelihood <- sum( log(eta.t))

                # Prior for Ri, uniform on the support (feasible set)
                log.Ri.prior    <- dunif(R.grid[j], min=lower.bound, max=upper.bound, log=TRUE) 

                # Store in G.Ri
                G.Ri[j] <- log.likelihood + log.Ri.prior
            }
                    
            # Rescale with constant (for exponential to compute) and transform with exponential
            G.Ri    <- G.Ri + abs(max(G.Ri))
            G.Ri    <- exp(G.Ri)

            # 4: Deterministic integration and normalization
            #----------------------------------------
            G.Phi <- array(0, aux$G)
            # Integration
            for(j in 2:G){
                G.Phi[j]    <- G.Phi[j-1] + (R.grid[j]-R.grid[j-1]) * ( abs(G.Ri[j-1]) + abs(G.Ri[j]-G.Ri[j-1])/2)
            }
            # Normalization to 1
            G.Phi   <- G.Phi / G.Phi[aux$G]
            
            # 5: Random draw of Si by numerical interpolation
            #----------------------------------------
            draw    <- runif(n=1, min=0, max=1)

            # Interpolation
            index.min   <- max(which(G.Phi<=draw))
            index.max   <- min(which(G.Phi>=draw))
            Ri.new      <- R.grid[index.min] + ( R.grid[index.max]-R.grid[index.min] ) * (draw-G.Phi[index.min])/(G.Phi[index.max]-G.Phi[index.min])

            # Store as an auxiliary variable
            R.regime    <- aux$R
            R.elements  <- R.regime[upper.tri(R.regime, diag=FALSE)]

            R.elements.i    <- R.elements
            R.elements.i[i] <- Ri.new
            R.regime[upper.tri(R.regime, diag=FALSE)] <- R.elements.i
            R.regime[lower.tri(R.regime, diag=FALSE)] <- R.elements.i 
            aux$R       <- R.regime

        }


    }


    
    # Update auxiliary Sigma matrix
    #----------------------------------------   
    aux$Sigma   <- diag(aux$S) %*% aux$R %*% diag(aux$S)
                                 

    # Output 
    #-----------------------------------------------------------------------------
    return(aux)
}
            
