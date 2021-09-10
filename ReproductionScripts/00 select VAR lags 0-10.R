# File name: 	    VAR-0.R
# By: 				md
# Purpose:			estimates unrestricted VAR


# Source EM
source("D:/YesBaby/Source-code/MSVAR/Caspur-MSVAR-source-EM.R")          
# Source Gibbs
source("D:/YesBaby/Source-code/Restricted-BVAR-Litterman/Caspur-Restricted-BVAR-Litterman-source-Gibbs.R")
  
# Data
load(file="D:/YesBaby/Data/Extended-dataset.RData")


# Set seed
# set.seed(62380)   # The first run
set.seed(123456)

# Model name
#---------------------------------------------------------------------------------------------------
model   <- "VAR-s03"


# Parameters
#---------------------------------------------------------------------------------------------------
Y   <- data
K   <- dim(Y)[2]

p.min   <- 0
p.max   <- 10

h           <- 1     # Forecast horizon
G           <- 50    # Grid size for griddy Gibbs

alpha <- 0.05        # Truncation for MHM
path.results    <- "D:/YesBaby/Results2013/ModelSelectionFull/"
# path.load    <- "./Results/ModelSelection/"
filename.csv    <- paste(path.results, model, ".csv",sep="")   

# For the prior: Standard deviations of the AR residuals of individual series, 17 lags
standard.deviations <- array(NA,K)
for(series in 1:K){
  AR                          <- ar(data[,series],aic=FALSE,order.max=17)
  AR.resid                    <- AR$resid
  standard.deviations[series] <- sd(AR.resid[!is.na(AR.resid)])
}

# Loop for Gibbs + MHM
#---------------------------------------------------------------------------------------------------
for(p in p.min:p.max){
    model.name  <- paste(model,"-p",p, sep="")  

    # OLS as starting values for Gibbs
    #-----------------------------------------------------------------------------------------------
    starting.values <- G.Initialization.BVAR(Y, p, G)
 
    # Priors - original: 
    #-----------------------------------------------------------------------------------------------
    priors  <- list(Beta.bar        = array(0, c(K*(1+K*p),1)), # Mean of AR coeffs
                    V.Beta.bar.m1   = solve(0.3*diag(K*(1+K*p))),   # Variance of AR coeffs
#                     V.Beta.bar.m1   = solve(10*diag(K*(1+K*p))),   # Variance of AR coeffs
                    S.mean          = 0,    # Mean of lognorm standard deviation coeffs
                    S.sd            = 2    # SD of lognorm standard deviation coeffs 
    ) 
    
    # Priors - shrinkage
    #----------------------------------------------------------------------------------------------      
#     Beta.lambda.1   <- 0.3                  # Litterman priors         
#     Beta.lambda.2   <- 1                    # Litterman priors    
#     Beta.lambda.4   <- 2                    # Litterman priors
#     Beta.c          <- -0.13412             # Litterman priors
#     Beta.sd         <- standard.deviations  # Litterman priors 
#     
#     # Form the variance of the priors
#     # Intercepts
#     priors.beta <- (Beta.sd^2) * (Beta.lambda.4^2)
#     # Matrix common to all lags
#     litterman.sd    <- array(1, c(K,K))
#     for(j in 1:K){
#       litterman.sd[j,]    <- Beta.sd^2 / Beta.sd[j]^2
#     }
#     # Loop for lags
#     if (p!=0){
#       for(lag in 1:p){
#         priors.lag  <- litterman.sd * Beta.lambda.1^2 * Beta.lambda.2^2 * exp( Beta.c * lag - Beta.c)
#         diag(priors.lag)    <- diag(priors.lag) / Beta.lambda.2^2
#         priors.beta <- rbind(priors.beta, priors.lag)
#       }
#     }
#     
#     # Invert and form matrix
#     V.Beta.bar.m1   <-diag(as.vector(1/matrix(priors.beta, ncol=1)))
#     
#     priors  <- list(Beta.bar        = array(0, c(K*(1+K*p),1)), # Mean of AR coeffs
#                     V.Beta.bar.m1   = V.Beta.bar.m1, # Variance of AR coeffs
#                     S.mean          = 0,    # Mean of lognorm standard deviation coeffs
#                     S.sd            = 2     # SD of lognorm standard deviation coeffs 
#     ) 
#     
    # Restrictions
    #-----------------------------------------------------------------------------------------------
    restrictions    <- list(Beta = array(1, c(1+K*p, K)),
                            Sigma= array(1, c(K,K))
                            ) 
    # Gibbs: burn in
    #-----------------------------------------------------------------------------------------------
    N.burnin        <- 10000
    Burnin.output   <- G.BVAR.Griddy(N=N.burnin, h=h, priors=priors, restrictions=restrictions, starting.values=starting.values, debug=FALSE, print.iterations=100)
                                     


    # Gibbs
    #-----------------------------------------------------------------------------------------------
    N.griddy        <- 50000
    starting.values <- Burnin.output$last.draws
    Griddy.output   <- G.BVAR.Griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=starting.values, debug=FALSE, print.iterations=100)

    # Save
    save(Y, K, N.burnin, N.griddy, p, h, G, priors, restrictions, Burnin.output, Griddy.output, file = paste(path.results, model.name, ".RData",sep=""))


    # MHM
    #-----------------------------------------------------------------------------------------------
    MHM         <- G.MHM.BVAR.Geweke(Gibbs.output=Griddy.output, priors=priors, restrictions=restrictions, alpha=alpha, debug=FALSE, print.iterations=250)
    text.string <- paste("\np: ",p,"\tMHM: ", MHM$MHM, "\tacceptance.rate: ", MHM$acceptance.rate, sep="")
    cat(text.string,file=filename.csv,append=TRUE,sep="",fill=FALSE)

    # Save
    
    # First run save
    #     save(Y, K, N.burnin, N.griddy, p, h, G, priors, restrictions, Burnin.output, Griddy.output, MHM, file = paste(path.results, model.name, ".RData",sep=""))
    # Second run save
    # save(Y, K, N.burnin, N.griddy, p, h, G, priors, restrictions, Burnin.output, Griddy.output, MHM, file = paste(path.results, model.name, "a.RData",sep=""))
    # 100 prior run
    # save(Y, K, N.burnin, N.griddy, p, h, G, priors, restrictions, Burnin.output, Griddy.output, MHM, file = paste(path.results, model.name, "b.RData",sep=""))
    # shrinkage prior
    # save(Y, K, N.burnin, N.griddy, p, h, G, priors, restrictions, Burnin.output, Griddy.output, MHM, file = paste(path.results, model.name, "s.RData",sep=""))
    # overall-shrinkage prior
    save(Y, K, N.burnin, N.griddy, p, h, G, priors, restrictions, Burnin.output, Griddy.output, MHM, file = paste(path.results, model.name, ".RData",sep=""))
    
}
