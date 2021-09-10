# File name: 	    VARp14-restricted.R
# By: 				md
# Purpose:			estimates restricted VAR 


# Setting the stage
rm(list=ls())
setwd("D:/Dropbox/YesBaby-maybe")

# Source Gibbs
source("Source-code/Restricted-BVAR-Litterman/Caspur-Restricted-BVAR-Litterman-source-Gibbs.R")

# Data
load("Data/Extended-dataset.RData") 


# Parameters
#---------------------------------------------------------------------------------------------------
Y   <- data
K   <- dim(Y)[2]

p   <- 14


h   <- 1    # Forecast horizon
G   <- 50   # Grid size for griddy Gibbs  

alpha <- 0.05 # Truncation for MHM
path.results    <- "Results/M3-restrictions/"
model.name      <- "VARp14-restricted"  
filename.csv    <- paste(path.results, model.name, "-alpha", alpha,".csv",sep="") 

t0 = proc.time()     

# OLS as starting values for Gibbs
#-----------------------------------------------------------------------------------------------
starting.values   <- G.Initialization.BVAR(Y, p, G)


# Priors: 
#-----------------------------------------------------------------------------------------------

priors  <- list(Beta.bar        = array(0, c(K*(1+K*p),1)), # Mean of AR coeffs
                V.Beta.bar.m1   = solve(.3*diag(K*(1+K*p))),   # Variance of AR coeffs
                S.mean          = 0,    # Mean of lognorm standard deviation coeffs
                S.sd            = 2    # SD of lognorm standard deviation coeffs 
                ) 

# Restrictions
#-----------------------------------------------------------------------------------------------
restrictions    <- list(Beta = array(1, c(1+K*p, K)),
                        Sigma= array(1, c(K,K))
                        ) 

# Granger non causality from money to income
seq.AR  <- seq(3,29,2)
restrictions$Beta[seq.AR,1] <- 0

# showtime!
cat("Beta")
print(t(restrictions$Beta))    
cat("Sigma")
print(restrictions$Sigma) 

# Gibbs:
#-----------------------------------------------------------------------------------------------

# Gibbs algorithm: Burnin
N.burnin <- 10000
Burnin.output <- G.BVAR.Griddy(N=N.burnin, h=h, priors=priors, restrictions=restrictions, starting.values=starting.values, debug=FALSE, print.iterations=100)

# Gibbs algorithm: Second run
N.griddy        <- 50000
Griddy.output   <- G.BVAR.Griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Burnin.output$last.draws, debug=FALSE, print.iterations=100)

# MHM
#-----------------------------------------------------------------------------------------------
MHM <- G.MHM.BVAR(Gibbs.output=Griddy.output, priors=priors, restrictions=restrictions, alpha=alpha, debug=TRUE, print.iterations=250)
        
text.string <- paste("\np: ",p,"\tMHM: ", MHM$MHM, "\talpha: ", MHM$alpha, "\tacceptance.rate:", MHM$acceptance.rate, sep="")
cat(text.string,file=filename.csv,append=TRUE,sep="",fill=FALSE)  

t1 <- proc.time()
time.exe <- (t1-t0)/3600
cat("\ntime of execution: ",time.exe[3], " [h]")       
                                                   

# Save Gibbs outputs
save(Y, K, N.burnin, N.griddy, p, h, G, priors, restrictions, Burnin.output, Griddy.output, MHM, error.code, file = paste(path.results, model.name, ".RData",sep=""))
