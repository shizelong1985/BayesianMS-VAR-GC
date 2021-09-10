# File name: 	    M3-A2-3-1.R
# By: 				md
# Purpose:			estimates restricted MSIAH 

# Setting the stage
rm(list=ls())
setwd("D:/Dropbox/YesBaby-maybe/")
# setwd("/Users/tomaszwozniak/Dropbox/YesBaby-maybe")

# Source EM
source("Source-code/MSVAR/Caspur-MSVAR-source-EM.R")          
# Source Gibbs
source("Source-code/Restricted-MSBVAR-Litterman/Caspur-Restricted-MSBVAR-Litterman-source-Gibbs.R")

# Data
load("Data/Extended-dataset.RData")

# set.seed(12345678) #MC1
set.seed(87654) #MC2

# Parameters
#---------------------------------------------------------------------------------------------------
Y   <- data
K   <- dim(Y)[2]

M   <- 3
p   <- 3

h   <- 1    # Forecast horizon
G   <- 50   # Grid size for griddy Gibbs  

alpha <- 0.05 # Truncation for MHM
path.results    <- "Results2014/Hypotheses/"
model.name      <- "M3-A2-3-1"
filename.csv    <- paste(path.results, model.name, "-alpha", alpha,".csv",sep="") 

# EM
#-----------------------------------------------------------------------------------------------
EM.output <- EM.MSIAH(Y=Y, M=M, p=p,  max.iter=1000, convergence=0.0000001, diagonal.PR_TR=0.9, plot=FALSE)


# Priors: 
#-----------------------------------------------------------------------------------------------
priors  <- list(Beta.bar        = array(0, c(K*(1+K*p),1)), # Mean of AR coeffs
                V.Beta.bar.m1   = solve(.3*diag(K*(1+K*p))),   # Variance of AR coeffs
                S.mean          = 0,    # Mean of lognorm standard deviation coeffs
                S.sd            = 2,    # SD of lognorm standard deviation coeffs
                w               = matrix(9 * diag(M) + matrix(1,M,M),ncol=1)    # Priors: Hidden Markov Chain, transition probabilities       
    )  


# Restrictions
#-----------------------------------------------------------------------------------------------
restrictions    <- G.set.restrictions(model.type="MSIAH", EM.output)

# i) p_11 = p21
restrictions$M  <- rbind(diag(M),diag(M),diag(M))
restrictions$dj <- 3


# showtime!
cat("Beta 0")
print(t(restrictions$Beta.0))
cat("Beta St")
for(regime in 1:M) print(t(restrictions$Beta.St[,,regime]))
cat("Sigma 0")
print(restrictions$Sigma.0)
cat("Sigma St")
print(restrictions$Sigma.St) 


# Gibbs
#-----------------------------------------------------------------------------------------------
# Gibbs algorithm: Burnin
N.burnin <- 10000
starting.values <- G.EM.to.Gibbs(model.type="MSIAH", EM.output, priors, restrictions, G) 
Burnin.output <- G.MSBVAR.griddy(N=N.burnin, h=h, priors=priors, restrictions=restrictions, starting.values=starting.values, debug=FALSE, print.iterations=1000)

# Gibbs algorithm: Second run
N.griddy <- 100000
t0 = proc.time()
Griddy.output <- G.MSBVAR.griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Burnin.output$last.draws, debug=FALSE, print.iterations=1000)
t1 <- proc.time()
time.exe <- (t1-t0)/3600
cat("\ntime of execution: ",time.exe[3], " [h]")       

# MHM
#-----------------------------------------------------------------------------------------------
MHM <- G.MHM.MSVAR(Griddy.output, priors, restrictions, alpha, debug=FALSE, print.iterations=250)

text.string <- paste("\nM: ",M,"\tp: ",p,"\tMHM: ", MHM$MHM, "\talpha: ", MHM$alpha, "\tacceptance.rate: ", MHM$acceptance.rate, " MC2", sep="")
cat(text.string,file=filename.csv,append=TRUE,sep="",fill=FALSE)  

# Save
#-----------------------------------------------------------------------------------------------
# save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, MHM, time.exe, file=paste(path.results, model.name, ".RData",sep=""))  #MC1
save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, MHM, time.exe, file=paste(path.results, model.name, "-MC2.RData",sep=""))  #MC2

# load(paste(path.results, model.name, ".RData",sep=""))
# N.griddy <- 200000
# Griddy.output.f <- G.MSBVAR.griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Griddy.output$last.draws, debug=FALSE, print.iterations=1000)
# MHM <- G.MHM.MSVAR(Griddy.output, priors, restrictions, alpha, debug=FALSE, print.iterations=250)
# save(Y, K, N.griddy, M, p, h, G, priors, restrictions, Griddy.output.f, MHM, file=paste(path.results, model.name, "-f.RData",sep=""))  

load(paste(path.results, model.name, "-MC2.RData",sep="")); ls()
MHM$MHM; MHM$acceptance.rate

plot.ts(t(Griddy.output$posteriors$Beta[1:7,1,]))
plot.ts(t(Griddy.output$posteriors$Beta[8:14,1,]))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$Beta[1,1,]))

plot.ts(t(Griddy.output$posteriors$Beta[1:7,2,]))
plot.ts(t(Griddy.output$posteriors$Beta[8:14,2,]))

plot.ts(t(Griddy.output$posteriors$Beta[1:7,3,]))
plot.ts(t(Griddy.output$posteriors$Beta[8:14,3,]))

plot.ts(t(Griddy.output$posteriors$S[,1,]))
plot.ts(t(Griddy.output$posteriors$S[,2,]))
plot.ts(t(Griddy.output$posteriors$S[,3,]))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$S[1,1,]))

plot.ts(t(Griddy.output$posteriors$R[1,,]))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$R[1,1,]))

plot.ts(t(Griddy.output$posteriors$w))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$w[1,]))

states.y    = matrix(NA,dim(Griddy.output$posteriors$S.t)[1],M)
for (m in 1:M){
   states.y[,m]    = apply(Griddy.output$posteriors$S.t==m,1,sum)/dim(Griddy.output$posteriors$S.t)[2]
}
states.y    = ts(states.y, start=c(1959,5), frequency=12)
plot(states.y)
