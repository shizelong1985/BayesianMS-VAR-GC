# File name: 	    M3-A3-rank2.R
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

# set.seed(23456) #MC1
# set.seed(987654) #MC2
set.seed(456765) #MC3


# Parameters
#---------------------------------------------------------------------------------------------------
Y   <- data
K   <- dim(Y)[2]

M   <- 3
p   <- 3

h   <- 1    # Forecast horizon
G   <- 100   # Grid size for griddy Gibbs  

alpha <- 0.05 # Truncation for MHM
# path.results    <- "/Users/tomaszwozniak/Dropbox/YesBaby-maybe/Results2013/Hypotheses-full/"
path.results    <- "./Results2014/Hypotheses/"
model.name      <- "M3-A3-rank2"
filename.csv    <- paste(path.results, model.name, "-alpha", alpha,".csv",sep="") 

# EM
#-----------------------------------------------------------------------------------------------
EM.output <- EM.MSIAH(Y=Y, M=M, p=p,  max.iter=1000, convergence=0.0000001, diagonal.PR_TR=0.9, plot=FALSE)


# Priors: 
#-----------------------------------------------------------------------------------------------
# For the prior: Standard deviations of the AR residuals of individual series, 17 lags
standard.deviations <- array(NA,K)
for(series in 1:K){
  AR                          <- ar(data[,series],aic=FALSE,order.max=17)
  AR.resid                    <- AR$resid
  standard.deviations[series] <- sd(AR.resid[!is.na(AR.resid)])
} 

priors  <- list(Beta.bar        = array(0, c(K*(1+K*p),1)), # Mean of AR coeffs
                V.Beta.bar.m1   = solve(.3*diag(K*(1+K*p))), # Variance of AR coeffs
                S.mean          = 0, # Mean of lognorm standard deviation coeffs
                S.sd            = 2, # SD of lognorm standard deviation coeffs
                w               = matrix(9 * diag(M) + matrix(1,M,M),ncol=1) # Hidden Markov Chain, transition probabilities
    		)


# Restrictions
#-----------------------------------------------------------------------------------------------
restrictions    <- G.set.restrictions(model.type="MSIAH", EM.output)

# ii) restrict the autoregressive parameters in the first regime and in the first equation to 0
restrictions$Beta.St[,1,1] = 0


# showtime!
cat("Beta 0")
print(t(restrictions$Beta.0))
cat("Beta St")
for(regime in 1:M) print(t(restrictions$Beta.St[,,regime]))
cat("Sigma 0")
print(restrictions$Sigma.0)
cat("Sigma St")
print(restrictions$Sigma.St) 

# source("Source-code/Restricted-MSBVAR-Litterman/Caspur-Restricted-MSBVAR-Litterman-source-Gibbs.R")

# Gibbs
#-----------------------------------------------------------------------------------------------
# # Gibbs algorithm: Burnin
N.burnin <- 1000
starting.values <- G.EM.to.Gibbs(model.type="MSIAH", EM.output, priors, restrictions, G)
starting.values$PR_TR <- matrix(1/3,M,M)
Burnin.output <- G.MSBVAR.griddy.Rank2(N=N.burnin, h=h, priors=priors, restrictions=restrictions, starting.values=starting.values, debug=FALSE, print.iterations=1000)


# Gibbs algorithm: Second run
N.griddy <- 100000
t0 = proc.time()
Griddy.output <- G.MSBVAR.griddy.Rank2(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Burnin.output$last.draws, debug=FALSE, print.iterations=1000)
t1 <- proc.time()
time.exe <- (t1-t0)/3600
cat("\ntime of execution: ",time.exe[3], " [h]")       

# MHM
#-----------------------------------------------------------------------------------------------
# MHM <- G.MHM.MSVAR.reducedRankP(Griddy.output, priors, restrictions, rankP=2, alpha, debug=FALSE, print.iterations=250)
# 
# 
# text.string <- paste("\nM: ",M,"\tp: ",p,"\tMHM: ", MHM$MHM, "\talpha: ", MHM$alpha, "\tacceptance.rate: ", MHM$acceptance.rate, "MC2", sep="")
# cat(text.string,file=filename.csv,append=TRUE,sep="",fill=FALSE)  

# save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, time.exe, file=paste(path.results, model.name, ".RData",sep=""))  #MC1
# save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, MHM, time.exe, file=paste(path.results, model.name, "-MC3.RData",sep=""))  #MC2
save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, time.exe, file=paste(path.results, model.name, "-MC4.RData",sep=""))  #MC2


# Transform the draws to obtain the draws from the restricted model (sic!)
load(paste(path.results, model.name, "-MC4.RData",sep=""))
qqq      = Griddy.output$posteriors
seq.AR   = seq(3,(1+K*p),K)
Beta     = array(NA,c(1+K*p,K,M))
for (s in 1:dim(qqq$Beta)[3]){
   
   for (m in 1:M){
      Beta[,,m]     = qqq$Beta[,m,s]
   }
   PR_TR    = matrix(qqq$PR_TR[,s],ncol=3)
   
   
   # Restriction A3 rank2 (iii): 
   for(i in seq.AR){
      Beta[i,1,1]   <- Beta[i,1,3]*((PR_TR[1,3]*PR_TR[2,2]-PR_TR[2,3]*PR_TR[1,2])/(PR_TR[1,1]*PR_TR[2,2]-PR_TR[2,1]*PR_TR[1,2]))
   }
   
   # Restriction A3 rank2 (i) and (ii): 
   seq.ARR  = (1:(1+K*p))[-seq.AR]
   for(i in seq.ARR){
      Beta[i,1,1]   <- - Beta[i,1,2]*((PR_TR[1,2]-PR_TR[2,2])/(PR_TR[1,1]-PR_TR[1,2])) - Beta[i,1,3]*((PR_TR[1,3]-PR_TR[2,3])/(aux$PR_TR[1,1]-PR_TR[2,1]))
   }
   
   for (m in 1:M){
      Griddy.output$posteriors$Beta[,m,s] = as.vector(Beta[,,m])
   }   
}


# Here select the observations...
stationary  = vector("logical",dim(qqq$Beta)[3])
for (s in 1:dim(qqq$Beta)[3]){
   sum.a    = matrix(0,K,K)
   for (pin in 1:p){
      sum.a    = sum.a + matrix(Griddy.output$posteriors$Beta[,1,s], ncol=K)[(pin-1)*K+1:K+1,] 
   }
   if (max(abs(eigen(sum.a)$values))<1){
      stationary[s] = TRUE
   }
}
length(stationary[stationary==TRUE])/length(stationary)
length(stationary[stationary==TRUE])

Griddy.output$posteriors$Beta    = Griddy.output$posteriors$Beta[,,stationary]
Griddy.output$posteriors$Sigma   = Griddy.output$posteriors$Sigma[,,stationary]
Griddy.output$posteriors$S       = Griddy.output$posteriors$S[,,stationary]
Griddy.output$posteriors$PR_TR   = Griddy.output$posteriors$PR_TR[,stationary]
Griddy.output$posteriors$w       = Griddy.output$posteriors$w[,stationary]
Griddy.output$posteriors$S.t     = Griddy.output$posteriors$S.t[,stationary]
Griddy.output$posteriors$c       = Griddy.output$posteriors$c[stationary]
Griddy.output$posteriors$R       = array(Griddy.output$posteriors$R[,,stationary],c(1,3,75564))

MHM <- G.MHM.MSVAR.reducedRankP(Griddy.output, priors, restrictions, rankP=2, alpha=0.05, debug=FALSE, print.iterations=250)
MHM$MHM; MHM$acceptance.rate

save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, MHM, time.exe, file=paste(path.results, model.name, "-MC5.RData",sep=""))


# compute c from PR_TR
# ccc   = function(x){return((x[3]-x[2])/(x[1]-x[2]))}
# cc    = apply(Burnin.output$posteriors$PR_TR,2,ccc)

# path.results    <- "./Results2014/Hypotheses/"
# model.name      <- "M3-A3-rank2"


load(paste(path.results, model.name, "-MC4.RData",sep=""))
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

plot.ts(t(Griddy.output$posteriors$PR_TR))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$PR_TR[1,]))

states.y    = matrix(NA,dim(Griddy.output$posteriors$S.t)[1],M)
for (m in 1:M){
   states.y[,m]    = apply(Griddy.output$posteriors$S.t==m,1,sum)/dim(Griddy.output$posteriors$S.t)[2]
}
states.y    = ts(states.y, start=c(1959,5), frequency=12)
plot(states.y)
