# File name: 	    M3-A1-3-1.R
# By: 				md
# Purpose:			estimates restricted MSIAH 

# Setting the stage
rm(list=ls())
# setwd("/Users/tomaszwozniak/Dropbox/YesBaby-maybe")
setwd("D:/Dropbox/YesBaby-maybe/")

# Source EM
source("Source-code/MSVAR/Caspur-MSVAR-source-EM.R")          
# Source Gibbs
source("Source-code/Restricted-MSBVAR-Litterman/Caspur-Restricted-MSBVAR-Litterman-source-Gibbs.R")

# Data
load("Data/Extended-dataset.RData")

# set.seed(1234567) # MC1
set.seed(987655) # MC2

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
model.name      <- "M3-A1-3-1"
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
                V.Beta.bar.m1   = solve(.3*diag(K*(1+K*p))),   # Variance of AR coeffs
                S.mean          = 0,    # Mean of lognorm standard deviation coeffs
                S.sd            = 2,    # SD of lognorm standard deviation coeffs
                w               = matrix(9 * diag(M) + matrix(1,M,M),ncol=1)    # Priors: Hidden Markov Chain, transition probabilities       
    )  


# Restrictions
#-----------------------------------------------------------------------------------------------
restrictions    <- G.set.restrictions(model.type="MSIAH", EM.output)

# i) mu_{2,s_2t} = mu_{2}
restrictions$Beta.0[1,2]	<- 1
restrictions$Beta.St[1,2,]	<- 0

# ii) A^{k}_{21,s_2t} = A^{k}_{21}
seq.AR  <- seq(2,7,K)
restrictions$Beta.0[seq.AR,2]	<- 1
restrictions$Beta.St[seq.AR,2,]	<- 0     

# iii) A^{k}_{22,s_2t} = A^{k}_{22}  
seq.AR  <- seq(3,7,K)
restrictions$Beta.0[seq.AR,2]	<- 1
restrictions$Beta.St[seq.AR,2,]	<- 0 

# iv) sigma_{22,s_2t} = sigma_{11}
restrictions$Sigma.0[2,2]	<- 1
restrictions$Sigma.St[2,2,]	<- 0  

# v) sigma_{12,st} = 0 
restrictions$Sigma.0[1,2]	<- 0
restrictions$Sigma.0[2,1]	<- 0
restrictions$Sigma.St[1,2,]	<- 0
restrictions$Sigma.St[2,1,]	<- 0

# vi) A^{k}_{12,s_1t} = 0
seq.AR  <- seq(3,7,K)
restrictions$Beta.0[seq.AR,1]	<- 0
restrictions$Beta.St[seq.AR,1,]	<- 0   
                                   
# Hidden Markov Chain
restrictions$M  <- diag(M^2)
restrictions$dj <- seq(from=M, to=M, length.out=M)

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
N.burnin <- 50000
starting.values <- G.EM.to.Gibbs(model.type="MSIAH", EM.output, priors, restrictions, G) 
Burnin.output <- G.MSBVAR.griddy(N=N.burnin, h=h, priors=priors, restrictions=restrictions, starting.values=starting.values, debug=FALSE, print.iterations=1000)

# Gibbs algorithm: Second run
N.griddy <- 100000
t0 = proc.time()
Griddy.output <- G.MSBVAR.griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Burnin.output$last.draws, debug=FALSE, print.iterations=1000)
t1 <- proc.time()
time.exe <- (t1-t0)/3600
cat("\ntime of execution: ",time.exe[3], " [h]")       

# 2nd try
#-----------------------------------------------------------------------------------------------
# load(paste(path.results, model.name, "-MC2.RData",sep=""))
# qqq   = Griddy.output$last.draws
# qqq$S = apply(qqq$U,2:3,sd)
# priors$w[c(1,5,9),] = 30
# 
# N.griddy <- 10000
# t0 = proc.time()
# Griddy.output <- G.MSBVAR.griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Griddy.output$last.draws, debug=FALSE, print.iterations=1000)
# t1 <- proc.time()
# time.exe <- (t1-t0)/3600
# cat("\ntime of execution: ",time.exe[3], " [h]")       

MHM <- G.MHM.MSVAR(Griddy.output, priors, restrictions, alpha, debug=FALSE, print.iterations=250)

text.string <- paste("\nM: ",M,"\tp: ",p,"\tMHM: ", MHM$MHM, "\talpha: ", MHM$alpha, "\tacceptance.rate: ", MHM$acceptance.rate, " MC3", sep="")
cat(text.string,file=filename.csv,append=TRUE,sep="",fill=FALSE)  

# Save
#-----------------------------------------------------------------------------------------------
# save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, MHM, time.exe, file=paste(path.results, model.name, ".RData",sep=""))  # MC1
save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Griddy.output, MHM, time.exe, file=paste(path.results, model.name, "-MC3.RData",sep=""))  #MC3

# load(paste(path.results, model.name, ".RData",sep=""))
# N.griddy <- 200000
# Griddy.output.f <- G.MSBVAR.griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Griddy.output$last.draws, debug=FALSE, print.iterations=1000)
# MHM <- G.MHM.MSVAR(Griddy.output, priors, restrictions, alpha, debug=FALSE, print.iterations=250)
# save(Y, K, N.griddy, M, p, h, G, priors, restrictions, Griddy.output.f, MHM, file=paste(path.results, model.name, "-f.RData",sep=""))  



load(paste(path.results, model.name, "-MC3.RData",sep="")); ls()
MHM$MHM; MHM$acceptance.rate

plot.ts(t(Griddy.output$posteriors$Beta[1:7,1,]))
plot.ts(t(Griddy.output$posteriors$Beta[8:14,1,]))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$Beta[1,1,]))

plot.ts(t(Griddy.output$posteriors$Beta[1:7,2,]))
plot.ts(t(Griddy.output$posteriors$Beta[8:14,2,]))

plot.ts(t(Griddy.output$posteriors$Beta[1:7,3,]))
plot.ts(t(Griddy.output$posteriors$Beta[8:14,3,]))

plot.ts(t(rbind(Griddy.output$posteriors$Beta[c(1,2,4,6),2,],Griddy.output$posteriors$Beta[c(1,2,4,6),3,])))


plot.ts(t(Griddy.output$posteriors$S[,1,]))
plot.ts(t(Griddy.output$posteriors$S[,2,]))
plot.ts(t(Griddy.output$posteriors$S[,3,]))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$S[1,1,]))

plot.ts(t(Griddy.output$posteriors$R[1,,]))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$R[1,1,]))

plot.ts(t(Griddy.output$posteriors$w))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$w[1,]))
plot.ts(t(Griddy.output$posteriors$PR_TR))

states.y    = matrix(NA,dim(Griddy.output$posteriors$S.t)[1],M)
for (m in 1:M){
   states.y[,m]    = apply(Griddy.output$posteriors$S.t==m,1,sum)/dim(Griddy.output$posteriors$S.t)[2]
}
states.y    = ts(states.y, start=c(1959,5), frequency=12)
plot(states.y)

apply(Griddy.output$posteriors$Beta,1:2,mean)
apply(Griddy.output$posteriors$Sigma,1:2,mean)
apply(Griddy.output$posteriors$S,1:2,mean)
apply(Griddy.output$posteriors$R,1:2,mean)
matrix(apply(Griddy.output$posteriors$PR_TR,1,mean),ncol=M)

mu1   = density(Griddy.output$posteriors$Beta[1,1,])
mu2   = density(Griddy.output$posteriors$Beta[1,2,])
mu3   = density(Griddy.output$posteriors$Beta[1,3,])

plot(mu1, xlim=c(min(mu1$x,mu2$x,mu3$x),max(mu1$x,mu2$x,mu3$x)),ylim=c(0,max(mu1$y,mu2$y,mu3$y)))
lines(mu2)
lines(mu3)

plot(Griddy.output$posteriors$Beta[1,1,],Griddy.output$posteriors$Beta[1,2,])
abline(a=0,b=1,col="red",lwd=2)

plot(Griddy.output$posteriors$Beta[1,1,],Griddy.output$posteriors$Beta[1,3,])
abline(a=0,b=1,col="red",lwd=2)

plot(Griddy.output$posteriors$Beta[1,2,],Griddy.output$posteriors$Beta[1,3,])
abline(a=0,b=1,col="red",lwd=2)


S1   = density(Griddy.output$posteriors$S[1,1,])
S2   = density(Griddy.output$posteriors$S[1,2,])
S3   = density(Griddy.output$posteriors$S[1,3,])

plot(S1, xlim=c(min(S1$x,S2$x,S3$x),max(S1$x,S2$x,S3$x)),ylim=c(0,max(S1$y,S2$y,S3$y)))
lines(S2)
lines(S3)

plot(Griddy.output$posteriors$S[1,1,],Griddy.output$posteriors$S[1,2,])
abline(a=0,b=1,col="red",lwd=2)

plot(Griddy.output$posteriors$S[1,1,],Griddy.output$posteriors$S[1,3,])
abline(a=0,b=1,col="red",lwd=2)

plot(Griddy.output$posteriors$S[1,2,],Griddy.output$posteriors$S[1,3,])
abline(a=0,b=1,col="red",lwd=2)

qqq=Griddy.output
plot.ts(t(rbind(Griddy.output$posteriors$Beta[c(1,2,4,6),2,],Griddy.output$posteriors$Beta[c(1,2,4,6),3,])))
plot.ts(t(rbind(qqq$posteriors$Beta[c(1,2,4,6),2,],qqq$posteriors$Beta[c(1,2,4,6),3,])))

draws    = 49545:50000
plot.ts(t(rbind(qqq$posteriors$Beta[2,2,draws],qqq$posteriors$Beta[2,3,draws])))
plot.ts(t(rbind(Griddy.output$posteriors$Beta[2,2,],Griddy.output$posteriors$Beta[2,3,])))

ss1      = c(7920:14000,14670:14870,15780:18720,19120:20340,21250:22020,25630:27380,27760:29840,30110:31610,34130:36750,37600:38300,40400:41050,41360:43610,44075:46870,49545:50000)

Beta2    = rep(NA,50000)
Beta2[ss1]    = Griddy.output$posteriors$Beta[2,2,ss1]
Beta2[-ss1]    = Griddy.output$posteriors$Beta[2,3,-ss1]


Beta3    = rep(NA,50000)
Beta3[ss1]    = Griddy.output$posteriors$Beta[2,3,ss1]
Beta3[-ss1]    = Griddy.output$posteriors$Beta[2,2,-ss1]
beta22   = cbind(Beta2,Beta3)
plot.ts(beta22)
apply(beta22,2,mean)

qqq$posteriors$Beta[,2,ss1]   = Griddy.output$posteriors$Beta[,3,ss1]
qqq$posteriors$Beta[,3,ss1]   = Griddy.output$posteriors$Beta[,2,ss1]

qqq$posteriors$Sigma[,2,ss1]   = Griddy.output$posteriors$Sigma[,3,ss1]
qqq$posteriors$Sigma[,3,ss1]   = Griddy.output$posteriors$Sigma[,2,ss1]

qqq$posteriors$S[,2,ss1]   = Griddy.output$posteriors$S[,3,ss1]
qqq$posteriors$S[,3,ss1]   = Griddy.output$posteriors$S[,2,ss1]

qqq$posteriors$PR_TR[4,ss1]   = Griddy.output$posteriors$PR_TR[7,ss1]
qqq$posteriors$PR_TR[7,ss1]   = Griddy.output$posteriors$PR_TR[4,ss1]

qqq$posteriors$PR_TR[2,ss1]   = Griddy.output$posteriors$PR_TR[3,ss1]
qqq$posteriors$PR_TR[3,ss1]   = Griddy.output$posteriors$PR_TR[2,ss1]

qqq$posteriors$PR_TR[5,ss1]   = Griddy.output$posteriors$PR_TR[9,ss1]
qqq$posteriors$PR_TR[9,ss1]   = Griddy.output$posteriors$PR_TR[5,ss1]

qqq$posteriors$PR_TR[6,ss1]   = Griddy.output$posteriors$PR_TR[8,ss1]
qqq$posteriors$PR_TR[8,ss1]   = Griddy.output$posteriors$PR_TR[6,ss1]

qqq$posteriors$w[4,ss1]   = Griddy.output$posteriors$w[7,ss1]
qqq$posteriors$w[7,ss1]   = Griddy.output$posteriors$w[4,ss1]
qqq$posteriors$w[2,ss1]   = Griddy.output$posteriors$w[3,ss1]
qqq$posteriors$w[3,ss1]   = Griddy.output$posteriors$w[2,ss1]
qqq$posteriors$w[5,ss1]   = Griddy.output$posteriors$w[9,ss1]
qqq$posteriors$w[9,ss1]   = Griddy.output$posteriors$w[5,ss1]
qqq$posteriors$w[6,ss1]   = Griddy.output$posteriors$w[8,ss1]
qqq$posteriors$w[8,ss1]   = Griddy.output$posteriors$w[6,ss1]

change.state = function(x){
   new   = x
   which2   = x==2
   which3   = x==3
   new[which2] = 3
   new[which3] = 2
   return(new)
}

qqq$posteriors$S.t[,ss1] = apply(Griddy.output$posteriors$S.t[,ss1],2,change.state)
states.y    = matrix(NA,dim(qqq$posteriors$S.t)[1],M)
for (m in 1:M){
   states.y[,m]    = apply(qqq$posteriors$S.t==m,1,sum)/dim(qqq$posteriors$S.t)[2]
}
states.y    = ts(states.y, start=c(1959,5), frequency=12)
plot(states.y)


plot.ts(t(qqq$posteriors$Beta[1:7,1,]))
plot.ts(t(qqq$posteriors$Beta[8:14,1,]))
1-rejectionRate(as.mcmc(qqq$posteriors$Beta[1,1,]))

plot.ts(t(qqq$posteriors$Beta[1:7,2,]))
plot.ts(t(qqq$posteriors$Beta[8:14,2,]))

plot.ts(t(qqq$posteriors$Beta[1:7,3,]))
plot.ts(t(qqq$posteriors$Beta[8:14,3,]))

plot.ts(t(qqq$posteriors$S[,1,]))
plot.ts(t(qqq$posteriors$S[,2,]))
plot.ts(t(qqq$posteriors$S[,3,]))

plot.ts(t(qqq$posteriors$PR_TR))

apply(qqq$posteriors$Beta,1:2,mean)
apply(qqq$posteriors$Sigma,1:2,mean)
apply(qqq$posteriors$S,1:2,mean)
apply(qqq$posteriors$R,1:2,mean)
matrix(apply(qqq$posteriors$PR_TR,1,mean),ncol=M)

# Save
#-----------------------------------------------------------------------------------------------
save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output,qqq, MHM, time.exe, file=paste(path.results, model.name, ".RData",sep=""))  

MHM$MHM
MHM <- G.MHM.MSVAR(qqq, priors, restrictions, alpha, debug=FALSE, print.iterations=250)
text.string <- paste("\nM: ",M,"\tp: ",p,"\tMHM: ", MHM$MHM, "\talpha: ", MHM$alpha, "\tacceptance.rate: ", MHM$acceptance.rate, sep="")
MHM$MHM
