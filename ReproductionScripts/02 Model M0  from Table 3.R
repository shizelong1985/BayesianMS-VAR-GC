# File name: 	    MSVAR-0-M2.R
# By: 				md
# Purpose:			estimates MSIAH 

rm(list=ls())
setwd("D:/Dropbox/YesBaby-maybe/")
setwd("/Users/szeridan/Dropbox/YesBaby-maybe")

# Source EM
source("Source-code/MSVAR/Caspur-MSVAR-source-EM.R")          
# Source Gibbs
source("Source-code/Restricted-MSBVAR-Litterman/Caspur-Restricted-MSBVAR-Litterman-source-Gibbs.R")
 
# Data
load(file="Data/Extended-dataset.RData")

# set.seed(62380) #MC1
set.seed(123456) #MC2

# Model name
#---------------------------------------------------------------------------------------------------
 

# Parameters
#---------------------------------------------------------------------------------------------------
Y   <- data
K   <- dim(Y)[2]

M  <- 3
p  = 3
h   <- 1    # Forecast horizon
G   <- 50   # Grid size for griddy Gibbs

alpha <- 0.05 # Truncation for MHM
path.results    <- "./Results2014/Hypotheses/"
model   <- "MSVAR-M3-"
model.name      <- paste(model,"p", p,sep="")
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

    # 2/3 Sigma
    restrictions$Sigma.0 <- array(0,c(K,K))
    restrictions$Sigma.St<- array(1,c(K,K,M))
  
    # 3/3 Hidden Markov Chain
    restrictions$M  <- diag(M^2)
    restrictions$dj <- seq(from=M, to=M, length.out=M)   

    # Gibbs:
    #-----------------------------------------------------------------------------------------------

    # Gibbs algorithm: Burnin
    N.burnin        <- 10000
    starting.values <- G.EM.to.Gibbs(model.type="MSIAH", EM.output, priors, restrictions, G) 
    Burnin.output   <- G.MSBVAR.griddy(N=N.burnin, h=h, priors=priors, restrictions=restrictions, starting.values=starting.values, debug=FALSE, print.iterations=1000)
    
    # Gibbs algorithm: Second run
    N.griddy        <- 100000     
    t0 = proc.time()
    Griddy.output   <- G.MSBVAR.griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Burnin.output$last.draws, debug=FALSE, print.iterations=1000)
    t1         = proc.time()
    time.exe   = (t1-t0)/3600
    cat("\ntime of execution: ",time.exe[3], " [h]")

    save(Y, K, N.burnin, N.griddy, M, p, h, G, EM.output, Burnin.output, Griddy.output, priors, restrictions, file=paste(path.results, model.name, ".RData",sep=""))

    # MHM
    MHM <- G.MHM.MSVAR(Griddy.output, priors, restrictions, alpha, debug=FALSE, print.iterations=250)

    text.string <- paste("\nM: ",M,"\tp: ",p,"\tMHM: ", MHM$MHM, "\talpha: ", MHM$alpha, "\tacceptance.rate: ", MHM$acceptance.rate, " MC2", sep="")
    cat(text.string,file=filename.csv,append=TRUE,sep="",fill=FALSE)  
   
#     save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, MHM, time.exe, file=paste(path.results, model.name, ".RData",sep="")) #MC1
save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, MHM, time.exe, file=paste(path.results, model.name, "-MC2.RData",sep="")) #MC2

# load(paste(path.results, model.name, ".RData",sep=""))
# N.griddy <- 200000
# Griddy.output.f <- G.MSBVAR.griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Griddy.output$last.draws, debug=FALSE, print.iterations=1000)
# MHM <- G.MHM.MSVAR(Griddy.output, priors, restrictions, alpha, debug=FALSE, print.iterations=250)
# save(Y, K, N.griddy, M, p, h, G, priors, restrictions, Griddy.output.f, MHM, file=paste(path.results, model.name, "-f.RData",sep=""))  

model   <- "MSVAR-M3-"
model.name      <- paste(model,"p", p,sep="")

load(paste(path.results, model.name, "-MC2.RData",sep=""))
load(paste(path.results,"MDD-MSVAR-M3-s03p3.RData",sep=""))
MHM$MHM; MHM$acceptance.rate



####################################################################################
# Prepare for the results to the paper
# First: label switching
####################################################################################
# WARNING! By mistake you did label switching for model "M3-M1-3"
####################################################################################
# NOTHING IS NEEDED UP TO:
# plot.ts(t(Griddy.output$posteriors$S[1,,]))
# plot.ts(Griddy.output$posteriors$S[1,1,])
# abline(h=15)
# plot.ts(t(Griddy.output$posteriors$Beta[6,,]))
# plot.ts(t(Griddy.output$posteriors$PR_TR[c(1,4,7),]))
# 
# s13 = (Griddy.output$posteriors$S[1,1,])<15
# qqq = Griddy.output$posteriors
# 
# 
# qqq$S[1,1,s13] = Griddy.output$posteriors$S[1,3,s13]
# qqq$S[1,3,s13] = Griddy.output$posteriors$S[1,1,s13]
# 
# qqq$Beta[1:7,1,s13] = Griddy.output$posteriors$Beta[1:7,3,s13]
# qqq$Beta[1:7,3,s13] = Griddy.output$posteriors$Beta[1:7,1,s13]
# 
# qqq$Sigma[1:3,1,s13] = Griddy.output$posteriors$Sigma[1:3,3,s13]
# qqq$Sigma[1:3,3,s13] = Griddy.output$posteriors$Sigma[1:3,1,s13]
# 
# qqq$R[,1,s13] = Griddy.output$posteriors$R[,3,s13]
# qqq$R[,3,s13] = Griddy.output$posteriors$R[,1,s13]
# 
# qqq$PR_TR[1,s13] = Griddy.output$posteriors$PR_TR[9,s13]
# qqq$PR_TR[9,s13] = Griddy.output$posteriors$PR_TR[1,s13]
# qqq$PR_TR[3,s13] = Griddy.output$posteriors$PR_TR[7,s13]
# qqq$PR_TR[7,s13] = Griddy.output$posteriors$PR_TR[3,s13]
# qqq$w[1,s13] = Griddy.output$posteriors$w[9,s13]
# qqq$w[9,s13] = Griddy.output$posteriors$w[1,s13]
# qqq$w[3,s13] = Griddy.output$posteriors$w[7,s13]
# qqq$w[7,s13] = Griddy.output$posteriors$w[3,s13]
# 
# qqq$S.t[,s13][which(Griddy.output$posteriors$S.t[,s13]==1)] = 3
# qqq$S.t[,s13][which(Griddy.output$posteriors$S.t[,s13]==3)] = 1
# plot.ts(t(qqq$Beta[1:7,1,]))
# plot.ts(t(qqq$Beta[8:14,1,]))
# 
# plot.ts(t(qqq$Beta[1:7,2,]))
# plot.ts(t(qqq$Beta[8:14,2,]))
# 
# plot.ts(t(qqq$Beta[1:7,3,]))
# plot.ts(t(qqq$Beta[8:14,3,]))
# 
# plot.ts(t(qqq$S[,1,]))
# plot.ts(t(qqq$S[,2,]))
# plot.ts(t(qqq$S[,3,]))
# 
# plot.ts(t(qqq$Sigma[,1,]))
# plot.ts(t(qqq$Sigma[,2,]))
# plot.ts(t(qqq$Sigma[,3,]))
# 
# plot.ts(t(qqq$R[1,,]))
# plot.ts(t(qqq$PR_TR))
# 
# states.y    = matrix(NA,dim(qqq$S.t)[1],M)
# for (m in 1:M){
#    states.y[,m]    = apply(qqq$S.t==m,1,sum)/dim(qqq$S.t)[2]
# }
# colnames(states.y) = c("State 1","State 2","State 3")
# states.y    = ts(states.y, start=c(1959,5), frequency=12)
# plot(states.y)
# 
# apply(qqq$Beta,1:2,mean)
# apply(qqq$Sigma,1:2,mean)
# apply(qqq$S,1:2,mean)
# apply(qqq$R,1:2,mean)
# matrix(apply(qqq$PR_TR,1,mean),ncol=M)
# HERE
####################################################################################












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

plot.ts(t(Griddy.output$posteriors$Sigma[,1,]))
plot.ts(t(Griddy.output$posteriors$Sigma[,2,]))
plot.ts(t(Griddy.output$posteriors$Sigma[,3,]))

plot.ts(t(Griddy.output$posteriors$R[1,,]))

plot.ts(t(Griddy.output$posteriors$R[1,,]))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$R[1,1,]))

plot.ts(t(Griddy.output$posteriors$PR_TR))
1-rejectionRate(as.mcmc(Griddy.output$posteriors$PR_TR[1,]))

states.y    = matrix(NA,dim(Griddy.output$posteriors$S.t)[1],M)
for (m in 1:M){
   states.y[,m]    = apply(Griddy.output$posteriors$S.t==m,1,sum)/dim(Griddy.output$posteriors$S.t)[2]
}
colnames(states.y) = c("State 1","State 2","State 3")
states.y    = ts(states.y, start=c(1959,5), frequency=12)
plot(states.y)

pdf(file="/Users/tomaszwozniak/Dropbox/Yes Baby/Matchik/Latex/13- paper/Plots/Regimes.pdf",width=10,height=7)
plot(states.y, main="", xlab="")
dev.off()

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




####################################################
# Report the parameters to the papaer
####################################################
sink("/Users/tomaszwozniak/Dropbox/YesBaby-maybe/Results2014/RSLTS-MS3VAR3-post.txt")
round(t(matrix(apply(Griddy.output$posteriors$Beta,1:2,mean)[,1],ncol=2)),3)
round(t(matrix(apply(Griddy.output$posteriors$Beta,1:2,sd)[,1],ncol=2)),3)

round(t(matrix(apply(Griddy.output$posteriors$Beta,1:2,mean)[,2],ncol=2)),3)
round(t(matrix(apply(Griddy.output$posteriors$Beta,1:2,sd)[,2],ncol=2)),3)

round(t(matrix(apply(Griddy.output$posteriors$Beta,1:2,mean)[,3],ncol=2)),3)
round(t(matrix(apply(Griddy.output$posteriors$Beta,1:2,sd)[,3],ncol=2)),3)

round(apply(Griddy.output$posteriors$S,1:2,mean),3)
round(apply(Griddy.output$posteriors$S,1:2,sd),3)

round(apply(Griddy.output$posteriors$R,1:2,mean),3)
round(apply(Griddy.output$posteriors$R,1:2,sd),3)

round(matrix(apply(Griddy.output$posteriors$PR_TR,1,mean),ncol=M),3)
round(matrix(apply(Griddy.output$posteriors$PR_TR,1,sd),ncol=M),3)
sink()









