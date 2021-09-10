
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@pm.me

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (2017) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



rm(list=ls())
# Source EM
source("BayesianMSVAR/MSVAR/Caspur-MSVAR-source-EM.R")          
# Source Gibbs
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Caspur-Restricted-MSBVAR-Litterman-source-Gibbs.R")
 
# Data
load(file="data.RData")

set.seed(123456)

# Parameters of the simulation
#---------------------------------------------------------------------------------------------------
Y           = data
K           = dim(Y)[2]

M           = 1
p           = 1
h           = 1    # Forecast horizon
G           = 50   # Grid size for griddy Gibbs

alpha       = 0.05 # Truncation for MHM
path.results= "./Results/"
model       = "MSVAR-M2-"
model.name  = paste(model,"p", p,sep="")
filename.csv= paste(path.results, model.name, "-alpha", alpha,".csv",sep="") 
    
# EM estimation of the model (results are later used for the starting values)
#-----------------------------------------------------------------------------------------------
EM.output <- EM.MSIAH(Y=Y, M=M, p=p,  max.iter=1000, convergence=0.0000001, diagonal.PR_TR=0.9, plot=FALSE)

# Prior distributions setup
#-----------------------------------------------------------------------------------------------
priors  <- list(Beta.bar        = array(0, c(K*(1+K*p),1)), # Mean of AR coeffs
                V.Beta.bar.m1   = solve(.3*diag(K*(1+K*p))),   # Variance of AR coeffs
                S.mean          = 0,    # Mean of lognorm standard deviation coeffs
                S.sd            = 2,    # SD of lognorm standard deviation coeffs
                w               = matrix(9 * diag(M) + matrix(1,M,M),ncol=1)    # Priors: Hidden Markov Chain, transition probabilities       
)

# Restrictions setup (unrestricted MS(2)-VAR(1))
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

# Gibbs algorithm: Burn-in sample
N.burnin        <- 10000
starting.values <- G.EM.to.Gibbs(model.type="MSIAH", EM.output, priors, restrictions, G) 
Burnin.output   <- G.MSBVAR.griddy(N=N.burnin, h=h, priors=priors, restrictions=restrictions, starting.values=starting.values, debug=FALSE, print.iterations=1000)

# Gibbs algorithm: Final sample
N.griddy        <- 100000     
t0 = proc.time()
Griddy.output   <- G.MSBVAR.griddy(N=N.griddy, h=h, priors=priors, restrictions=restrictions, starting.values=Burnin.output$last.draws, debug=FALSE, print.iterations=1000)
t1         = proc.time()
time.exe   = (t1-t0)/3600
cat("\ntime of execution: ",time.exe[3], " [h]")

save(Y, K, N.burnin, N.griddy, M, p, h, G, EM.output, Burnin.output, Griddy.output, priors, restrictions, file=paste(path.results, model.name, ".RData",sep=""))

# Marginal data density estimated with Modified Harmonic Mean estimator by Geweke (1999,2005)
MHM <- G.MHM.MSVAR(Griddy.output, priors, restrictions, alpha, debug=FALSE, print.iterations=250)

text.string <- paste("\nM: ",M,"\tp: ",p,"\tMHM: ", MHM$MHM, "\talpha: ", MHM$alpha, "\tacceptance.rate: ", MHM$acceptance.rate, " MC2", sep="")
cat(text.string,file=filename.csv,append=TRUE,sep="",fill=FALSE)  
   
save(Y, K, N.burnin, N.griddy, M, p, h, G, priors, restrictions, EM.output, Burnin.output, Griddy.output, MHM, time.exe, file=paste(path.results, model.name, "-MC2.RData",sep="")) #MC2
