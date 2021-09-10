
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Gibbs algorithms
#---------------------------------------------------------------------------------------------------

# Libraries
#-----------------------------------------------------------------------------------------------
### install.packages("mvtnorm", dependencies=TRUE)
### install.packages("MCMCpack", dependencies=TRUE)
### install.packages("tmvtnorm", dependencies=TRUE)
### install.packages("DEoptim", dependencies=TRUE)
library("mvtnorm")
library("MCMCpack")
library("tmvtnorm")
library("DEoptim")

# Functions
#-----------------------------------------------------------------------------------------------
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR.R")    
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-Filtering-Smoothing.R")
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-Hidden-Markov-Chain.R")
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-Standard-Deviation.R")
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-Correlation.R") 
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-Regression.R")
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-Forecast.R")
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-Initialization.R")   
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-MHM.R")
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-MHM-reducedRankP.R")
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-HMNR.R")   
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-other-Rank2.R")
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-other-A3-rank1.R")   
source("BayesianMSVAR/Restricted-MSBVAR-Litterman/Functions/G-Restricted-MSBVAR-other-A3-rank2.R")
