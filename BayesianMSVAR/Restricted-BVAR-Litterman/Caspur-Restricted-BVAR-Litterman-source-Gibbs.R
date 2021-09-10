
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
#install.packages("mvtnorm", dependencies=TRUE)
#install.packages("MCMCpack", dependencies=TRUE)
#install.packages("tmvtnorm", dependencies=TRUE)

library("mvtnorm")
library("MCMCpack")
library("tmvtnorm") 

                                                         
source("BayesianMSVAR/Restricted-BVAR-Litterman/Functions/G-Restricted-BVAR.R")
source("BayesianMSVAR/Restricted-BVAR-Litterman/Functions/G-Restricted-BVAR-Initialization.R")   
source("BayesianMSVAR/Restricted-BVAR-Litterman/Functions/G-Restricted-BVAR-Correlation.R")   
source("BayesianMSVAR/Restricted-BVAR-Litterman/Functions/G-Restricted-BVAR-Standard-Deviation.R")   
source("BayesianMSVAR/Restricted-BVAR-Litterman/Functions/G-Restricted-BVAR-Regression.R")
source("BayesianMSVAR/Restricted-BVAR-Litterman/Functions/G-Restricted-BVAR-Forecast.R")
source("BayesianMSVAR/Restricted-BVAR-Litterman/Functions/G-Restricted-BVAR-MHM-Geweke.R")
source("BayesianMSVAR/Restricted-BVAR-Litterman/Functions/G-Restricted-BVAR-MHM-SWZ.R")
source("BayesianMSVAR/Restricted-BVAR-Litterman/Functions/G-Restricted-BVAR-HMNR.R")
