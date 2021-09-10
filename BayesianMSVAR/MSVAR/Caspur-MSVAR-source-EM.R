
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# EM algorithms
#---------------------------------------------------------------------------------------------------



# Libraries
#---------------------------------------------------------------------------------------------------
#install.packages("mvtnorm", dependencies=TRUE)

library("mvtnorm")



# Functions common to all MSVAR models
#---------------------------------------------------------------------------------------------------
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSVAR-Wrapper.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSVAR-Expectation.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSVAR-Information-Criteria.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSVAR-QPS-LPS-WER.R") 


source("BayesianMSVAR/MSVAR/Functions/Model-Selection/IC-Model-Selection.R") 
source("BayesianMSVAR/MSVAR/Functions/Model-Selection/PS-Model-Selection.R") 
source("BayesianMSVAR/MSVAR/Functions/Model-Selection/LR-Model-Selection.R") 


                                                         
# Model specific functions
#---------------------------------------------------------------------------------------------------

# MSAH 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSAH-Initialization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSAH-Maximization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSAH-VAR.R") 
 
 
# MSH 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSH-Initialization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSH-Maximization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSH-VAR.R") 


# MSI 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSI-Initialization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSI-Maximization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSI-VAR.R") 


# MSIA  
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSIA-Initialization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSIA-Maximization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSIA-VAR.R") 


# MSIAH 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSIAH-Initialization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSIAH-Maximization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSIAH-VAR.R") 


# MSIH 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSIH-Initialization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSIH-Maximization.R") 
source("BayesianMSVAR/MSVAR/Functions/EM-algorithm/MSIH-VAR.R") 
