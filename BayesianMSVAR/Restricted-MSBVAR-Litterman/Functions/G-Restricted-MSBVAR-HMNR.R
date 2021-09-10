
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.




G.HMNR.MSVAR = function(likeli){
      ###########################################################################
      # The function computes the marginal density of data from the MCMCM output
      # using Newton, Raftery (1994) Harmonic Mean estimator	
      ###########################################################################
      # Inputs: likeli - vector of the values of the log-likelihood functions evaluated at the draws from the posterior distribution
      ###########################################################################
      # Output:	ln(p(y|M)) - logarithm of the marginal density of data
      ###########################################################################
      
      # Computing the value of the marginal density of data as in eq. (15) of SWZ
      constant = mean(likeli) 
      estimate = 1/(mean(exp( - likeli + constant)))
            
      ml.log = log(estimate) + constant
      
      # TO DO: Numerical Standard error (How to do that?)
      
      return(ml.log)
}