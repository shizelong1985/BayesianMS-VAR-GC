
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Functions for computing information criterion
#---------------------------------------------------------------------------------------------------


# Wrapper function, calculating all of the criteria
get.IC <- function(model.type,K,p,M,T,max.log.likelihood,print.results=TRUE)    {
    
    if(model.type == "MSAH"){
        number.parameters   <- M * (M-1) + K * (1+M*K*p) + K*M * (K+1) / 2  
    }else if(model.type == "MSH"){
        number.parameters   <- M * (M-1) + K * (1+K*p) + K*M * (K+1) / 2
    }else if(model.type == "MSI"){
        number.parameters   <- M * (M-1) + K * (M+K*p) + K * (K+1) / 2
    }else if(model.type == "MSIA"){
        number.parameters   <- M * (M-1) + K*M * (1+K*p) + K * (K+1) / 2
    }else if(model.type == "MSIH"){
        number.parameters   <- M * (M-1) + K * (M+K*p) + K*M * (K+1) / 2
    }else if(model.type == "MSIAH"){
        number.parameters   <- M * (M-1) + K*M * (1+K*p) + K*M * (K+1) / 2 
    }else{
        stop("\n\n\t The model entered as parameter is not implemented.")
    }

    out.AIC <- get.AIC(T,max.log.likelihood,number.parameters)
    out.BIC <- get.BIC(T,max.log.likelihood,number.parameters)
    out.HQC <- get.HQC(T,max.log.likelihood,number.parameters)

    if(print.results == TRUE){
        cat("Model:", model.type, ", M=", M, ", p=", p, "\n")
        cat("Number of observations=",T , "\nNumber of parameters=",number.parameters, "\nMax log likelihood=", max.log.likelihood, "\n")
        cat("AIC=", out.AIC, "\tBIC=", out.BIC, "\tHQC=", out.HQC, "\n")
    }

    output  = list( AIC = out.AIC,
                    BIC = out.BIC,
                    HQC = out.HQC)

    return(output)
}


# Akaike information criterion (AIC) 
#---------------------------------------------------------------------------------------------------
 
get.AIC <- function(T,max.log.likelihood,number.parameters)    {
    ### out.AIC <- (max.log.likelihood - number.parameters)         # Psaradakis / Spagnolo
    out.AIC <- (-2*max.log.likelihood + 2*number.parameters) / T    # Krolzig
    return(out.AIC)
}


# Bayesian information criterion (BIC)
#---------------------------------------------------------------------------------------------------
 
get.BIC <- function(T,max.log.likelihood,number.parameters)    {
    ### out.BIC <- max.log.likelihood - 0.5 * log(T) * number.parameters    # Psaradakis / Spagnolo
    out.BIC <- (-2*max.log.likelihood + number.parameters * log(T)) / T     # Krolzig
    return(out.BIC)
}


# Hannan–Quinn criterion (HQC)
#---------------------------------------------------------------------------------------------------
 
get.HQC <- function(T,max.log.likelihood,number.parameters)    {
    ### out.HQC <- max.log.likelihood - log(log(T)) * number.parameters   # Psaradakis / Spagnolo
    out.HQC <- (-2*max.log.likelihood + 2 * number.parameters * log(log(T))) / T    # Krolzig
    return(out.HQC)
}


      
