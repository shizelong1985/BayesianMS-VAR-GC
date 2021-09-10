
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Model selection based on information criteria
# For comparison of AIC, BIC, HQC within a given number of regimes
#---------------------------------------------------------------------------------------------------


# Altug & Bildirici (2010): 
#   Kapetanios (2001): AIC tends to choose longer lag lengths in MS-AR models whereas BIC tends to
#                      select more parsimonious ones.
#   Ivanov & Killian (2005): HQC is the most accurate criterion for selecting lag lengths in a
#                            quarterly VAR situation (metric is MSE).



# return IC results of model estimation for different lags, from 1 to p.max

Model.Selection.IC  <- function(model.type, Y, M, p.max, diag.PR_TR, convergence, max.iter, print, plot){
    Y   <- as.matrix(Y) # Cast into matrix, useful in unidimensional case

    # Storage: p rows, 3 columns (AIC, BIC, HQC)
    #--------------------------------------------------------------------------------------------
    IC.output   <- array(NA, c(p.max,3))


    for(p in 1:p.max){
        # Resize sample in order to estimate equally big samples over variations of lags
        # Estimate with EM
        #----------------------------------------------------------------------------------------
        ### cat("Estimation of model", p, "/", p.max, "\n")
        cat(p, "")
       
        to.remove   <- p.max - p  
        if(to.remove == 0){
            EM.output   <- MSVAR.EM(model.type, Y, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
        }else{
            Y.reduced   <- Y[-(1:to.remove),]
            EM.output   <- MSVAR.EM(model.type, Y.reduced, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
        }
        
        # Store results
        #----------------------------------------------------------------------------------------   
        if(EM.output$convergence.code == 0) {
           IC.output[p,1]   <- EM.output$IC$AIC
           IC.output[p,2]   <- EM.output$IC$BIC
           IC.output[p,3]   <- EM.output$IC$HQC
        }
    }

    # Model selection for each criterion
    #--------------------------------------------------------------------------------------------  
    min.AIC <- which.min(IC.output[,1])
    min.BIC <- which.min(IC.output[,2])
    min.HQC <- which.min(IC.output[,3])
    

    # Print results
    #--------------------------------------------------------------------------------------------  
    cat("\n\nLags \tAIC \t BIC \t HQC\n")
    print(IC.output)
    cat("\nModel selection:\n")
    cat("\tAIC, p=", min.AIC, "\n")
    cat("\tBIC, p=", min.BIC, "\n")
    cat("\tHQC, p=", min.HQC, "\n")
    

    # Output
    #-------------------------------------------------------------------------------------------- 
    output  <- list(IC.output   = IC.output,
                    min.AIC     = min.AIC,
                    min.BIC     = min.BIC,
                    min.HQC     = min.HQC)
    return(output)
}



#---------------------------------------------------------------------------------------------------
# Loop for different type of models and specifications
#---------------------------------------------------------------------------------------------------

IC.model.selection.loop <- function(models, Y, p.max, M,max.iter, convergence, diag.PR_TR, print, plot, series){

    number.models   <- length(models)

    output.table    <- array(NA, c(number.models * p.max, 6))

    cat("Series", series, "\n")
    for(model.no in 1:length(models)){
        model.type      <- models[model.no]
        print(model.type)

        out.selection   <- Model.Selection.IC(model.type=model.type,Y=Y, p.max=p.max, M=M,max.iter=max.iter, convergence=convergence, diag.PR_TR=diag.PR_TR, print=print, plot=plot)

        # row indexes
        row.begin   <- (model.no-1)*p.max + 1
        row.end     <- row.begin + p.max - 1
        # model name
        output.table[(row.begin:row.end),1]     <- model.type
        # models number of regimes
        output.table[(row.begin:row.end),2]     <- M
        # model lags
        output.table[(row.begin:row.end),3]     <- seq(1:p.max)
        # AIC, BIC, HQC
        output.table[(row.begin:row.end),4:6]   <- out.selection$IC.output 

    }
      
   filename <- paste("Model-selection-IC-series-", series, ".csv", sep="") 
   write.table(output.table, file = filename, sep = ",", col.names = c("Model", "regimes", "lags", "AIC", "BIC", "HQC"),row.names = FALSE)


   return(output.table)
}                              
