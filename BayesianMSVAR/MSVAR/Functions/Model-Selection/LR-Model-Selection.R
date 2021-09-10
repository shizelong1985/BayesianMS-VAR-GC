
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Model selection based on the Likelihood Ratio test
#---------------------------------------------------------------------------------------------------

# About likelihood ratio test: Cameron, Trivedi textbook, p. 234

# Altug & Bildirici (2010): 
#   Ang & Bekaert (2002): 


# return LR results of model estimation for 1 to M regimes

Model.Selection.LR.2vs3  <- function(model.type, Y, p.max, diag.PR_TR, convergence, max.iter, print, plot){

    # Storage: Rows: p. Cols: maxL2, maxlL3, LR (2vs3), corrected degrees of freedom
    #--------------------------------------------------------------------------------------------
    LR.output   <- array(NA, c(p.max,4))


    for(p in 1:p.max){
        # Resize sample in order to estimate equally big samples over variations of lags
        # Estimate with EM
        #----------------------------------------------------------------------------------------
        ### cat("Estimation of model", p, "/", p.max, "\n")
        cat(p, "")

        to.remove       <- p.max - p  
        start.Y.resize  <- c(start(Y)[1], start(Y)[2]+to.remove) 

        if(to.remove == 0){
            Y.tmp       <- as.matrix(Y)
            Y.resize    <- Y.tmp
        }else{
            Y.tmp       <- as.matrix(Y)
            Y.resize    <- Y.tmp[-(1:to.remove),]
        }
        for(M in c(2,3)){
            EM.output   <- MSVAR.EM(model.type, Y.resize, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
            # Store results
            max.lik     <- EM.output$likelihood[length(EM.output$likelihood)]
            if(M==2){
                LR.output[p,1] <- max.lik
            }else if(M==3){
                LR.output[p,2] <- max.lik
            }
        }
        LR.output[p,3]  <- 2*(LR.output[p,2] - LR.output[p,1])
        # Corrected degrees of freedom
            # p11+p12=1, p21+p22=1
            c.df    <- 2
            # nuisances unidentified under the null
            # p31, p32
            c.df    <- c.df + 2
            # regime 3 parameters
            K       <- dim(Y.tmp)[2]
            if(model.type == "MSAH"){
                r.3.parameters   <- K^2 * p + K*(K+1) / 2  
            }else if(model.type == "MSH"){
                r.3.parameters   <- K * (K+1)/ 2
            }else if(model.type == "MSI"){
                r.3.parameters   <- K
            }else if(model.type == "MSIA"){
                r.3.parameters   <- K * (1+K*p)
            }else if(model.type == "MSIH"){
                r.3.parameters   <- K + K * (K+1) / 2
            }else if(model.type == "MSIAH"){
                r.3.parameters   <- K * (1+K*p) + K * (K+1) / 2 
            }else{
                stop("\n\n\t The model entered as parameter is not implemented.")
            }
        c.df    <- c.df + r.3.parameters
            
        LR.output[p,4]  <- c.df
      
    }
    

    # Output
    #-------------------------------------------------------------------------------------------- 
    output  <- list(LR.output   = LR.output)
    return(output)
}
 



#---------------------------------------------------------------------------------------------------
# Loop for different type of models and specifications
#---------------------------------------------------------------------------------------------------

LR.model.selection.loop <- function(models,Y, p.max, max.iter, convergence, diag.PR_TR, print, plot, series){

    number.models   <- length(models)

    output.table    <- array(NA, c(number.models * p.max, 5))

    cat("Series", series, "\n")
    for(model.no in 1:length(models)){
        model.type      <- models[model.no]
        print(model.type)

        out.selection   <- Model.Selection.LR.2vs3(model.type=model.type,Y=Y, p.max=p.max,max.iter=max.iter, convergence=convergence, diag.PR_TR=diag.PR_TR, print=print, plot=plot)

        # row indexes
        row.begin   <- (model.no-1)*p.max + 1
        row.end     <- row.begin + p.max - 1
        # model name
        output.table[(row.begin:row.end),1]     <- model.type
        # LR maxL2, maxlL3, LR (2vs3), corrected degrees of freedom 
        output.table[(row.begin:row.end),2:5]   <- out.selection$LR.output

    }
      
   # Save results in CSV file
   filename <- paste("Model-selection-LR-series-", series, ".csv", sep="") 
   write.table(output.table, file = filename, sep = ",", col.names = c("Model", "2 regimes max. lik.", "3 regimes max. lik.", "LR", "df"),row.names = FALSE)


   return(output.table)
}
