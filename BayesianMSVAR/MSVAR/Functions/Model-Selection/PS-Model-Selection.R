
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Model selection based on the probability scores
#---------------------------------------------------------------------------------------------------

# QPS, LPS

# return PS results of model estimation for 1 to M regimes

Model.Selection.PS  <- function(model.type, prob.type, Y, M, p.max, diag.PR_TR, convergence, max.iter, NBER.data, print, plot){

    # Storage: p rows, 6 columns (QPS, LPS, WER, AIC, BIC, HQC)
    #--------------------------------------------------------------------------------------------
    PS.output   <- array(NA, c(p.max,6))


    for(p in 1:p.max){
        # Resize sample in order to estimate equally big samples over variations of lags
        # Estimate with EM
        #----------------------------------------------------------------------------------------
        ### cat("Estimation of model", p, "/", p.max, "\n")
        ### cat("Lags: p=", p, "")
        cat(p, " ")

        to.remove       <- p.max - p
        start.Y.resize  <- c(start(Y)[1], start(Y)[2]+to.remove)

        if(to.remove == 0){
            Y.tmp       <- as.matrix(Y)
            Y.resize    <- Y.tmp
        }else{
            Y.tmp       <- as.matrix(Y)
            Y.resize    <- Y.tmp[-(1:to.remove),]
        }
        EM.output   <- MSVAR.EM(model.type, Y.resize, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
        
        # Compute LPS, QPS
        start.date.sp <- c(start.Y.resize[1], start.Y.resize[2]+p)
        Xi            <- NBER.data[-(1:p.max)]  # NBER.data series are equal to 1 when a recessions occurs

        out.QPS.LPS.WER <- QPS.LPS.WER.plot(prob.type, EM.output, Xi, start.date.sp, frequency=frequency(Y), print=print)           

        # Store results
        #----------------------------------------------------------------------------------------   
        if(EM.output$convergence.code == 0) {
           PS.output[p,1]   <- EM.output$IC$AIC
           PS.output[p,2]   <- EM.output$IC$BIC
           PS.output[p,3]   <- EM.output$IC$HQC   
           PS.output[p,4]   <- out.QPS.LPS.WER$qps
           PS.output[p,5]   <- out.QPS.LPS.WER$lps
           PS.output[p,6]   <- out.QPS.LPS.WER$wer
        }
    }

    # Model selection for each criterion
    #--------------------------------------------------------------------------------------------  
    min.AIC <- which.min(PS.output[,1])
    min.BIC <- which.min(PS.output[,2])
    min.HQC <- which.min(PS.output[,3])
    min.QPS <- which.min(PS.output[,4])
    min.LPS <- which.min(PS.output[,5])
    min.WER <- which.min(PS.output[,6])  
    

    # Print results
    #--------------------------------------------------------------------------------------------  
    if(print){
        cat("\n\nLags \tQPS \tLPS \tWER \tAIC  \tBIC  \tHQC \n")
        print(PS.output)
        cat("\nModel selection:\n")
        cat("\tQPS, p=", min.QPS, "\n")
        cat("\tLPS, p=", min.LPS, "\n")
        cat("\tWER, p=", min.WER, "\n")
        cat("\tAIC, p=", min.AIC, "\n")
        cat("\tBIC, p=", min.BIC, "\n")
        cat("\tHQC, p=", min.HQC, "\n")
    }

    # Output
    #-------------------------------------------------------------------------------------------- 
    output  <- list(PS.output   = PS.output,
                    min.AIC     = min.AIC,
                    min.BIC     = min.BIC,
                    min.HQC     = min.HQC,
                    min.QPS     = min.QPS,
                    min.LPS     = min.LPS,
                    min.WER     = min.WER)
    return(output)
}







#---------------------------------------------------------------------------------------------------
# LOOP with series name
#---------------------------------------------------------------------------------------------------
 
loop.model.selection.name    <- function(models, prob.type, Y, p.max, M,max.iter, convergence, diag.PR_TR, NBER.data, print, plot, series){

    number.models   <- length(models)

    output.table    <- array(NA, c(number.models * p.max, 9))

    cat("\nSeries", series, "\n")
    for(model.no in 1:length(models)){
        model.type      <- models[model.no]
        print(model.type)

        out.selection   <- Model.Selection.PS(model.type=model.type, prob.type=prob.type, Y=Y, p.max=p.max, M=M,max.iter=max.iter, convergence=convergence, diag.PR_TR=diag.PR_TR, NBER.data=NBER.data, print=print, plot=plot)

        # row indexes
        row.begin   <- (model.no-1)*p.max + 1
        row.end     <- row.begin + p.max - 1
        # model name
        output.table[(row.begin:row.end),1]     <- model.type
        # models number of regimes
        output.table[(row.begin:row.end),2]     <- M
        # model lags
        output.table[(row.begin:row.end),3]     <- seq(1:p.max)
        # AIC, BIC, HQC, QPS, LPS, WER
        output.table[(row.begin:row.end),4:9]   <- out.selection$PS.output 
    }
      
    # Write results into CSV
    filename <- paste("MS-", series, "-M", M, "-", prob.type, ".csv", sep="") 
    write.table(output.table, file = filename, sep = ",", col.names = c("Model", "regimes", "lags", "AIC", "BIC", "HQC", "QPS", "LPS", "WER"),row.names = FALSE)


    return(output.table)
}




#---------------------------------------------------------------------------------------------------
# LOOP
#---------------------------------------------------------------------------------------------------

### loop.model.selection    <- function(models, prob.type, Y, p.max, M,max.iter, convergence, diag.PR_TR, NBER.data, print, plot){

    ### number.models   <- length(models)

    ### output.table    <- array(NA, c(number.models * p.max, 9))
   ###  
    ### for(model.no in 1:length(models)){
        ### model.type      <- models[model.no]
        ### out.selection   <- Model.Selection.PS(model.type=model.type, prob.type=prob.type, Y=Y, p.max=p.max, M=M,max.iter=max.iter, convergence=convergence, diag.PR_TR=diag.PR_TR, NBER.data=NBER.data, print=print, plot=plot)

        ### # row indexes
        ### row.begin   <- (model.no-1)*p.max + 1
        ### row.end     <- row.begin + p.max - 1
        ### # model name
        ### output.table[(row.begin:row.end),1]     <- model.type
        ### # models number of regimes
        ### output.table[(row.begin:row.end),2]     <- M
        ### # model lags
        ### output.table[(row.begin:row.end),3]     <- seq(1:p.max)
        ### # AIC, BIC, HQC, QPS, LPS, WER
        ### output.table[(row.begin:row.end),4:9]   <- out.selection$PS.output 

    ### }
     ###  
   ### # Save results in CSV file
   ### write.table(output.table, file = "Model-selection.csv", sep = ",", col.names = c("Model", "regimes", "lags", "AIC", "BIC", "HQC", "QPS", "LPS", "WER"),row.names = FALSE)


   ### return(output.table)
### }
