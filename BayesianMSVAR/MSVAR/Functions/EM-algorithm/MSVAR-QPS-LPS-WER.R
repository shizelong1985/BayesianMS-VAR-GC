
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.




# Wrapper for  stats + plot
#---------------------------------------------------------------------------------------------------
# Choose QPS, LPS, WER computation either with smoothed or filtered probabilities
QPS.LPS.WER.plot  <- function(prob.type, EM.output, Xi, start.date, frequency, print=TRUE){

    if(prob.type == "Sp"){
        Est.P   <- EM.output$sp
    }else if(prob.type == "Fp"){
        Est.P   <- EM.output$fp
    }

    p   <- dim(EM.output$Y)[1] -  dim(Est.P)[1]

    T.Est.P <- dim(Est.P)[1]
    T.Xi    <- length(Xi)
    if(T.Est.P != T.Xi) { stop("\n\n\t -- Length of regime vector must equal the dimension of the series after estimation, i.e. think about removing p lags to the regime vector.") }  


    # plot smoothed probabilities against known Xi
    #-----------------------------------------------------------------------------------------------
    op <- par(mfrow = c(M,1), mar=c(2,1,3,0.5), oma=c(1.5,2,1,1)) # M graphs on one plot, and margins settings

    # Plot of smoothed probabilities
    for(regime in 1:M){
        ### series.plot <- ts(cbind(Est.P[,regime], Xi),  start=start.date, frequency=frequency)
        ### title   <- paste("Estimated probabilities, regime ", regime, sep="")
        ### plot(series.plot, main=title ,ylim=c(0, 1), plot.type="single",col=c(1,4,2),lty=c(1,3,2))
        series.plot <- ts(cbind(Est.P[,regime], Xi),  start=start.date, frequency=frequency)
        title   <- paste("Estimated probabilities, regime ", regime, sep="")
        plot(series.plot, main=title ,ylim=c(0, 1), plot.type="single",col=c(1,2),lty=c(1,2))  
    }
    # par(ask=TRUE)
    par(op) # At end of plotting, reset to previous settings


    # QPS, LPS, WER
    #-----------------------------------------------------------------------------------------------
    out.QPS <- QPS.mapping(Est.P, Xi)
    regime1 <- out.QPS$regime1
    qps     <- out.QPS$qps
    lps     <- LPS.post.mapping(Est.P, Xi, regime1)
    wer     <- WER.post.mapping(Est.P, Xi, regime1)

    if(print==TRUE) cat("QPS =", qps, " , LPS =", lps, " , WER =", wer, "\n")


    # Output
    #-----------------------------------------------------------------------------------------------
    output  <- list(    qps = qps,
                        lps = lps,
                        wer = wer)
    return(output)
}








#---------------------------------------------------------------------------------------------------
# QPS and LPS with detection of regime 1 in QPS function, plus use of this info in LPQ function
#---------------------------------------------------------------------------------------------------

# Performs mapping of regime 1 in Xi in the smoothed probabilities of the EM estimated model
# Mapping corresponds to the lowest QPS calculated over all the estimated regimes
#---------------------------------------------------------------------------------------------------
QPS.mapping     <- function(Est.P, Xi)   {
    
    T   <- length(Xi)
    M   <- dim(Est.P)[2]

    qps <- array(NA,c(1,M))

    # St: Dummy taking the value 1 if the regime is the regime-th one
    St                 <- matrix(0,1,T)
    St[which(Xi==1)]   <- 1

    for(regime in 1:M){
        aaa             <- (Est.P[,regime] - St[,])^2 
        qps[,regime]    <- sum(aaa) * 2 / T
    }
    
    # Select which estimated regime produces less QPS when compared to regime 1 in Xi
    best    <- which.min(qps)
    ### print(qps)
    ### print(best)
    
    output  <- list(regime1 = best,
                    qps     = qps[,best])
    return(output)
}


#---------------------------------------------------------------------------------------------------
LPS.post.mapping     <- function(Est.P, Xi, regime1)   {

    T   <- length(Xi)
    
    # St: Dummy taking the value 1 if the regime is the first one
    St                  <- matrix(0,T,1)
    St[which(Xi==1),1]  <- 1
    
    ### sum.st.1        <- (1-St) * apply(matrix(Est.P[,regime1],ncol=1), 1, function(x) log(max(1-x,0)))
    sum.st.1        <- (1-St) * apply(matrix(Est.P[,regime1],ncol=1), 1, function(x) log(1-x))
    ### sum.st.other    <- St * apply(matrix(Est.P[,regime1],ncol=1), 1, function(x) log(min(x,1))) 
    sum.st.other    <- St * apply(matrix(Est.P[,regime1],ncol=1), 1, function(x) log(x)) 
    sum.st.1        <- replace(sum.st.1, which(sum.st.1=="NaN"),0)
    sum.st.other    <- replace(sum.st.other, which(sum.st.other=="NaN"),0)
    sum.both        <- sum.st.1 + sum.st.other

    lps <- -1 * sum(sum.both) / T  
    return(lps)
}
     


#---------------------------------------------------------------------------------------------------
# Percentage of wrongly estimated regimes, more intuitive measure than LPS and QPS
#---------------------------------------------------------------------------------------------------
WER.post.mapping     <- function(Est.P, Xi, regime1)   {

    T   <- length(Xi)

    # St: Dummy taking the value 1 if the regime is the first one
    St                  <- matrix(0,T,1)
    St[which(Xi==1),1]  <- 1

    # Estimated regimes
    Est.Xi  <- max.col(Est.P)
    # Correct for mapping: construct dummy where Est.Xi = 1 is estimated regime is regime1
    Est.St                          <- matrix(0,T,1)
    Est.St[which(Est.Xi==regime1),1]<- 1

    # Compute stat
    Xi.diff         <- St - Est.St
    Xi.deviations   <- 100 * sum(Xi.diff != 0) / T  

    wer             <- Xi.deviations
    return(wer)
}
     
     



 




#---------------------------------------------------------------------------------------------------
# Quadratic and Logarithmic Probability losses 

# QPS and LPS here only measure when the regime 1 is correctly identified, not the accuracy in the estimation of other regimes!!   
#---------------------------------------------------------------------------------------------------

QPS     <- function(Est.P, Xi)   {
    
    T   <- length(Xi)
    
    # St: Dummy taking the value 1 if the regime is the first one
    St                  <- matrix(0,1,T)
    St[which(Xi==1)]    <- 1

    aaa <- (Est.P[1,] - St[,])^2 
    
    qps <- sum(aaa) * 2 / T
    return(qps)
}





LPS     <- function(Est.P, Xi)   {
    
    T   <- length(Xi)
    
    # St: Dummy taking the value 1 if the regime is the first one
    St                  <- matrix(0,T,1)
    St[which(Xi==1)]    <- 1
    
    ### sum.st.1        <- (1-St) * apply(matrix(Est.P[,1],ncol=1), 1, function(x) log(max(1-x,0)))
    sum.st.1        <- (1-St) * apply(matrix(Est.P[,1],ncol=1), 1, function(x) log(1-x))
    ### sum.st.other    <- St * apply(matrix(Sp[,1],ncol=1), 1, function(x) log(min(x,1))) 
    sum.st.other    <- St * apply(matrix(Est.P[,1],ncol=1), 1, function(x) log(x)) 
    sum.st.1        <- replace(sum.st.1, which(sum.st.1=="NaN"),0)
    sum.st.other    <- replace(sum.st.other, which(sum.st.other=="NaN"),0)
    sum.both        <- sum.st.1 + sum.st.other

    lps             <- -1 * sum(sum.both) / T
    return(lps)
}

