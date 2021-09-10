
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Wrapper: call all EM algos depending on the selected type
#---------------------------------------------------------------------------------------------------

MSVAR.EM <- function(model.type, Y, M, p, max.iter=100, convergence=0.0001, diagonal.PR_TR=0.5, print=TRUE, plot=TRUE) {

    if(model.type == "MSAH") {
        EM.output   <- EM.MSAH(Y, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
    }else if(model.type == "MSH") {
        EM.output   <- EM.MSH(Y, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
    }else if(model.type == "MSI") {
        EM.output   <- EM.MSI(Y, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
    }else if(model.type == "MSIAH") {
        EM.output   <- EM.MSIAH(Y, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
    }else if(model.type == "MSIA") {
        EM.output   <- EM.MSIA(Y, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
    }else if(model.type == "MSIH") {
        EM.output   <- EM.MSIH(Y, M, p, max.iter=max.iter, convergence=convergence, diagonal.PR_TR=diag.PR_TR, print=print, plot=plot)
    }else{
        stop("\n\n\t The model entered as parameter is not implemented.")
    }

    return(EM.output)
}
