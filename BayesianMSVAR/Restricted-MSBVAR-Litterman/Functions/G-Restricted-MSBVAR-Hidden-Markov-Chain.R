
# Copyrighted 2016 Matthieu Droumaguet and Tomasz Woźniak

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 
#    Tomasz Woźniak:        wozniak.tom@gmail.com

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
# Gibbs: hidden Markov step, drawing from the Dirichlet distribution
# Restricted matrix of transition probabilities, P, Sims et al.(2008)
# 
# Metropolis-Hastings algorithm for stationary Markov Chains (Fruehwirth-Schnatter 2006, p.341)
# For stationary Markov chains, the initial distribution is equal to the ergodic distribution. The rows of the transition matrix are not independent a posteriori. One can use a Metropolis-Hastings algorithm.
#---------------------------------------------------------------------------------------------------
  
G.hidden.markov.chain <- function(aux, priors, restrictions){
    
    # Setup constants 
    M   <- dim(aux$PR_TR)[1]
    p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    T   <- dim(aux$Y)[1]
    TT  <- T - p
    K   <- dim(aux$Y)[2]

    Q   <- length(restrictions$dj)
    d   <- sum(restrictions$dj)


    PR_TR.old   <- aux$PR_TR

    # count transitions from xi 
    transitions.matrix  <- count.regime.transitions(aux$xi)
    transitions.vector  <- matrix(transitions.matrix, ncol=1)

    w.new           <- array(NA, c(d,1))
    wj.transitions  <- array(NA, c(d,1)) 
    for(j in 1:Q){
        index.begin <- sum(restrictions$dj[1:j-1]) + 1
        dj          <- restrictions$dj[j]
        index.end   <- index.begin + dj - 1
        
        for(i in index.begin:index.end){
            wj.transitions[i]   <- matrix(restrictions$M[,i],nrow=1) %*% transitions.vector 
        }

        alpha       <- wj.transitions[index.begin:index.end] + priors$w[index.begin:index.end]
        w.new[index.begin:index.end]    <- rdirichlet(n=1, alpha=alpha)
    }
    

    PR_TR.new   <- matrix(restrictions$M %*% w.new, nrow=M,byrow=TRUE)

    # Ergodic probabilities
    ergodic.old <- G.Ergodic.PR_TR(PR_TR.old)
    ergodic.new <- G.Ergodic.PR_TR(PR_TR.new)
    # Draw S0 states from previous ergodic probabilities:
    S0 <- sample(1:M, 1, prob=ergodic.old)
	
    # MH acceptance ratio
    AA  <- ergodic.new[S0] / ergodic.old[S0]

     # Metropolis-Hastings acceptance or rejection
    if(runif(1) <= AA ){
        aux$PR_TR   <- PR_TR.new
		aux$w		<- w.new
    }else{
        aux$PR_TR   <- PR_TR.old
    }

    # Output 
    return(aux)

}




#---------------------------------------------------------------------------------------------------
# Counts the transitions between regimes, from xi 
#---------------------------------------------------------------------------------------------------

count.regime.transitions <- function(xi){
	M   <- dim(xi)[1]
    TT  <- dim(xi)[2]

    count  <- matrix(0, M, M) 

    s       <- max.col(t(xi))
    lag.s   <- diff(s)

    for(row in 1:M){
        in.regime <- which(s == row) - 1 # -1 because we work with lags
        transitions <- lag.s[in.regime]
        for(col in 1:M){
            look.at <-  row - col
            count[row,col]  <- sum(transitions == look.at)
        }
    }

    return(count)
}      

 #---------------------------------------------------------------------------------------------------
# Gibbs: hidden Markov step, drawing from the Dirichlet distribution
# Metropolis-Hastings algorithm for stationary Markov Chains (Fruehwirth-Schnatter 2006, p.341)
# For stationary Markov chains, the initial distribution is equal to the ergodic distribution. The rows of the transition matrix are not independent a posteriori. One can use a Metropolis-Hastings algorithm.
#---------------------------------------------------------------------------------------------------

### G.hidden.markov.dirichlet.priors <- function(aux, priors){
   ###  
    ### # Setup constants 
    ### M   <- dim(aux$PR_TR)[1]
    ### p   <- (dim(aux$Y)[1] - dim(aux$U)[1])
    ### T   <- dim(aux$Y)[1]
    ### TT  <- T - p
    ### K   <- dim(aux$Y)[2]   

    ### PR_TR.old   <- aux$PR_TR
       ###  
    ### # count transitions from xi 
    ### transitions <- count.regime.transitions(aux$xi)
    ### # incorporate priors
    ### alpha       <- transitions + priors$alpha
    ### # Draw from Dirichlet distribution      
    ### PR_TR.new   <- t(apply(alpha, 1, rdirichlet, n=1))  

    ### # Ergodic probabilities
    ### ergodic.old <- G.Ergodic.PR_TR(PR_TR.old)
    ### ergodic.new <- G.Ergodic.PR_TR(PR_TR.new)
    ### # Draw S0 states:
    ### S0 <- sample(1:M, 1, prob=ergodic.new)
    ### # MH acceptance ratio
    ### AA  <- ergodic.new[S0] / ergodic.old[S0]

     ### # Metropolis-Hastings acceptance or rejection
    ### if(runif(1) <= AA ){
        ### aux$PR_TR   <- PR_TR.new
    ### }else{
        ### aux$PR_TR   <- PR_TR.old
    ### }

    ### # Output 
    ### return(aux)

### }







#---------------------------------------------------------------------------------------------------
# Counts the transitions between regimes, from xi 
#---------------------------------------------------------------------------------------------------

### count.regime.transitions <- function(xi){

	### M   <- dim(xi)[1]
    ### TT  <- dim(xi)[2]
   ###  
    ### s   <- max.col(t(xi))

    ### count  <- matrix(0, M, M)
    ### for (t in 2:TT){ 
        ### count[s[t-1],s[t]]  <- count[s[t-1], s[t]] + 1
    ### }
 ###  
	### return(count)
### } 
