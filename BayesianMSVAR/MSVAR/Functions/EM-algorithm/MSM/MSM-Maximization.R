
# Copyrighted 2016 Matthieu Droumaguet

# Contact:
#    Matthieu Droumaguet:   matthieu.droumaguet@gmail.com 

# The codes are available under the GNU General Public License v3.0. 
# To refer to the codes in publications, please, cite the following paper:
# Droumaguet, M., Warne, A., Woźniak, T. (in press) Granger Causality and Regime Inference in Markov-Switching VAR Models with Bayesian Methods, Journal of Applied Econometrics.



#---------------------------------------------------------------------------------------------------
#  MSM: M-step
#---------------------------------------------------------------------------------------------------
 
M.step.MSM  <- function(Y, p, M, sp, fp, mu, A, Sigma, F)    {

    T   <- dim(Y)[1] - p
    K   <- dim(Y)[2]
    N   <- M^(1+p)
    y   <- vec(t(Y[-(1:p),]))

    # Sigma in (K x K) format
    Sigma   <- Sigma[,,1]
    xi      <- t(sp)
     
    ### cat("\nM.step\n___________________________\n")

    # Storage
    alpha   <- matrix(NA, K*p, K)
    eta.t 	<- matrix(NA, N, T) 
    
    # stack A into alpha
    for(lag in 1:p){
        row.begin   <- (lag-1)*K + 1
        row.end     <- row.begin + K - 1
        alpha[row.begin:row.end,]   <- A[,,lag]
    }
    
    
    parameters  <- rbind(matrix(alpha,ncol=1),matrix(mu,ncol=1),matrix(Sigma,ncol=1)) 
    max.change  <- 1
    iter        <- 1
    while(max.change > 1e-3) {


        if((iter %% 10)==0) cat(iter," ")
        if(iter > 200) stop("Too many iterations in the maximization loop")

        ### cat("Sigma\n")
        ### print(Sigma)  
        ### cat("mu\n")
        ### print(mu)
        ### cat("alpha\n")
        ### print(alpha)



        #____________________________________________
        #
        out.mu      <- mu.tilda(Y,M,p,xi,Sigma,alpha)
        mu          <- out.mu$mu
        ### cat("OUTPUT mu\n")
        ### print(out.mu$mu.array)

        
        #____________________________________________
        #
        out.alpha   <- alpha.tilda(Y,M,p,xi,Sigma,mu)
        alpha       <- out.alpha$alpha
        ### cat("OUTPUT alpha\n")
        ### print(out.alpha$A)

        
        #____________________________________________
        #
        # Residuals
        U           <- residuals.max.step(Y, M, p, xi, alpha, mu)
        ### cat("OUTPUT residuals\n") 
        ### print(U)

        
        #____________________________________________
        #
        out.Sigma   <- sigma.tilda(U,p,M,xi)
        Sigma       <- out.Sigma
        ### cat("OUTPUT Sigma\n")
        ### print(Sigma)     

        
        #____________________________________________
        #
        # Convergence criteria: 
        parametersm1<- parameters
        parameters  <- rbind(matrix(alpha,ncol=1),matrix(mu,ncol=1),matrix(Sigma,ncol=1))
        change      <- parameters - parametersm1
        max.change  <- max(abs(change))
        ### cat("Iteration:", iter, "   Convergence criterion:", max.change, "\n")

        iter        <- iter + 1  
        ### print("Press enter")
        ### scan()     
    }


    cat("\t\tMax. step converged in", iter, "iterations\n")


    # Densities  
    #_____________________________________________________________________________
    #

    # Densities: Using R's normal density function >> faster
    for(regime in 1:N){
        res.tmp         <- matrix(U[,,regime], T, K)
        eta.tmp         <- dmvnorm(res.tmp, sigma=matrix(Sigma,ncol=K)) 
        eta.t[regime,]  <- matrix(eta.tmp, nrow=1)   
    }

    # likelihood  
    #_____________________________________________________________________________
    #  
    # xi.tp1.t and log.lik
    xi.tp1.t    <- Ergodic.PR_TR(F)
    log.lik     <- log(t(eta.t[,1]) %*% xi.tp1.t)
    for(t in 1:(T-1)){
        xi.tp1.t    <- t(F) %*% fp[,t]
        log.lik     <- log.lik +  log(t(eta.t[,t+1]) %*% xi.tp1.t)
    }   


    ### cat("Max. step Sigma\n")
    ### print(Sigma)        
    ### cat("Max. step mu\n")
    ### print(mu)
    ### cat("Max. step alpha\n")
    ### print(alpha)


    ### print("Press enter")
    ### scan()          

    # Sigma : (K x K x M) format for common expectation function
    #______________________________________________________________________________________________
    #
    sigma.out   <- array(NA, c(K,K,N))
    for(regime in 1:N){
     	sigma.out[,,regime] <- Sigma
    }    


	output <- list( mu          = out.mu$mu.array, 
                    A           = out.alpha$A, 
                    Sigma       = sigma.out, 
                    U           = U, 
                    likelihood  = log.lik )
    return(output)
} 






#-------------------------------------------------------------------------------------------------
mu.tilda    <- function(Y, M, p, xi, Sigma, alpha)  {
    
    cat("---------- mu.tilda IN \n")
    
    K       <- dim(Y)[2]
    T       <- dim(Y)[1] - p
    y       <- vec(t(Y[-(1:p),]))    # y, (TK x 1)
    A       <- array(NA, c(K, K, p))
    mu.array<- array(NA, c(K, 1, M))

    Sigma.inv   <- solve(Sigma)
    

    # Get A from stacked format
    for(j in 1:p){
        row.begin   <- (j-1)*K + 1
        row.end     <- row.begin + K - 1
        A[,,j]      <- alpha[row.begin:row.end,] 
    }

    N       <- M^(p+1)
    N.step  <- N / M

    LHS <- matrix(0, M*K, M*K)
    MHS <- matrix(0, M*K, T*K)


    for(m in 1:M){

        n.begin     <- (m-1) * N.step + 1
        n.end       <- n.begin + N.step - 1
        n.seq.m     <- n.begin:n.end 


        for(n in n.seq.m){
            T.n     <- drop(matrix(1,1,T) %*% xi[,n]) # drop() deletes dimensions
            X.n.mu  <- X.n.mu(A,p,M,n)
        
            tmp     <- t(X.n.mu) %*% Sigma.inv %*% X.n.mu
            LHS     <- LHS + tmp * T.n
            ### print(LHS)
           
            MHS     <- MHS + kronecker(t(xi[,n]),t(X.n.mu) %*% Sigma.inv)
            ### print(MHS)
        }
    }

    # Construction of the alpha vector = (alpha1', ..., alpham')'
    # !!!!! NOT SURE ABOUT taking transpose or not...
    alpha.vec   <- vec(alpha)
    # !!!!!!
    
    ### cat("alpha.vec\n")
    ### print(alpha.vec)  
    ### cat("kronecker(X.bar(Y,p), diag(K))\n")
    ### print(kronecker(X.bar(Y,p), diag(K)))  
    
    RHS <- y - kronecker(X.bar(Y,p), diag(K)) %*%  alpha.vec

    mu  <- solve(LHS) %*% MHS %*% RHS

    ### cat("mu\n")
    ### print(mu)

    for(regime in 1:M){
        row.begin           <- (regime-1) * K + 1
        row.end             <- row.begin + K - 1
        mu.array[,,regime]  <- mu[row.begin:row.end,]
    }
    

    cat("---------- mu.tilda OUT \n\n\n")



    output  <-list( mu      = mu,
                    mu.array= mu.array)
    return(output)
} 
   


#-------------------------------------------------------------------------------------------------
# RHS is problematic: need to take the series minus the intercept corresponding to the current regime. Smoothed or highest probability?
# Highest probability first!
alpha.tilda <- function(Y, M, p, xi, Sigma, mu) {

    cat("---------- alpha.tilda IN \n")

    K   <- dim(Y)[2]
    T   <- dim(Y)[1] - p
    y   <- vec(t(Y[-(1:p),]))    # y, (TK x 1)  
    A   <- array(NA, c(K, K, p))
    N       <- M^(p+1)
    N.step  <- N / M
    
    M.mat   <- t(mu) 

    # build the (TK x 1) matrix of intercepts for each observation
    highest.sp  <- max.col(xi)
    states      <- states(M,p)
    mu.m        <- matrix(NA, T*K, 1)
    for(t in 1:T){
        row.begin       <- (t-1)*K + 1
        row.end         <- row.begin + K - 1
        regime.base.n   <- highest.sp[t]
        regime.base.m   <- states[regime.base.n, 1]
        mu.begin        <- (regime.base.m-1)*K + 1
        mu.end          <- mu.begin -1 + K
 
        mu.m[row.begin:row.end, ]   <- mu[mu.begin:mu.end, ] 
    }

    ### print(highest.sp)
    ### print(mu.m)

    LHS <- matrix(0, K*p, K*p)
    MHS <- matrix(0, K*p, T)
 
    for(m in 1:M){

        n.begin <- (m-1) * N.step + 1
        n.end   <- n.begin + N.step - 1
        n.seq.m <- n.begin:n.end

        ### cat("M:",m,"\n")
        ### cat("n.seq.m:",n.seq.m,"\n")

        for(n in n.seq.m){
            ### T.n     <- matrix(1, 1, T) %*% xi[, n]
            Xi.n    <- diag(xi[, n])
            X.n.bar.star<- X.n.bar.star(Y, p, M, mu , n)
        
            tmp     <- t(X.n.bar.star) %*% Xi.n
            
            LHS     <- LHS +  tmp %*% X.n.bar.star
            MHS     <- MHS + tmp
        }
    }
    ### print(dim(LHS))
    ### print(dim(MHS))
    ### print(dim(y))
    ### print(mu)
    ### print(dim(kronecker(matrix(1,T,1), mu[1,])))
    ### print(dim(kronecker(solve(LHS) %*% MHS, diag(K))))
    ### RHS     <- y - kronecker(matrix(1,T,1), mu)

    RHS     <- y - mu.m     # kronecker notation does not work in here, because mu is regime m dependent

    alpha   <- kronecker(solve(LHS) %*% MHS, diag(K)) %*% RHS

    ### print(alpha)
    
    # alpha is a vector, convert it to a matrix
    alpha.mat   <- vec.inv(alpha,K)
    for(lag in 1:p){
        col.begin   <- (lag-1) * K + 1
        col.end     <- col.begin + K - 1
        A[,,lag]    <- alpha.mat[,col.begin:col.end]
    } 

    ### cat("alpha\n")
    ### print(alpha) 
    ### cat("A\n")
    ### print(A)
    cat("---------- alpha.tilda OUT \n\n\n")
 
    output  <- list(    alpha   = alpha,
                        A       = A)
    return(output)
}






#-------------------------------------------------------------------------------------------------
sigma.tilda <- function(U, p, M, xi)    {

    cat("---------- sigma.tilda IN \n") 

    T       <- dim(xi)[1]
    K       <- dim(U)[2]
    N       <- M^(p+1)
    N.step  <- N / M  
    Sigma   <- matrix(0,K,K) 

    ### print(U)

    for(m in 1:M){

        n.begin <- (m-1) * N.step + 1
        n.end   <- n.begin + N.step - 1
        n.seq.m <- n.begin:n.end
        
        ### cat("n.seq.m:", n.seq.m, "\n")    
        ### tmp <- matrix(1,1,T) %*% xi[,n.seq.m]
        ### T.m <- sum(tmp)

        for(n in n.seq.m){
            ### T.n     <- matrix(1,1,T) %*% xi[,n]
            Xi.n    <- diag(xi[,n])
            
            # which T.m ?????
            # Sigma being regime independent, I scale by T degrees of freedom.
            Sigma   <- Sigma +  T^(-1) * (t(U[,,n]) %*% Xi.n %*% U[,,n])
        }
    }
    
    cat("---------- sigma.tilda OUT \n\n\n")

    return(Sigma)
} 







            

#_________________________________________________________________________________________________
# returns residuals

residuals.max.step <- function(Y, M, p, xi, alpha, mu)  {

	cat("---------- residuals IN \n") 
    
    TT	<- dim(Y)[1]
    T   <- TT - p
    K	<- dim(Y)[2]
    N   <- M^(1+p)

    U   <- array(NA,c(T,K,N))
    
    states  <- states(M,p) 

    # alpha is a vector, change it to a matrix
    alpha.mat   <- vec.inv(alpha,K)
    ### print(alpha.mat)

    for(n in 1:N){
        regime.Y    <- states[n,1]
        mu.begin    <- (regime.Y-1)*K + 1
        mu.end      <- mu.begin -1 + K
        mu.Y        <- matrix(mu[mu.begin:mu.end,],ncol=1)

        X.bar.star  <- X.bar.star(Y,p,mu,states[n,-1])

        ### cat("X.bar.star\n")
        ### print(X.bar.star)
        ### cat("t alpha.mat\n")  
        ### print(dim(t(alpha.mat)))
        ### print(t(alpha.mat))
        ### cat("kronecker\n")  
        ### print(kronecker(matrix(1,T,1), t(mu.Y))) 

        ### print("Press enter")
        ### scan()      
        ### cat("dimensions\n")
        ### print(states[n,-1])
        ### print(dim(X.bar.star))
        ### print(dim(alpha))
        ### print(dim(t(alpha)))
        ### print(dim(kronecker(matrix(1,T,1), t(mu.Y))))

        U[,,n]  <- Y[-(1:p),] - X.bar.star %*% t(alpha.mat) - kronecker(matrix(1,T,1), t(mu.Y))
    }
          
    ### cat("U\n")
    ### print(U)
    ### print("Press enter")
    ### scan()          

    ### cat("---------- residuals OUT \n\n\n") 

    output 	<- U
	return(output)
}





#_________________________________________________________________________________________________
#



# L.j, (M x M^(1+p))
L.j <- function(M,p,j)  {
    
    ### cat("\t\t\t---- L.j IN \n")
 
    L.j <-  kronecker(kronecker(matrix(1,1,M^j),diag(M)),matrix(1,1,M^(p-j)))

    ### cat("\t\t\t---- L.j OUT \n")
     
    return(L.j)
}





# X.bar, (T x Kp)
X.bar   <- function(Y,p)    {
    
    ### cat("\t\t---- X.bar IN \n")

    TT      <- dim(Y)[1]
    T       <- TT - p
    K       <- dim(Y)[2]
    X_bar   <- matrix(0,T,(K*p))

    for(i in 1:p){
        col.begin   <- (i - 1) * K + 1
        col.end     <- col.begin - 1 + K

        row.end     <- TT - i
        row.begin   <- row.end - T + 1

        ### cat("col.begin", col.begin, "col.end", col.end, "row.begin", row.begin, "row.end", row.end, "\n")
        X_bar[,col.begin:col.end]   <- Y[row.begin:row.end,]
    }

    ### cat("\t\t---- X.bar OUT \n")
     
    return(X_bar)
}





# X.n.bar.star, (T x Kp)
X.n.bar.star    <- function(Y,p,M,mu,n) {
    # Where is n used exactly?????
    # in M.mat maybe???? History of regimes???
    # does not seem so

    # 1T X 1N' seems to be the problem: no kronecker.
    # We replace by a TxN matrix of 0 except the n th column of ones...
   
    # cat("\t---- X.n.bar.star IN \n")

    TT      <- dim(Y)[1]
    T       <- TT - p
    K       <- dim(Y)[2]
    N       <- M^(1+p)
    L       <- matrix(NA, N, M*p)

    X.bar   <- X.bar(Y,p)
    for(j in 1:p){
        col.begin               <- (j-1)*M + 1
        col.end                 <- col.begin + M - 1
        tmp                     <- L.j(M,p,j)
        L[,col.begin:col.end]   <- t(tmp)
    }

    M.mat   <- kronecker(diag(p), matrix(mu,nrow=K)) 

    # Could this BE ?
    EUREKA  <- diag(N)
    EUREKA  <- matrix(EUREKA[n,],nrow=1)
    # cat("EUREKA\n")
    # print(EUREKA)
    # print(kronecker(matrix(1,T,1),EUREKA))
    # print(L)
    # print(t(M.mat))
    # print(kronecker(matrix(1,T,1),EUREKA) %*% L %*% t(M.mat))
    # print("Press enter")
    # scan()        
   
    X.n.bar.star<- X.bar - kronecker(matrix(1,T,1),EUREKA) %*% L %*% t(M.mat)
    # !!!!!!

    # cat("\t---- X.n.bar.star OUT \n")
    
    return(X.n.bar.star)
}

 




# X.mu, (M^(1+p)K x MK)
X.mu    <- function(A,p,M)  {
    
    ### cat("\t---- X.mu IN \n")
 
    K       <- dim(A)[1]
    # j = 0
    X.mu    <- - kronecker(t(L.j(M,p,0)), -diag(K))    # seems one must take the transpose of L
                            
    for(j in 1:p){
        X.mu   <- X.mu - kronecker(t(L.j(M,p,j)), A[,,j])
    }

    ### cat("\t---- X.mu OUT \n")
  
    return(X.mu)
}  






# X.n.mu, (K x MK)
X.n.mu  <- function(A,p,M,n)    {
    
    ### cat("\t---- X.n.mu IN \n")
    
    K       <- dim(A)[1]
    N       <- M^(1+p)

    iota.trans  <- matrix(diag(N)[,n], nrow=1)
    tmp         <- kronecker(iota.trans,diag(K)) %*% X.mu(A,p,M)
    X.n.mu      <- matrix(tmp,nrow=K)   # X.n.mu can be a vector, force type to matrix
    
    ### cat("\t---- X.n.mu OUT \n")
    
    return(X.n.mu)
}






#_________________________________________________________________________________________________
# 

 
# X.bar.star, (T x Kp)
X.bar.star   <- function(Y,p,mu,states) {
    
    ### cat("\t---- X.bar.star IN \n")
    ### cat(states)

    TT      <- dim(Y)[1]
    T       <- TT - p
    K       <- dim(Y)[2]

    X.bar   <- X.bar(Y,p)

    ### RHS     <- matrix(NA,1,K*p)
    RHS     <- NULL
    for(j in 1:p){
        regime.Y    <- states[j]
        mu.begin    <- (regime.Y-1)*K + 1
        mu.end      <- mu.begin -1 + K
        mu.j        <- matrix(mu[mu.begin:mu.end,],ncol=1)   
        
        RHS         <- cbind(RHS,t(mu.j))
    }
    ### print(RHS)
    ### X.bar.star  <- X.bar - matrix(1,T,1) %*% kronecker(matrix(1,1,p),t(mu))
    ### cat("X.bar \n")
    ### print(X.bar)

    X.bar.star  <- X.bar - matrix(1,T,1) %*% RHS 

    ### print(dim(X.bar.star))
    ### print(X.bar.star)
   
    ### cat("\t---- X.bar.star OUT \n")
    
    return(X.bar.star)
}
