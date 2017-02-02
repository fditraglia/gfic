#-----------------------------------------------------------------#
#SCRIPT:  GFIC_simulation_v1.R
#AUTHOR:  Frank DiTraglia
#EMAIL:   fditra@sas.upenn.edu
#DATE:    February 7th, 2013
#-----------------------------------------------------------------#
#DESCRIPTION:   
#
#     This script contains functions to carry out a simulation studying the properties GFIC model and moment selection in a linear, dynamic panel model. The correct model is:
#
#     y[i,t] = gamma * y[i,t-1] + theta * x[i,t] + eta[i] + v[i,t]
#
#The random variables x[i,t], eta[i] and v[i,t] are normal with mean zero and constant variances: s.x.sq, s.eta.sq, s.v.sq. The covariance between x[i,t] and eta[i] is s.x.eta; the covariance between x[i,t] and v[i,t-1] is s.x.v while that between x[i,t] and v[i,k] for (k != t-1) is zero. That is, x is weakly exogenous (predetermined) but is correlated with the one-period lagged error.
#     The purpose of this simulation is to calculate the RMSE of four estimators of the paramters for varying values of gamma and s.x.v. The four estimators are as follows:
#
#     Name        Lagged y?     Assume x Strictly Exogenous?
#     --------    ---------     ----------------------------
#     b.LW        YES           NO (assume weakly exogenous)
#     b.LS        YES           YES
#     b.W         NO            NO (assume weakly exogenous)
#     b.S         NO            YES
#
#The first two estimators, b.LW and b.LS, use the correct model and hence estimate (gamma, theta). The second two estimators omit the lagged y, estimating only theta.
#     The idea here is to compare bias versus variance. By omitting the lagged effect, we gain a time period for estimation and have one less paramter to estimate; by assuming strict exogeneity, we gain additional instruments. Both of these decrease the variance of the estimator. However, if the assumptions are wrong they introduce bias. The question is, how small must gamma be for it to make sense to omit the dynamic effect from an RMSE perspective? Similarly, how small must s.x.v be for assuming strict exogneity to make sense?
#     To answer this question, this script compares the RMSE of the 2SLS estimator for b.LW, b.LS, b.W and b.S. The simulation setup described above is similar to that of Andrews and Lu (2001) with two main differences. First, I do not include a constant term. Second, I do not create extra pre-sample observations to accomodate the lagged y. That is, in my simulation omitting the lag gives an extra time period. This is important because in a short panel additional time periods could considerably reduce the variance of an estimator. Third, I set the initial value of y, y.0, equal to the mean of its stationary distribution (zero) as in Arellano and Bond (1991), unlike Andrews and Lu. Fourth, I use 2SLS and only strong instruments rather than all instruments available under the assumptions. My estimator is essentially the Anderson and Hsiao Estimator (1982).
#-----------------------------------------------------------------#


#Load packages used below
library(MASS)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("/Users/fditra/Dropbox/PhD_dissertation/Generalized FIC Paper/simulations/functions/fast2SLS_v2.cpp")

#-----------------------------------------------------------------#
# Test Parameter Values (for debugging)
#-----------------------------------------------------------------#
# gamma <- 0.85
# r.x.v <- 0.5
# N.sims <- 2000
# N.i <- 500
# N.t <- 5
# theta <- 0.5
# s.x.eta <- 0.2
# s.eta.sq <- 1
# s.v.sq <- 1
# s.x.sq <- 1
#-----------------------------------------------------------------#

#-----------------------------------------------------------------#
# FUNCTION TO CARRY OUT THE SIMULATION A FIXED PARAMETER VALUES
#-----------------------------------------------------------------#
RMSE.sim.FAST <- function(gamma, r.x.v, N.sims = 1000, 
                          N.i = 250, N.t = 4, theta = 0.5, s.x.eta = 0.2, 
                          s.eta.sq = 1, s.v.sq = 1, s.x.sq = 1){
  
  #Function argument is correlation r.x.v, but code below is in terms of covariance
  s.x.v <- r.x.v/(sqrt(s.v.sq)*sqrt(s.x.sq))
  
  #Set up the covariance matrix for (x.i1, ..., x.iT, eta.i, v.i1, ... v.iT)
  I.T <- diag(N.t)
  one.T <- rep(1, N.t)
  zero.T <- rep(0, N.t)
  G <- diag(N.t - 1) #G sets up the predeterminedness of x
  G <- cbind(G, rep(0, N.t - 1))
  G <- rbind(rep(0, N.t), G)
  S.1 <- cbind(s.x.sq * I.T, s.x.eta * one.T, s.x.v * G)
  S.2 <- cbind(s.x.eta * t(one.T), s.eta.sq, t(zero.T))
  S.3 <- cbind(s.x.v * t(G), zero.T, s.v.sq * I.T)
  S <- rbind(S.1, S.2, S.3)
  colnames(S) <- NULL
  
  
  #Generate all of the simulations AT ONCE
  #Each row is an individual: a variate from the distribution of the vector
  #(x.i1, ..., x.iT, eta.i, v.i1, ... v.iT)
  sims <- mvrnorm(n = N.sims * N.i, mu =rep(0, ncol(S)), Sigma = S)
  
  
  #Extract the components corresponding to x, eta, and u
  x <- sims[,1:N.t]
  eta <- sims[,(N.t + 1)]
  v <- sims[,-(1:(N.t + 1))]
  
  #Clean up
  rm(sims, G, I.T, S, S.1, S.2, S.3, one.T, zero.T)
  
  
  
  #Initialize matrix, each of whose rows corresponds to an initial time period and each of whose rows corresponds to an individual
  y <- matrix(NA, nrow = N.i * N.sims, ncol = N.t)
  
  #Set y.0 to zero, the mean of its stationary distribution
  y.0 <- 0
  
  #Generate the first value of y using y.0
  y[,1] <- gamma * y.0 + theta * x[,1] + eta + v[,1]
  
  
  #Generate the remaining values of y
  for(j in 2:N.t){
    
    y[,j] <- gamma * y[,j-1] + theta * x[,j] + eta + v[,j]
    
  }#End for(i in 2:Time)
  
  #Clean up
  rm(eta, v)
  
  
  #Model matrices for model without a lag:
  #Construct the differenced LHS y vector for the model without a lag:
  #----------------------#
  #   y[i,2] - y[i,1]
  #   y[i,3] - y[i,2]
  #      .        .
  #      .        .
  #      .        .
  #   y[i,T] - y[i,T-1]
  #----------------------#
  y.tilde <- c(t(y[, -1])) - c(t(y[, -N.t]))

  
  #Construct the RHS variables for the model without a lag:
  #----------------------#
  #   x[i,2] - x[i,1]
  #   x[i,3] - x[i,2]
  #      .        .
  #      .        .
  #      .        .
  #   x[i,T] - x[i,T-1]
  #----------------------#
  X.tilde <- c(t(x[, -1])) - c(t(x[, -N.t]))

  
  #Construct the differenced LHS y vector for the model with a lag:
  #----------------------#
  #   y[i,3] - y[i,2]
  #   y[i,4] - y[i,3]
  #      .        .
  #      .        .
  #      .        .
  #   y[i,T] - y[i,T-1]
  #----------------------#
  y.tilde.L <- c(t(y[, -c(1, 2)])) - c(t(y[, -c(1, N.t)]))

  
  #Construct the RHS variables for the model with a lag:
  #---------------------------------------------------#
  #           X1                     X2
  #---------------------------------------------------#
  #   y[i,2]  -  y[i,1]        x[i,3] - x[i,2]
  #   y[i,3]  -  y[i,2]        x[i,4] - x[i,3]
  #      .         .              .        .
  #      .         .              .        .
  #      .         .              .        .
  #   y[i,T-1] - y[i,T-2]      x[i,T] - x[i,T-1]
  #---------------------------------------------------#
  X1.L <- c(t(y[, -c(1, N.t)])) - c(t(y[,-c(N.t - 1, N.t)]))
  X2.L <- c(t(x[,-c(1,2)])) - c(t(x[,-c(1, N.t)]))
  X.tilde.L <- cbind(X1.L, X2.L)
  rm(X1.L, X2.L)

  
  
  #Construct the Instrument Blocks for the model without a lag
  #---------------------------------------------------#
  #     Z.x.minus.plus       |         Z.x.plus
  #---------------------------------------------------#
  #   x[i,1]          0      |   x[i,2]        0
  #         .                |         .
  #          .               |          .
  #           .              |           .
  #   0        x[i,T-1]      |   0        x[i,T]
  #---------------------------------------------------#
  rows <- seq_len((N.t - 1) * N.i * N.sims)
  cols <- rep(1:(N.t - 1), times = N.i * N.sims)
  Z.x.minus.plus <- sparseMatrix(i = rows, j = cols, x = c(t(x[,-N.t])))
  Z.x.plus <- sparseMatrix(i = rows, j = cols, x = c(t(x[,-1])))

  
  
  #Construct the Instrument Blocks for the model with a lag
  #----------------------------------------------------#
  #      Z.y        |    Z.x.minus     |    Z.x
  #----------------------------------------------------#
  # y[i,1]       0  |  x[i,2]       0  |  x[i,3]     0
  #      .          |       .          |       .
  #       .         |        .         |        .
  #        .        |         .        |         .
  # 0     y[i,T-2]  |  0     x[i,T-1]  |  0     x[i,T]
  #----------------------------------------------------#
  rows <- seq_len((N.t - 2) * N.i * N.sims)
  cols <- rep(1:(N.t - 2), times = N.i * N.sims)
  Z.y <- sparseMatrix(i = rows, j = cols, x = c(t(y[,-c(N.t, N.t - 1)])))
  Z.x.minus <- sparseMatrix(i = rows, j = cols, x = c(t(x[,-c(1, N.t)])))
  Z.x <- sparseMatrix(i = rows, j = cols, x = c(t(x[,-c(1, 2)])))
  
  
  #------------------------------------------
  
  #Coerce to ordinary matrices
  Z.y <- as.matrix(Z.y)
  Z.x.minus <- as.matrix(Z.x.minus)
  Z.x <- as.matrix(Z.x)
  Z.x.minus.plus <- as.matrix(Z.x.minus.plus)
  Z.x.plus <- as.matrix(Z.x.plus)
  
  Z.LS <- cbind(Z.y, Z.x.minus, Z.x)
  Z.LW <- cbind(Z.y, Z.x.minus)
  Z.S <- cbind(Z.x.minus.plus, Z.x.plus)
  Z.W <- Z.x.minus.plus
  
  #Clean up
  rm(rows, cols, x, y, Z.x, Z.x.minus, Z.y, Z.x.plus, Z.x.minus.plus)
  
  
  #Now we'll break apart the instrument and model matrices for each replication
  #What is the size of each block? (A block is a replication)
  block.size.L <- length(y.tilde.L)/N.sims
  block.size <- length(y.tilde)/N.sims
  
  seq.block <- seq_len(block.size)
  seq.block.L <- seq_len(block.size.L)
  
  #Initialize dataframe to store estimates from each replication
  theta.hat <- matrix(NA, nrow = N.sims, ncol = 4)
  theta.hat <- data.frame(theta.hat)
  names(theta.hat) <- c('LW', 'LS', 'W', 'S')
  
  #Initialize dataframe to store GFIC and J-test from each replication
  GFIC <- J.test <- theta.hat
  
  
  #Loop through simulation replications
  for(k in 1:N.sims)
  {
    
    #Extract the block corresponding to the current replication for the unlagged model
    block.k <- (k - 1) * block.size + seq.block
    X.tilde.k <- X.tilde[block.k] #This is a vector, not a matrix!
    y.tilde.k <- y.tilde[block.k] #Also a vector, not a matrix!
    Z.W.k <- Z.W[block.k,]
    Z.S.k <- Z.S[block.k,]
    
    #Extract the block corresponding to the current replication for the lagged model
    block.L.k <- (k - 1) * block.size.L + seq.block.L
    X.tilde.L.k <- X.tilde.L[block.L.k,]
    y.tilde.L.k <- y.tilde.L[block.L.k] #A vector
    Z.LW.k <- Z.LW[block.L.k,]
    Z.LS.k <- Z.LS[block.L.k,]
    
    #Fit 2SLS with centered covariance matrices for everything except LW
    reg.LW <- fast2SLS(y.tilde.L.k, X.tilde.L.k, Z.LW.k, FALSE)
    reg.LS <- fast2SLS(y.tilde.L.k, X.tilde.L.k, Z.LS.k, TRUE)
    reg.W <- fast2SLS(y.tilde.k, as.matrix(X.tilde.k), as.matrix(Z.W.k), TRUE)
    reg.S <- fast2SLS(y.tilde.k, as.matrix(X.tilde.k), Z.S.k, TRUE)
    
    #Store the 2SLS estimates for each specification
    theta.LW <- reg.LW$b[2]
    theta.LS <- reg.LS$b[2]
    theta.W <- reg.W$b[1]
    theta.S <- reg.S$b[1]
    theta.hat[k,] <- cbind(theta.LW, theta.LS, theta.W, theta.S)
    
    #Store the J-test statistics for each specification, but use a centered version for LW in this instance (These are for the Andrews & Lu, 2001, procedure and they specify that the same method for estimating the variance matrix should be used for each specification.)
    J.LW <- fast2SLS(y.tilde.L.k, X.tilde.L.k, Z.LW.k, TRUE)$Jtest
    J.LS <- reg.LS$Jtest
    J.W <- reg.W$Jtest
    J.S <- reg.S$Jtest
    J.test[k,] <- cbind(J.LW, J.LS, J.W, J.S)
    
    
    #Estimates of bias parameters (tau, delta) and their covariance matrix
    delta <- sqrt(N.i) * reg.LW$b[1]
    Z.x.k <- Z.LS.k[,(2 * (N.t - 2) + 1):(3 * (N.t - 2))]
    tau.vec <- t(Z.x.k) %*% (y.tilde.L.k - X.tilde.L.k %*% reg.LW$b)/sqrt(N.i)
    tau <- mean(tau.vec)
    
    #-------------------------------------------------------------
    # Now we estimate the covariance matrix of (delta, tau). 
    #-------------------------------------------------------------
    
    #-------------------------------------------------------------
    # To begin we'll collect the quantities needed for the block
    # corresponding to tau.
    #-------------------------------------------------------------
    
    #-------------------------------------------------------------
    # We first need an estimate of the matrix I call Psi in my
    # derivations. Psi is formed from K.LW, which we've already
    # estimated above, and two other quantities:
    #
    #     xi.1 = E( x[i,t] * (y[i,t-1] - y[i,t-2]) )
    #     xi.2 = E( x[i,t] * (x[i,t] - x[i,t-1]) )
    #
    # Where the expectations are understood to be calculated in
    # the limit, i.e. with the mis-specification equal to zero.
    # We estimate these using their sample analogues, averaging 
    # over individuals and time by stationarity. For xi.1 we
    # can only use the time periods 3 through T whereas we can
    # use the periods 2 through T.
    #
    # An alternative way to proceed, without stationarity. In this
    # case we would calculate:
    #
    #     t(Z.x.k) %*% X.tilde.L.k/N.i 
    #
    # and us this quantity rather than filling in xi.1 and xi.2 
    # repeatedly, leading to the following estimator for Psi:
    #
    #       Psi <- -t(Z.x.k) %*% (X.tilde.L.k/N.i) %*% reg.LW$K
    #
    # Assuming that stationarity holds, however, this estimator is
    # less efficient.
    #-------------------------------------------------------------
    xi.1 <- mean(rowSums(Z.x.k) * X.tilde.L.k[,1])
    Z.x.plus.k <- Z.S.k[,((N.t - 1) + 1):(2 * (N.t - 1))]
    xi.2 <- mean(rowSums(Z.x.plus.k) * X.tilde.k)
    xi <- matrix(c(xi.1, xi.2), nrow = 1)
    i.Tminus2 <- matrix(1, nrow = (N.t - 2))
    Psi <- (-xi %x% i.Tminus2) %*% reg.LW$K #%x% is kronecker
    B <- (t(i.Tminus2)/(N.t - 2)) %*% cbind(Psi, diag(N.t - 2))
    
    #-------------------------------------------------------------
    # The quantitites needed for the block corresponding to delta
    # are easy: the first row of K.LW and some zeros.
    #
    #     A = c(K.LW[1,], zeros)
    #-------------------------------------------------------------
    A <- c(reg.LW$K[1,], rep(0, N.t - 2))
    
    #-------------------------------------------------------------
    # Finally we compute the variance matrix of our estimated bias
    # parameters as follows. The ordering is as follows:
    #
    #   [ delta^2         delta * tau ]
    #   [ tau * delta     tau^2       ]
    #-------------------------------------------------------------
    V.bias <- rbind(A, B) %*% reg.LS$Omega %*% t(rbind(A, B))
    
    #-------------------------------------------------------------
    # Using V.bias we construct the following asymptotically
    # unbiased estimators of the *squared* and *interacted* bias
    # parameters. (These are the quantities that appear in AMSE.)
    #-------------------------------------------------------------
    delta.squared <- delta^2 - V.bias[1,1]
    tau.squared <- tau^2 - V.bias[2,2]
    tau.delta <- tau * delta - V.bias[2,1]
    
    #-------------------------------------------------------------
    # Finally we're ready to calculate the GFIC, i.e. the AMSE 
    # estimate, for each specification.
    #-------------------------------------------------------------
    
    #-------------------------------------------------------------
    # GFIC for specification LW -- correct model, correct 
    # exogeneity assumption
    #-------------------------------------------------------------
    # By assumption the this specification is correct, hence there 
    # is no asymptotic bias. Thus, GFIC is simply an asymptotic
    # variance estimate. Again using the assumption of correct
    # specification, we use the uncentered, panel-robust variance
    # estimator. Since the variance matrix is arranged as follows:
    #
    #   [gamma^2         gamma * theta ]
    #   [theta * gamma   theta^2       ]
    #
    # We need the bottom right element.
    #-------------------------------------------------------------
    GFIC.LW <- reg.LW$V[2,2]
    
    #-------------------------------------------------------------
    # GFIC for specification LS -- correct model, incorrect 
    # exogeneity assumption
    #-------------------------------------------------------------
    # The estimtor of theta in this specification inherits an 
    # asymptotic bias from tau but not from delta. 
    #-------------------------------------------------------------
    tau.mat <- matrix(0, nrow = 3, ncol = 3)
    tau.mat[3,3] <- tau.squared
    inner <- tau.mat %x% (i.Tminus2 %*% t(i.Tminus2))
    GFIC.LS <- (reg.LS$K[2,] %*% inner) %*% t(reg.LS$K)[,2] + reg.LS$V[2,2]
    
    #-------------------------------------------------------------
    # GFIC for specification S -- incorrect model, incorrect 
    # exogeneity assumption
    #-------------------------------------------------------------
    # The estimtor of theta in this specification inherits an 
    # asymptotic bias from both delta and tau. The asymptotic bias
    # also involves two the following two quantities:
    #
    #     xi.3 = E( x[i,t] * (y[i,t] - y[i,t-1]) )
    #     xi.1 = E( x[i,t] * (y[i,t-1] - y[i,t-2]) )
    #
    # The second of these is simply xi.1 from above, but we need to
    # calculate the second, which we'll call xi.3. Again, the idea
    # is to use the assumption of stationarity so that we can
    # average over individuals and time to get a sample analogue.
    # These quantities emerge because the expression for the
    # scaled and centered estimator of theta under this particular
    # specification involves:
    #
    #   plim (Zw' %*% delta.y.lag.plus/n)
    #
    # Since we don't observe the zeroth time period, we can't
    # use the above directly, so we really do need stationarity.
    #-------------------------------------------------------------
    xi.3 <- mean(y.tilde.k * X.tilde.k)
    bias.mat <- matrix(0, nrow = 2, ncol = 2)
    bias.mat[1,1] <- delta.squared * xi.3^2
    bias.mat[2,2] <- delta.squared * xi.1^2 + 2 * tau.delta * xi.1 + tau.squared
    bias.mat[1,2] <- bias.mat[2,1] <- delta.squared * xi.1 * xi.3 + tau.delta * xi.3  
    i.Tminus1 <- matrix(1, nrow = (N.t - 1))
    inner <- bias.mat %x% (i.Tminus1 %*% t(i.Tminus1))
    GFIC.S <- reg.S$K %*% inner %*% t(reg.S$K) + reg.S$V #1x1 var matrix
    
    #-------------------------------------------------------------
    # GFIC for specification W -- incorrect model, exogeneity
    # assumption
    #-------------------------------------------------------------
    # The estimtor of theta in this specification inherits an 
    # asymptotic bias from both delta and tau. It also involves
    # the quantity xi.3 from above.
    #-------------------------------------------------------------
    inner <- xi.3^2 * delta.squared * i.Tminus1 %*% t(i.Tminus1)
    GFIC.W <- reg.W$K %*% inner %*% t(reg.W$K) + reg.W$V #1x1 var matrix
    
    #Store all of the GFIC values
    GFIC[k,] <- cbind(GFIC.LW, GFIC.LS, GFIC.W, GFIC.S)
  
    
    
  }#END for k
  
  

  
  
  #-------------------------------------------------------------
  # Selection by GFIC
  #-------------------------------------------------------------
  which.GFIC <- apply(GFIC, 1, which.min)
  theta.GFIC <- theta.hat[cbind(1:N.sims, which.GFIC)] 
  #The preceding command uses the fact that `[' can accept a matrix argument to extract a different column in each row of a dataframe.
  
  
  
  #-------------------------------------------------------------
  # Selection by Downward J-tests
  #
  # The idea is to pick a significance level and sequentially
  # test each specification (except LW) in order from most
  # to least restrictive: S, W, LS. A specification is
  # accepted if the test does not reject. If all three tests
  # reject, then we use specification LW.
  #
  # To carry out these tests, we first need the corresponding
  # critical values. The degrees of freedom for the chi-squared
  # distribution in this case equal the number of
  # overidentifying restrictions. Recall that:
  #   
  #   Spec   #Par   #Instr/t    N.t   #MCs      #OIR
  #   ----   ----   ---------   ---   ------    --------
  #    LW      2        2       T-2   2(T-2)    2(T-2)-2 = 2T-6
  #    LS      2        3       T-2   3(T-2)    3(T-2)-2 = 3T-8
  #    W       1        1       T-1   T-1       (T-1)-1  =  T-2
  #    S       1        2       T-1   2(T-1)    2(T-1)-1 = 2T-3
  #
  # Now, we want to check for rejection sequentially: S, W, LS
  # but this operation can in fact be vectorized. The idea is 
  # to create a condition that looks at whether each test 
  # rejects (TRUE) or fails to reject (FALSE) and returns a 
  # number from 1-4 corresponding to the column of theta.hat
  # that the sequence of tests indicates we should use as our
  # post-selection estimator. That is:
  #     
  #     1 = use LW
  #     2 = use LS
  #     3 = use W
  #     4 = use S
  #
  # If we let LS, W, and S denote, for the moment, indicators
  # for *rejection* by the corresponding J-test, what we need is
  # the following:
  #
  #     if(S == FALSE){return(4)}
  #     if((S == TRUE) && (W == FALSE)){return(3)}
  #     if((S == TRUE) && (W == TRUE) && (S = FALSE)){return(2)}
  #     else{return(1)}
  #
  # Using the equivalences TRUE = 1 and FALSE = 0, we can write
  # the above set of conditions succinctly as:
  #
  #     4 - S - S * W - S * W * LS
  #
  # Which is equal to:
  #
  #     4 - S * (1 + W * (1 + LS))
  #-------------------------------------------------------------
  df <- c(3 * N.t - 8, N.t - 2, 2 * N.t - 3)
  names(df) <- c("LS", "W", "S")
  c.10 <-qchisq(0.9, df)
  names(c.10) <- c("LS", "W", "S")
  c.5 <-qchisq(0.95, df)
  names(c.5) <- c("LS", "W", "S")
  reject.10 <- J.test[,-1] > matrix(c.10, byrow = TRUE, nrow = N.sims, ncol = 3)
  reject.5 <- J.test[,-1] > matrix(c.5, byrow = TRUE, nrow = N.sims, ncol = 3)
  reject.5 <- as.data.frame(reject.5)
  reject.10 <- as.data.frame(reject.10)
  select.J10 <- 4 - reject.10$S * (1 + reject.10$W * (1 + reject.10$LS)) 
  select.J5 <- 4 - reject.5$S * (1 + reject.5$W * (1 + reject.5$LS)) 
  
  theta.J10 <- theta.hat[cbind(1:N.sims, select.J10)]
  theta.J5 <- theta.hat[cbind(1:N.sims, select.J5)]
  
  #-------------------------------------------------------------
  # Selection by Andrews & Lu (2001) Criteria
  #
  # In the notation of this paper, |b| denotes the number of
  # parameters estimated in a given specification and |c| the
  # number of moment conditions used in estimation. Here I'll 
  # just call these b and c respectively. Thus, we have:
  #
  #   Spec   b = #Par   c = #MCs      c - b = #OIR
  #   ----   --------   ---------     ------------
  #    LW      2        2(T-2)        2(T-2)-2 = 2T-6
  #    LS      2        3(T-2)        3(T-2)-2 = 3T-8
  #    W       1        T-1           (T-1)-1  =  T-2
  #    S       1        2(T-1)        2(T-1)-1 = 2T-3
  #
  # The criteria from the paper are as follows, where J is the 
  # J-test statistic:
  #
  #   BIC-type:   J - (|c| - |b|) * log(n)
  #   AIC-type:   J - 2 * (|c| - |b|) 
  #   HQ-type:    J - Q * (|c| - |b|) * log(log(n))
  #
  # In the Hannan-Quinn-type criterion, Q should be > 2. I'll 
  # use Q = 2.01 as is often done in practice. The sample size
  # in this case, n, is the dimension that goes to infinity in
  # the limit. Hence, in our case it should be the number of
  # individuals, N.i.
  #-------------------------------------------------------------
  cminusb <- c(2 * N.t - 6, 3 * N.t - 8, N.t - 2, 2 * N.t - 3)
  names(cminusb) <- c("LW", "LS", "W", "S")
  AIC.penalty <- 2 * cminusb
  HQ.penalty <- 2.01 * cminusb * log(log(N.i)) 
  BIC.penalty <- cminusb * log(N.i)
  
  BIC.penalty <- matrix(BIC.penalty, byrow = TRUE, nrow = N.sims, ncol = length(BIC.penalty))
  BIC <- J.test - BIC.penalty
  select.BIC <- apply(BIC, 1, which.min)
  theta.BIC <- theta.hat[cbind(1:N.sims, select.BIC)]
  
  HQ.penalty <- matrix(HQ.penalty, byrow = TRUE, nrow = N.sims, ncol = length(HQ.penalty))
  HQ <- J.test - HQ.penalty
  select.HQ <- apply(HQ, 1, which.min)
  theta.HQ <- theta.hat[cbind(1:N.sims, select.HQ)]
  
  AIC.penalty <- matrix(AIC.penalty, byrow = TRUE, nrow = N.sims, ncol = length(AIC.penalty))
  AIC <- J.test - AIC.penalty
  select.AIC <- apply(AIC, 1, which.min)
  theta.AIC <- theta.hat[cbind(1:N.sims, select.AIC)]
  
  
  #-------------------------------------------------------------
  # Calculate RMSE of *everything at once*
  #
  # The ultimate goal is to compare the RMSE of each estimator
  # and the RMSE of the post-selection estimator arising from:
  # GFIC, Andrews-Lu procedures, and Downward J-test.
  #-------------------------------------------------------------
  everything <- cbind(theta.hat, theta.GFIC, theta.J10, theta.J5, 
                      theta.BIC, theta.HQ, theta.AIC)
  RMSE <- sqrt(colMeans(everything - theta)^2 
                          + apply(everything, 2, var))
  names(RMSE) <- c('LW', 'LS', 'W', 'S', 'GFIC', 'J10', 'J5', 'BIC', 'HQ', 'AIC')
  return(RMSE)
  
}#END RMSE.sim.FAST

#-----------------------------------------------------------------#
# END FUNCTION TO CARRY OUT THE SIMULATION A FIXED PARAMETERS
#-----------------------------------------------------------------#

#Test run to get a sense of how long things will take
#set.seed(1527)
#system.time(RMSE.sim.FAST(gamma = 0.85, r.x.v = 0.5, N.sims = 2000, N.t = 5, N.i = 500))
#  user  system elapsed 
#15.178   2.267  17.440 

#Look for bottlenecks
#Rprof()
#test <- RMSE.sim.FAST(gamma = 0.85, r.x.v = 0.5, N.sims = 2000, N.t = 5, N.i = 500)
#Rprof(NULL)
#summaryRprof()



#-----------------------------------------------------------------#
# FUNCTION TO CARRY OUT THE SIMULATION OVER A GRID
#-----------------------------------------------------------------#
RMSE.grid.FAST <- function(g, r, N.sims = 1000, 
                           N.i = 250, N.t = 3, theta = 0.5, s.x.eta = 0.2, 
                           s.eta.sq = 1, s.v.sq = 1, s.x.sq = 1)
{
  
  g.vector <- rep(g, each = length(r))
  r.vector <- rep(r, times = length(g))
  
  result <- mapply(RMSE.sim.FAST, gamma = g.vector, r.x.v = r.vector, 
                   N.sims = N.sims, N.i = N.i, N.t = N.t, theta = theta, 
                   s.x.eta = s.x.eta, s.eta.sq = s.eta.sq, s.v.sq = s.v.sq, 
                   s.x.sq = s.x.sq, SIMPLIFY = 'matrix')
  
  parameters <- cbind(g.vector, r.vector)
  colnames(parameters) <- c('gamma', 'r.x.v')
  cbind(parameters, t(result))
  
  
}#END RMSE.grid.FAST
#-----------------------------------------------------------------#
# END FUNCTION TO CARRY OUT THE SIMULATION OVER A GRID
#-----------------------------------------------------------------#




