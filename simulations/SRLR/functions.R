#-----------------------------------------------------------------#
#DESCRIPTION:   
# This code contains functions to carry out a simulation studying the properties GFIC model and moment selection in a linear, dynamic panel model. The correct model is:
#
#     y[i,t] = theta * x[i,t] + gamma * y[i,t-1] +  eta[i] + v[i,t]
#
#The random variables x[i,t], eta[i] and v[i,t] are normal with mean zero and constant variances: s.x.sq, s.eta.sq, s.v.sq. The covariance between x[i,t] and eta[i] is s.x.eta; the covariance between x[i,t] and v[i,t-1] is s.x.v while that between x[i,t] and v[i,k] for (k != t-1) is zero. That is, x is weakly exogenous (predetermined) but is correlated with the one-period lagged error.
#     The purpose of this simulation is to calculate the RMSE of four estimators of the paramters for varying values of gamma and s.x.v. The four estimators are as follows:
#
#     Name        Lagged y?     Assume x Strictly Exogenous?
#     --------    ---------     ----------------------------
#     b.LW        YES           NO (assume weakly exogenous)
#     b.W         NO            NO (assume weakly exogenous)
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
sourceCpp("fast2SLS.cpp")


#-----------------------------------------------------------------#
# FUNCTION TO CARRY OUT THE SIMULATION A FIXED PARAMETER VALUES
#-----------------------------------------------------------------#
RMSE.sim.FAST <- function(alpha1=0.4, alpha2, r.x.v, N.sims = 1000, 
                          N.i = 250, N.t = 5, beta = 0.5, s.x.eta = 0.2, 
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
  
  #Set y.0, y.lag1 to zero, the mean of its stationary distribution
  y.0 <- 0
  y.lag1 <-0
  
  
  #Generate the first value of y using y.0
  y[,1] <- alpha1 * y.0 + alpha2 * y.lag1 + beta * x[,1] + eta + v[,1]
  y[,2] <- alpha1 * y[,1] + alpha2 * y.0 + beta * x[,2] + eta + v[,2]
  
  #Generate the remaining values of y
  for(j in 3:N.t){
    
    y[,j] <- alpha1 * y[,j-1] + alpha2 * y[,j-2] + beta * x[,j] + eta + v[,j]
    
  }#End for(i in 3:Time)
  
  #Clean up
  rm(eta, v)
  
  #Model matrices for model without a lag:
  #Construct the differenced LHS y vector for the model without a lag:
  #----------------------#
  #   y[i,3] - y[i,2]
  #      .        .
  #      .        .
  #      .        .
  #   y[i,T] - y[i,T-1]
  #----------------------#
  y.tilde <- c(t(y[, -c(1,2)])) - c(t(y[, -c(1,N.t)]))
  
  
  #Construct the RHS variables for the model without a lag:
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
  X1 <- c(t(y[, -c(1, N.t)])) - c(t(y[,-c(N.t - 1, N.t)]))
  X2 <- c(t(x[,-c(1,2)])) - c(t(x[,-c(1, N.t)]))
  X.tilde <- cbind(X1, X2)
  rm(X1, X2)
  
  
  #Construct the differenced LHS y vector for the model with a lag:
  #----------------------#
  #   y[i,4] - y[i,3]
  #      .        .
  #      .        .
  #      .        .
  #   y[i,T] - y[i,T-1]
  #----------------------#
  y.tilde.L <- c(t(y[, -c(1, 2, 3)])) - c(t(y[, -c(1, 2, N.t)]))
  
  
  #Construct the RHS variables for the model with a lag:
  #------------------------------------------------------------------------#
  #           X1                     X2                      x3
  #------------------------------------------------------------------------#
  #   y[i,3]  -  y[i,2]        y[i,2] - y[i,1]        x[i,4] - x[i,3]
  #      .         .              .        .            .        .
  #      .         .              .        .            .        .
  #      .         .              .        .            .        .
  #   y[i,T-1] - y[i,T-2]      y[i,T-2] - y[i,T-3]    x[i,T] - x[i,T-1] 
  #------------------------------------------------------------------------#
  X1.L <- c(t(y[, -c(1, 2, N.t)])) - c(t(y[,-c(1, N.t - 1, N.t)]))
  X2.L <- c(t(y[, -c(1, N.t - 1, N.t)])) - c(t(y[,-c(N.t - 2, N.t - 1, N.t)]))
  X3.L <- c(t(x[,-c(1,2,3)])) - c(t(x[,-c(1,2, N.t)]))
  X.tilde.L <- cbind(X1.L, X2.L, X3.L)
  rm(X1.L, X3.L)
  
  
  
  #Construct the Instrument Blocks for the model without a lag

  #----------------------------------------------------#
  #      Z.y.minus2.plus     |      Z.x.minus.plus     
  #----------------------------------------------------#
  #   y[i,1]          0      |   x[i,2]        0
  #         .                |         .
  #          .               |          .
  #           .              |           .
  #   0        y[i,T-2]      |   0        x[i,T-1]
  #---------------------------------------------------#
  rows <- seq_len((N.t - 2) * N.i * N.sims)
  cols <- rep(1:(N.t - 2), times = N.i * N.sims)
  Z.y.minus2.plus <- sparseMatrix(i = rows, j = cols, x = c(t(x[,-c(N.t - 1, N.t)])))
  Z.x.minus.plus <- sparseMatrix(i = rows, j = cols, x = c(t(x[,-c(1, N.t)])))
  
  
  
  #Construct the Instrument Blocks for the model with a lag
  #---------------------------------------------------------------#
  #      Z.y.minus2        |    Z.y.minus3     |    Z.x.minus
  #---------------------------------------------------------------#
  #    y[i,2]       0      |  y[i,1]       0   |  x[i,3]     0
  #        .               |       .           |       .
  #         .              |        .          |        .
  #          .             |         .         |         .
  #    0       y[i,T-2]    |  0       y[i,T-3] |  0       x[i,T-1]
  #---------------------------------------------------------------#
  rows <- seq_len((N.t - 3) * N.i * N.sims)
  cols <- rep(1:(N.t - 3), times = N.i * N.sims)
  Z.y.minus2 <- sparseMatrix(i = rows, j = cols, x = c(t(y[,-c(N.t, N.t - 1, 1)])))
  Z.y.minus3 <- sparseMatrix(i = rows, j = cols, x = c(t(y[,-c(N.t, N.t - 1, N.t - 2)])))
  Z.x.minus <- sparseMatrix(i = rows, j = cols, x = c(t(x[,-c(1, 2, N.t)])))
  
  
  #------------------------------------------
  
  #Coerce to ordinary matrices
  Z.y.minus2 <- as.matrix(Z.y.minus2)
  Z.y.minus2.plus <- as.matrix(Z.y.minus2.plus)
  Z.y.minus3 <- as.matrix(Z.y.minus3)
  Z.x.minus <- as.matrix(Z.x.minus)
  Z.x.minus.plus <- as.matrix(Z.x.minus.plus)

  
  Z.LW <- cbind(Z.y.minus2, Z.y.minus3, Z.x.minus)
  Z.W <- cbind(Z.y.minus2.plus, Z.x.minus.plus)
  
  #Clean up
  rm(rows, cols, x, y, Z.y.minus2.plus, Z.y.minus3, Z.x.minus.plus)
  
  
  
  #Now we'll break apart the instrument and model matrices for each replication
  #What is the size of each block? (A block is a replication)
  block.size.L <- length(y.tilde.L)/N.sims
  block.size <- length(y.tilde)/N.sims
  
  seq.block <- seq_len(block.size)
  seq.block.L <- seq_len(block.size.L)
  
  #Initialize dataframe to store estimates from each replication
  beta.hat <- matrix(NA, nrow = N.sims, ncol = 2)
  beta.hat <- data.frame(beta.hat)
  names(beta.hat) <- c('LW', 'W')
  
  longrun.hat <- matrix(NA, nrow = N.sims, ncol = 2)
  longrun.hat <- data.frame(longrun.hat)
  names(longrun.hat) <- c('LW', 'W')
  
  #Initialize dataframe to store GFIC and J-test from each replication
  GFIC <- GFIC.longrun <- J.test <- beta.hat
  
  
  #Loop through simulation replications
  for(k in 1:N.sims)
  {
    
    #Extract the block corresponding to the current replication for the unlagged model
    block.k <- (k - 1) * block.size + seq.block
    X.tilde.k <- X.tilde[block.k,] #This is a vector, not a matrix!
    y.tilde.k <- y.tilde[block.k] #Also a vector, not a matrix!
    Z.W.k <- Z.W[block.k,]
    
    #Extract the block corresponding to the current replication for the lagged model
    block.L.k <- (k - 1) * block.size.L + seq.block.L
    X.tilde.L.k <- X.tilde.L[block.L.k,]
    y.tilde.L.k <- y.tilde.L[block.L.k] #A vector
    Z.LW.k <- Z.LW[block.L.k,]
    
    #Fit 2SLS with centered covariance matrices for everything except LW
    reg.LW <- fast2SLS(y.tilde.L.k, X.tilde.L.k, Z.LW.k, center = TRUE, n = N.i)
    reg.W <- fast2SLS(y.tilde.k, as.matrix(X.tilde.k), as.matrix(Z.W.k), center = TRUE, n = N.i)
    
    #Store the 2SLS estimates for each specification
    alpha1.LW <- reg.LW$b[1]
    alpha2.LW <- reg.LW$b[2]
    beta.LW <- reg.LW$b[3]
    
    alpha1.W <-reg.W$b[1]
    beta.W <- reg.W$b[2]
    
    longrun.LW <- beta.LW/(1-alpha1.LW-alpha2.LW)
    longrun.W <- beta.W/(1-alpha1.W)

    beta.hat[k,] <- cbind(beta.LW, beta.W)
    longrun.hat[k,] <-cbind(longrun.LW, longrun.W)
    
    hprime.LW <- c(beta.LW/(1-alpha1.LW)^2, beta.LW/(1-alpha1.LW)^2, 1/(1-alpha1.LW))
    hprime.W <- c(beta.W/(1-alpha1.W)^2, 1/(1-alpha1.W))
    
    
    #Store the J-test statistics for each specification, but use a centered version for LW in this instance (These are for the Andrews & Lu, 2001, procedure and they specify that the same method for estimating the variance matrix should be used for each specification.)
    J.LW <- fast2SLS(y.tilde.L.k, X.tilde.L.k, Z.LW.k, center = TRUE, n = N.i)$Jtest
    J.W <- reg.W$Jtest
    J.test[k,] <- cbind(J.LW, J.W)
    
    
    #Estimates of bias parameter delta
    delta <- sqrt(N.i) * reg.LW$b[2]

  
    #-------------------------------------------------------------
    # We first need an estimate of the matrix I call Psi in my
    # derivations. Psi is formed from K.LW, which we've already
    # estimated above, and two other quantities:
    
    #    xi.0 = E( y[i,t-2]*(y[i,t-2]-y[i,t-3]) )
    #    xi.1 = E( x[i,t-1]*(y[i,t-2]-y[i,t-3]) )
    
  
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
  
    xi.0 <- mean(rowSums(Z.y.minus2)*X2.L)
    xi.1 <- mean(rowSums(Z.x.minus)*X2.L)
    
    
    xi <- matrix(c(xi.0, xi.1), nrow = 1)
    i.Tminus2 <- matrix(1, nrow = (N.t - 2))
    
    
    
    #-------------------------------------------------------------
    # Finally we compute the variance matrix of our estimated bias
    # parameters delta^2 as follows. 
    #-------------------------------------------------------------
    v.bias <- cbind(0,1,0)%*%reg.LW$K%*%reg.LW$Omega%*%t(reg.LW$K)%*%rbind(0,1,0)
    
    #-------------------------------------------------------------
    # We construct the following asymptotically
    # unbiased estimators of the *squared* bias
    # parameters. (These are the quantities that appear in AMSE.)
    #-------------------------------------------------------------
    delta.squared <- delta^2 - v.bias
    #If the squared bias estimate is negative, set it to zero
    
    if(delta.squared >= 0){
      delta.squared = delta.squared;
    } else {
      delta.squared = 0;
    }

    #-------------------------------------------------------------
    # Finally we're ready to calculate the GFIC, i.e. the AMSE 
    # estimate, for each specification.
    #-------------------------------------------------------------
    
    #-------------------------------------------------------------
    # GFIC for specification LW -- correct model
    #-------------------------------------------------------------
    GFIC.LW <- reg.LW$V[3,3]
    
    GFIC.LWlongrun <- hprime.LW%*%reg.LW$V%*%matrix(hprime.LW)
    
    
    #-------------------------------------------------------------
    # GFIC for specification W -- incorrect model
    #-------------------------------------------------------------
    
    bias.mat <- matrix(0, nrow = 2, ncol = 2)
    bias.mat[1,1] <- delta.squared * xi.0^2
    bias.mat[2,2] <- delta.squared * xi.1^2 
    bias.mat[1,2] <- bias.mat[2,1] <- delta.squared * xi.0 * xi.1 
    inner <- bias.mat %x% (i.Tminus2 %*% t(i.Tminus2))
  
      GFIC.W <- reg.W$K[2,] %*% inner %*% matrix(reg.W$K[2,])  + reg.W$V[2,2] #1x1 var matrix
      
      GFIC.Wlongrun <- hprime.W%*%reg.W$K%*%inner%*%matrix(hprime.W%*%reg.W$K) + hprime.W%*%reg.W$V%*%matrix(hprime.W)
    
    #Store all of the GFIC values
    GFIC[k,] <- cbind(GFIC.LW, GFIC.W)
    GFIC.longrun[k,] <-cbind(GFIC.LWlongrun, GFIC.Wlongrun)
    
    
    
  # Notation: V in the paper  = Omega in the code  
  #           K*V*K in the paper = V in the code  
    
    
  }#END for k
  
  
  
  
  
  #-------------------------------------------------------------
  # Selection by GFIC
  #-------------------------------------------------------------
  which.GFIC <- apply(GFIC, 1, which.min)
  beta.GFIC <- beta.hat[cbind(1:N.sims, which.GFIC)] 
  
  which.GFIClongrun <- apply(GFIC.longrun, 1, which.min)
  longrun.GFIC <- longrun.hat[cbind(1:N.sims, which.GFIClongrun)] 
  
  #summary(which.GFIC)
  tabulate(which.GFIC)
  tabulate(which.GFIClongrun)
  
  
  #The preceding command uses the fact that `[' can accept a matrix argument to extract a different column in each row of a dataframe.
  
  
  
 
  
  
  #-------------------------------------------------------------
  # Calculate RMSE of *everything at once*
  #
  # The ultimate goal is to compare the RMSE of each estimator
  # and the RMSE of the post-selection estimator arising from:
  # GFIC, Andrews-Lu procedures, and Downward J-test.
  #-------------------------------------------------------------
  everything <- cbind(beta.hat, beta.GFIC)
  
  ## consider MAD (median absolute error)
  abs_dev  = abs(everything - beta)
  RMSE <- apply(abs_dev, 2, median)
  
  
  names(RMSE) <- c('LW', 'W', 'GFIC')
  
  everything.longrun <-cbind(longrun.hat, longrun.GFIC)
 
 
 ## consider MAD (median absolute error) 
 abs_dev.long  = abs(everything.longrun - beta/(1-alpha1-alpha2))
 RMSE.longrun <- apply(abs_dev.long, 2, median)
 


  names(RMSE.longrun) <- c('LW_LR', 'W_LR', 'GFIC_LR')
  
  final = c(RMSE, RMSE.longrun)
  
  return(final)
  
}#END RMSE.sim.FAST


#-----------------------------------------------------------------#
# FUNCTION TO CARRY OUT THE SIMULATION OVER A GRID
#-----------------------------------------------------------------#
RMSE.grid.FAST <- function(alpha1 = 0.4, g, r, N.sims=1000, 
                           N.i = 250, N.t = 5, beta = 0.5, s.x.eta = 0.2, 
                           s.eta.sq = 1, s.v.sq = 1, s.x.sq = 1)
{
  
  g.vector <- rep(g, each = length(r))
  r.vector <- rep(r, times = length(g))
  
  result <- mapply(RMSE.sim.FAST, alpha1 = alpha1, alpha2 = g.vector, r.x.v = r.vector, 
                   N.sims = N.sims, N.i = N.i, N.t = N.t, beta = beta, 
                   s.x.eta = s.x.eta, s.eta.sq = s.eta.sq, s.v.sq = s.v.sq, 
                   s.x.sq = s.x.sq, SIMPLIFY = 'matrix')
  
  parameters <- cbind(g.vector, r.vector)
  colnames(parameters) <- c('alpha2', 'r.x.v')
  cbind(parameters, t(result))
  
  
}#END RMSE.grid.FAST
#-----------------------------------------------------------------#
# END FUNCTION TO CARRY OUT THE SIMULATION OVER A GRID
#-----------------------------------------------------------------#

