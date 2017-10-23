#-----------------------------------------------------------------#
#SCRIPT:  rmse_simulation_FAST.R
#AUTHOR:  Frank DiTraglia
#EMAIL:   fditraglia@gmail.com 
#DATE:    June 6th, 2012
#-----------------------------------------------------------------#
#DESCRIPTION:   
#
#     This script contains functions to carry out a simulation studying the effect of mis-specification on the MSE of a target parameter in a linear, dynamic panel model. The correct model is:
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
#     To answer this question, this script compares the RMSE of the 2SLS estimator for b.LW, b.LS, b.W and b.S. The simulation setup described above is similar to that of Andrews and Lu (2001) with two main differences. First, I do not include a constant term. Second, I do not create extra pre-sample observations to accomodate the lagged y. That is, in my simulation omitting the lag gives an extra time period. This is important because in a short panel additional time periods could considerably reduce the variance of an estimator. Third, I set the initial value of y, y.0, equal to the mean of its stationary distribution (zero) as in Arellano and Bond (1991), unlike Andrews and Lu. Fourth, I use only strong instruments rather than all instruments available under the assumptions, essentially the standard Arellano-Bond Moment Conditions.
#-----------------------------------------------------------------#




#Load packages used below
library(MASS)
library(Matrix)





#-----------------------------------------------------------------#
# FUNCTION TO CARRY OUT 2SLS ESTIMATION
#-----------------------------------------------------------------#
reg2SLS.FAST <- function(X, y, Z){
  
  g <- solve(crossprod(Z), crossprod(Z,X))
  X.hat <- Z %*% g
  solve(crossprod(X.hat), crossprod(X.hat, y))
  
}#END reg2SLS.FAST
#-----------------------------------------------------------------#
# END FUNCTION TO CARRY OUT 2SLS ESTIMATION
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
# 

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
  y.panel.2.T <- y[, -1]
  y.2.T <- c(t(y.panel.2.T))
  y.panel.1.Tminus1 <- y[, -N.t]
  y.1.Tminus1 <- c(t(y.panel.1.Tminus1))
  y.tilde <- y.2.T - y.1.Tminus1
  
  #Clean up
  rm(y.panel.2.T, y.2.T, y.panel.1.Tminus1, y.1.Tminus1)

  
  
  
  #Construct the RHS variables for the model without a lag:
  #----------------------#
  #   x[i,2] - x[i,1]
  #   x[i,3] - x[i,2]
  #      .        .
  #      .        .
  #      .        .
  #   x[i,T] - x[i,T-1]
  #----------------------#
  x.panel.2.T <- x[, -1]
  x.2.T <- c(t(x.panel.2.T))
  x.panel.1.Tminus1 <- x[, -N.t]
  x.1.Tminus1 <- c(t(x.panel.1.Tminus1))
  X.tilde <- x.2.T - x.1.Tminus1
  
  #Clean up
  rm(x.panel.2.T, x.2.T, x.panel.1.Tminus1, x.1.Tminus1)

  
  
  
  #Construct the differenced LHS y vector for the model with a lag:
  #----------------------#
  #   y[i,3] - y[i,2]
  #   y[i,4] - y[i,3]
  #      .        .
  #      .        .
  #      .        .
  #   y[i,T] - y[i,T-1]
  #----------------------#
  y.panel.3.T <- y[, -c(1, 2)]
  y.3.T <- c(t(y.panel.3.T))
  y.panel.2.Tminus1 <- y[, -c(1, N.t)]
  y.2.Tminus1 <- c(t(y.panel.2.Tminus1))
  y.tilde.L <- y.3.T - y.2.Tminus1
  
  #Clean up
  rm(y.panel.3.T, y.3.T, y.panel.2.Tminus1) #Don't remove y.2.Tminus1 yet! We need it below!

  
  
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
  y.panel.1.Tminus2 <- y[,-c(N.t - 1, N.t)]
  y.1.Tminus2 <- c(t(y.panel.1.Tminus2))
  X1.L <- y.2.Tminus1 - y.1.Tminus2
  x.panel.3.T <- x[,-c(1,2)]
  x.3.T <- c(t(x.panel.3.T))
  x.panel.2.Tminus1 <- x[,-c(1, N.t)]
  x.2.Tminus1 <- c(t(x.panel.2.Tminus1))
  X2.L <- x.3.T - x.2.Tminus1 
  X.tilde.L <- cbind(X1.L, X2.L)
  
  #Clean up
  rm(y.2.Tminus1, y.panel.1.Tminus2, y.1.Tminus2, X1.L, x.panel.3.T, x.3.T, x.panel.2.Tminus1, x.2.Tminus1, X2.L)
  
  
  
  #Instrument matrices for S: no lag, strictly exogenous
  x1.S <- c(t(x[,-N.t]))
  x2.S <- c(t(x[,-1]))
  values.S <- c(x1.S, x2.S)
  rm(x1.S, x2.S)#Clean up
  
  rows.S <- rep(seq_len(N.i * N.sims * (N.t - 1)), times = 2)
  index.S <- seq_len(N.t - 1)
  col1.S <- rep(index.S * 2 - 1, times = N.i * N.sims)
  col2.S <- rep(index.S * 2, times = N.i * N.sims)
  cols.S <- c(col1.S, col2.S)
  rm(index.S, col1.S, col2.S)
  
  Z.S <- sparseMatrix(i = rows.S, j = cols.S, x = values.S)
    #spMatrix(nrow = length(rows.S), ncol = 2*(N.t - 1), i = rows.S, j = cols.S, x = values.S)#sparseMatrix() causes memory issues
  
  #Clean up
  rm(values.S, rows.S, cols.S)
  
  #Instrument matrices for W: no lag, weakly exogenous
  values.W <- c(t(x[,-N.t]))
  rows.W <- seq_len(N.i * N.sims * (N.t - 1))
  cols.W <- rep(seq_len(N.t - 1), times = N.i * N.sims)
  Z.W <- sparseMatrix(i = rows.W, j = cols.W, x = values.W)
    
  #spMatrix(nrow = length(rows.W), ncol = (N.t - 1), i = rows.W, j = cols.W, x = values.W)#sparseMatrix() causes memory issues
  
  #Clean up
  rm(values.W, rows.W, cols.W)
  
  
  #Instrument matrices for LW: lag, weakly exogenous
  x.LW <- c(t(x[,-c(1, N.t)]))
  y.LW <- c(t(y[,-c(N.t, N.t - 1)]))
  values.LW <- c(y.LW, x.LW)
  rm(x.LW, y.LW)
  
  rows.LW <- rep(seq_len(N.i * N.sims * (N.t - 2)), times = 2)
  index.LW <- seq_len(N.t - 2)
  col1.LW <- rep(index.LW * 2 - 1, times = N.i * N.sims)
  col2.LW <- rep(index.LW * 2, times = N.i * N.sims)
  cols.LW <- c(col1.LW, col2.LW)
  rm(index.LW, col1.LW, col2.LW)
  
  Z.LW <- sparseMatrix(i = rows.LW, j = cols.LW, x = values.LW) 
  #spMatrix(nrow = length(rows.LW), ncol = 2*(N.t - 2), i = rows.LW, j = cols.LW, x = values.LW)#sparseMatrix() causes memory issues
  
  #Clean up
  rm(values.LW, rows.LW, cols.LW)
  
  
  #Instrument matrices for LS: lag, strictly exogenous
  x1.LS <- c(t(x[,-c(1, N.t)]))
  x2.LS <- c(t(x[,-c(1, 2)]))
  y.LS <- c(t(y[,-c(N.t, N.t - 1)]))
  values.LS <- c(y.LS, x1.LS, x2.LS)
  rm(x1.LS, x2.LS, y.LS)
  
  rows.LS <- rep(seq_len(N.i * N.sims * (N.t - 2)), times = 3)
  index.LS <- seq_len(N.t - 2)
  col1.LS <- rep(index.LS * 3 - 2, times = N.i * N.sims)
  col2.LS <- rep(index.LS * 3 - 1, times = N.i * N.sims)
  col3.LS <- rep(index.LS * 3, times = N.i * N.sims)
  cols.LS <- c(col1.LS, col2.LS, col3.LS)
  rm(index.LS, col1.LS, col2.LS, col3.LS)

  
  Z.LS <- sparseMatrix(i = rows.LS, j = cols.LS, x = values.LS)
    #spMatrix(nrow = length(rows.LS), ncol = 3*(N.t - 2), i = rows.LS, j = cols.LS, x = values.LS)#sparseMatrix() causes memory issues
  
  #Now that we've set up the model and instrument matrices, we no longer require the raw data
  rm(values.LS, rows.LS, cols.LS, x, y)

  
  #Coerce to ordinary matrices
  Z.LS <- as.matrix(Z.LS)
  Z.LW <- as.matrix(Z.LW)
  Z.S <- as.matrix(Z.S)
  Z.W <- as.matrix(Z.W)
  
  
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
    
    theta.LW <- reg2SLS.FAST(X = X.tilde.L.k, y = y.tilde.L.k, Z = Z.LW.k)[2]
    theta.LS <- reg2SLS.FAST(X = X.tilde.L.k, y = y.tilde.L.k, Z = Z.LS.k)[2]
    theta.W <- reg2SLS.FAST(X = X.tilde.k, y = y.tilde.k, Z = Z.W.k)[1]
    theta.S <- reg2SLS.FAST(X = X.tilde.k, y = y.tilde.k, Z = Z.S.k)[1]
    
    theta.hat[k,] <- cbind(theta.LW, theta.LS, theta.W, theta.S)
    
  }#END for k
  
  
  #Calculate RMSE of each estimator
  bias <- apply(theta.hat - theta, 2, mean)
  variance <- apply(theta.hat, 2, var)
  sqrt(bias^2 + variance)
  
  
}#END RMSE.sim.FAST

#-----------------------------------------------------------------#
# END FUNCTION TO CARRY OUT THE SIMULATION A FIXED PARAMETERS
#-----------------------------------------------------------------#

#set.seed(1527)
#memory.limit(size = 3000)
#system.time(RMSE.sim.FAST(gamma = 0.85, r.x.v = 0.5, N.sims = 2000, N.t = 5, N.i = 500))
#user  system elapsed 
#8.722   2.825  11.607 

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





