#This setup follows the "base design" of Arellano & Bond 1991 Table 1.
ABsim <- function(a){
  
  dgp_cpp(a1 = a, a2 = 0, g = 0.1, N_i = 100, N_t = 7, 
          burn_in = 10, b = 1, r = 0.9, theta = 0, 
          s_e = sqrt(0.9), s_eta = 1, s_v = 1)
}


#Simple function for the one-step Arellano Bond Estimator as employed in the first column of Table 1 of their paper. Corresponds to AMsim above.
ABfit <- function(x, y){
  
  N.i <- nrow(x)
  N.t <- ncol(x)
  
  #Shortcut to avoid running seq_len(N.i) over and over
  individuals <- seq_len(N.i)
  
  #Label rows and columns for sanity-checking
  rownames(y) <- rownames(x) <- paste('i=', 1:N.i, sep = '')
  colnames(y) <- colnames(x) <- paste('t=', 1:N.t, sep = '')
  
  #Take first differences but "pad" with NAs
  x.diff <- t(apply(x, 1, function(v) diff(c(NA, v))))
  y.diff <- t(apply(y, 1, function(v) diff(c(NA, v))))
  
  #Convert to lists in which each element is an individual
  x.diff <- lapply(individuals, function(i) x.diff[i,])
  y.diff <- lapply(individuals, function(i) y.diff[i,])
  y <- lapply(individuals, function(i) y[i,])
  
  #Instruments for each individual (list of lists)
  Z.i <- function(i){
    
    lapply(3:N.t, function(j) c(y[[i]][1:(j - 2)], x.diff[[i]][j]))
    
  }
  
  Z <- lapply(individuals, Z.i)
  
  #Convert to list of sparse matrices
  Z <- lapply(individuals, function(i) t(bdiag(Z[[i]])))
  
  
  #Regressors for each individual (list of matrices)
  X.tilde.i <- function(i){
    
    cbind(y.diff[[i]][2:(N.t - 1)], x.diff[[i]][3:N.t])
    
  }
  
  X.tilde <- lapply(individuals, X.tilde.i)
  
  #Outcomes for each individual (list of vectors)
  y.tilde.i <- function(i){
    
    y.diff[[i]][3:N.t]
    
  }
  y.tilde <- lapply(individuals, y.tilde.i)
  
  
  XZ <- lapply(individuals, function(i) crossprod(X.tilde[[i]], Z[[i]]))
  XZ <- Reduce('+', XZ) #Sum over all individuals
  
  Zy <- lapply(individuals, function(i) crossprod(Z[[i]], y.tilde[[i]]))
  Zy <- Reduce('+', Zy)
  
  
  #Probably better not to use "solve" here but I'm not sure how to handle a qr decomposition for sparse matrices...
  H <- bandSparse(5, 5, k = c(0,-1, 1), list(rep(2, 5), rep(-1, 5), rep(-1, 5)))
  W.inv <- lapply(individuals, function(i) t(Z[[i]]) %*% H %*% Z[[i]])
  W.inv <- Reduce('+', W.inv) / N.i
  W <- solve(W.inv)
  
  XZW <- XZ %*% W
  K.inv <- XZW %*% t(XZ)
  K <- chol2inv(qr.R(qr(K.inv)))
  b <- solve(K.inv) %*% XZW %*% Zy
  #Direct calculation is: solve(XZ %*% W %*% t(XZ)) %*% XZ %*% W %*% Zy  
  row.names(b) <- c('alpha', 'beta')
  return(b)
  
  
}#END ABfit



#R version of dgp.cpp - much slower but gives same results.
dgpR <- function(a1, a2, g, N.i, N.t, burn.in = 10, 
                 b = 1, r = 0, theta = 0.2, s.e = 1, 
                 s.eta = 1, s.v = 1){
  
  #Generate Exogenous Errors
  e <- rnorm(N.i * (N.t + burn.in), mean = 0, sd = s.e)
  e <- matrix(e, nrow = N.i)
  v <- rnorm(N.i * (N.t + burn.in), mean = 0, sd = s.v)
  v <- matrix(v, nrow = N.i)
  eta <- rnorm(N.i, mean = 0, sd = s.eta)
  
  #Matrices to Store Endogenous Variables
  y <- x <- xi <- matrix(NA, nrow = N.i, ncol = N.t + burn.in)
  
  #First Initialization Step - Set Presample Observations to 0
  xi[, 1] <- e[, 1]
  x[, 1] <- theta * eta + xi[, 1]
  y[, 1] <- b * x[, 1] + eta + v[, 1]

  #Second Initialization Step - Set Presample Observations to 0
  xi[, 2] <- r * xi[, 1] + e[, 2]
  x[, 2] <- theta * eta + g * v[, 1] + xi[, 2]
  y[, 2] <- a1 * y[, 1] + b * x[, 2] + eta + v[, 2]
  
  #Generate Endogenous Variables
  for(i in 3:(N.t + burn.in)){
    
    xi[, i] <- r * xi[, i - 1] + e[, i]
    x[, i] <- theta * eta + g * v[, i - 1] + xi[, i]
    y[, i] <- a1 * y[, i - 1] + a2 * y[, i - 2] + b * x[, i] + eta + v[, i]
    
  }#END for

  #Discard Burn-in Samples
  x <- x[, -c(1:burn.in)]
  y <- y[, -c(1:burn.in)]
  
  out <- list(x = x, y = y)
  return(out)
  
}#END dgpR


