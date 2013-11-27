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


#foo <- dgpR(a1 = 0.3, a2 = -0.1, g = 0.1, N.i = 25, N.t = 4,r = 0.1)
#foo$x
#foo$y


