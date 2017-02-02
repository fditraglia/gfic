#Load package for simulating from MV normal
library(mvtnorm)

#---------------------------------------------------------------------#
#                       FUNCTION Dpanel.sim()
#---------------------------------------------------------------------#
#This function generates one simulation replication from a simple dynamic panel model similar to that from Andrews and Lu (2001). Unlike Andrews and Lu, this model does not include a constant and only generates exactly as many observations (per i) as there are time periods: there are no ``pre-sample'' observations except for y.0 which is used to initialize the simulation. 
#---------------------------------------------------------------------#
#   a             Effect of lagged dependent variable y.[i,t-1]
#   b             Effect of regressor x
#   s.x.eta       Covariance between regressor x and individual effect eta.
#   s.x.v         Covariance between regressor x and lagged error term v
#   s.eta.sq      Variance of individual effect eta
#   s.v.sq        Variance of error term
#   s.x.sq        Variance of regressor x
#   N             Number of individuals
#   Time          Number of time periods
#---------------------------------------------------------------------#

Dpanel.sim <- function(N = 250, Time = 3, a = 0.85, b = 0.5, 
                       s.x.eta = 0.2, s.x.v = 0.5, s.eta.sq = 1, 
                       s.v.sq = 1, s.x.sq = 1){
  
  
  
  #Set up identity matrix and columns of ones and zeros used in covariance matrix
  I.T <- diag(Time)
  one.T <- rep(1, Time)
  zero.T <- rep(0, Time)
  
  #This matrix sets up the predeterminedness of x
  Gamma <- diag(Time - 1)
  Gamma <- cbind(Gamma, rep(0, Time - 1))
  Gamma <- rbind(rep(0, Time), Gamma)
  
  #Set up the covariance matrix for (x.i1, ..., x.iT, eta.i, v.i1, ... v.iT)
  S.1 <- cbind(s.x.sq * I.T, s.x.eta * one.T, s.x.v * Gamma)
  S.2 <- cbind(s.x.eta * t(one.T), s.eta.sq, t(zero.T))
  S.3 <- cbind(s.x.v * t(Gamma), zero.T, s.v.sq * I.T)
  S <- rbind(S.1, S.2, S.3)
  colnames(S) <- NULL
  
  #Generate the simulations: each row is an individual, i.e. a variate from the distribution of the vector:
  #(x.i1, ..., x.iT, eta.i, v.i1, ... v.iT)
  sims <- rmvnorm(N, sigma = S) 
  #Extract the components corresponding to x, eta, and u
  x <- sims[,1:Time]
  eta <- sims[,(Time + 1)]
  v<- sims[,-(1:(Time + 1))]
  
  
  #Initialize matrix, each of whose rows corresponds to an initial time period and each of whose rows corresponds to an individual
  y <- matrix(NA, nrow = N, ncol = Time)
  
    
  #Set y.0 to zero, the mean of its stationary distribution
  y.0 <- 0
  
  
  #Generate the first value of y using y.0
  y[,1] <- a * y.0 + b * x[,1] + eta + v[,1]
  
  
  #Generate the remaining values of y
  for(j in 2:Time){
    
    y[,j] <- a * y[,j-1] + b * x[,j] + eta + v[,j]
    
  }#End for(i in 2:Time)
  

  #Return the observables x, y
  time.indices <- paste('t=', 1:Time, sep = '')
  unit.indices <- paste('i=', 1:N, sep = '')
  colnames(x) <- time.indices
  rownames(x) <- unit.indices
  colnames(y) <- time.indices
  rownames(y) <- unit.indices
  
  out <- list(x = x, y)
  names(out) <- c('x', 'y')
  
  return(out)
  
  
}#END Dpanel.sim