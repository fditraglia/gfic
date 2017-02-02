#-----------------------------------------------------------------#
#FILENAME:  model_matrices.R
#AUTHOR:    Frank DiTraglia
#EMAIL:     fditraglia@gmail.com 
#DATE:      June 1st, 2012
#-----------------------------------------------------------------#
#DESCRIPTION:   
#
#       This script contains two functions to set up ``y'' vector and ``X'' matrix for the panel model:
#
#     y[i,t] = gamma * y[i,t-1] + theta * x[i,t] + eta[i] + v[i,t]
#
#where x is a scalar. Estimation is in first differences to remove the individual effect eta[i]. The function model.matrices excludes the lag (hence using an additional time period), while model.matrices.L includes the lag. Both functions take as their input two matrices, x.panel and y.panel, where columns are time periods and rows individuals.
#-----------------------------------------------------------------#


model.matrices <- function(x.panel, y.panel){
  
  N.t <- ncol(x.panel)
  
  #Construct the differenced LHS y vector for the model without a lag:
  #----------------------#
  #   y[i,2] - y[i,1]
  #   y[i,3] - y[i,2]
  #      .        .
  #      .        .
  #      .        .
  #   y[i,T] - y[i,T-1]
  #----------------------#
  y.panel.2.T <- y.panel[, -1]
  y.2.T <- c(t(y.panel.2.T))
  y.panel.1.Tminus1 <- y.panel[, -N.t]
  y.1.Tminus1 <- c(t(y.panel.1.Tminus1))
  y.tilde <- y.2.T - y.1.Tminus1
  
  
  
  #Construct the RHS variables for the model without a lag:
  #----------------------#
  #   x[i,2] - x[i,1]
  #   x[i,3] - x[i,2]
  #      .        .
  #      .        .
  #      .        .
  #   x[i,T] - x[i,T-1]
  #----------------------#
  x.panel.2.T <- x.panel[, -1]
  x.2.T <- c(t(x.panel.2.T))
  x.panel.1.Tminus1 <- x.panel[, -N.t]
  x.1.Tminus1 <- c(t(x.panel.1.Tminus1))
  X.tilde <- x.2.T - x.1.Tminus1
  
  return(list(X.tilde = X.tilde, y.tilde = y.tilde))
  
  
}#END model.matrices



model.matrices.L <- function(x.panel, y.panel){
  
  
  N.t <- ncol(x.panel)
  
  #Construct the differenced LHS y vector for the model with a lag:
  #----------------------#
  #   y[i,3] - y[i,2]
  #   y[i,4] - y[i,3]
  #      .        .
  #      .        .
  #      .        .
  #   y[i,T] - y[i,T-1]
  #----------------------#
  y.panel.3.T <- y.panel[, -c(1, 2)]
  y.3.T <- c(t(y.panel.3.T))
  y.panel.2.Tminus1 <- y.panel[, -c(1, N.t)]
  y.2.Tminus1 <- c(t(y.panel.2.Tminus1))
  y.tilde.L <- y.3.T - y.2.Tminus1
  
  
  
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
  y.panel.1.Tminus2 <- y.panel[,-c(N.t - 1, N.t)]
  y.1.Tminus2 <- c(t(y.panel.1.Tminus2))
  X1 <- y.2.Tminus1 - y.1.Tminus2
  x.panel.3.T <- x.panel[,-c(1,2)]
  x.3.T <- c(t(x.panel.3.T))
  x.panel.2.Tminus1 <- x.panel[,-c(1, N.t)]
  x.2.Tminus1 <- c(t(x.panel.2.Tminus1))
  X2 <- x.3.T - x.2.Tminus1 
  X.tilde.L <- cbind(X1, X2)
  
  return(list(X.tilde.L = X.tilde.L, y.tilde.L = y.tilde.L))
  
  
}#END model.matrices.L

