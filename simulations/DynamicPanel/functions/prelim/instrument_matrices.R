#-----------------------------------------------------------------#
#FILENAME:  instrument_matrices.R
#AUTHOR:    Frank DiTraglia
#EMAIL:     fditraglia@gmail.com 
#DATE:      May 31st, 2012
#-----------------------------------------------------------------#
#DESCRIPTION:   
#
#       This file contains several functions for setting up instrument matrices to estimate a simple, dynamic panel model with correlated individual effects:
#
#     y[i,t] = gamma * y[i,t-1] + theta * x[i,t] + eta[i] + v[i,t]
#
#where x is a scalar. Estimation is in first differences to remove the individual effect eta[i]. The instrument matrices correspond to four estimators:
#
# Estimator   Lag y?  Assume x Strict Exog.?  Time t instruments
# ---------   ------  ----------------------  ------------------
# b.LW        YES     NO (assume weak exog.)  y[i,t-2], x[i,t-1]
# b.LS        YES     YES                     y[i,t-2], x[i,t-1], x[i,t]     
# b.W         NO      NO (assume weak exog.)  x[i,t-1]    
# b.S         NO      YES                     x[i,t-1], x[i,t]
#
#and rather than using all valid instruments, include only those that are by definition strong. These instrument sets are roughly those proposed by Arellano and Bond (1991).
#-----------------------------------------------------------------#


#Load package for constructing and manipulating sparse matrices
library(Matrix)




#---------------------------------------------------------------------#
#                     FUNCTION instruments.S
#---------------------------------------------------------------------#
#DESCRIPTION:
#   
#     This function constructs the instrument matrix for the estimator b.S described in the documentation above.
#---------------------------------------------------------------------#
#ARGUMENTS:
#
#             x.panel     matrix containing values of regressor x
#                           rows are individuals, columns time periods
#---------------------------------------------------------------------#
#RETURNS:
#             sparse matrix of individual instrument matrices
#---------------------------------------------------------------------#
instruments.S <- function(x.panel)
{
  
  N.t <- ncol(x.panel) #Number of time periods
  N.i <- nrow(x.panel) #Number of individuals
  
  x1 <- c(t(x.panel[,-N.t]))
  x2 <- c(t(x.panel[,-1]))
  values <- c(x1, x2)
  
  row.index <- rep(seq_len(N.i * (N.t - 1)), times = 2)
  
  index <- seq_len(N.t - 1)
  col1 <- rep(index * 2 - 1, times = N.i)
  col2 <- rep(index * 2, times = N.i)
  col.index <- c(col1, col2)
  
  as.matrix(sparseMatrix(i = row.index, j = col.index, x = values))
  
  
}#END instruments.S

#---------------------------------------------------------------------#
#                   END FUNCTION instruments.S
#---------------------------------------------------------------------#






#---------------------------------------------------------------------#
#                     FUNCTION instruments.W
#---------------------------------------------------------------------#
#DESCRIPTION:
#   
#     This function constructs the instrument matrix for the estimator b.W described in the documentation above.
#---------------------------------------------------------------------#
#ARGUMENTS:
#
#             x.panel     matrix containing values of regressor x
#                           rows are individuals, columns time periods
#---------------------------------------------------------------------#
#RETURNS:
#             sparse matrix of individual instrument matrices
#---------------------------------------------------------------------#
instruments.W <- function(x.panel)
{

  N.t <- ncol(x.panel) #Number of time periods
  N.i <- nrow(x.panel) #Number of individuals
  
  
  #Vector of nonzero entries that to fill the full instrument matrix
  values <- c(t(x.panel[,-N.t]))
  
  
  #Row and column indices for entries of final sparse matrix
  row.index <- seq_len(N.i * (N.t - 1))
  col.index <- rep(seq_len(N.t - 1), times = N.i)
  
  #Finally, assemble the matrix of instruments
  as.matrix(sparseMatrix(i = row.index, j = col.index, x = values))
  
  
}#END instruments.W

#---------------------------------------------------------------------#
#                   END FUNCTION instruments.W
#---------------------------------------------------------------------#








#---------------------------------------------------------------------#
#                     FUNCTION instruments.LW
#---------------------------------------------------------------------#
#DESCRIPTION:
#   
#     This function constructs the instrument matrix for the estimator b.LW described in the documentation above.
#---------------------------------------------------------------------#
#ARGUMENTS:
#
#             x.panel     matrix containing values of regressor x.
#                           Rows are individuals, columns time periods.
#             y.panel     matrix containing values of outcome y.
#                           Rows are individuals, columns time periods.
#---------------------------------------------------------------------#
#RETURNS:
#             sparse matrix of individual instrument matrices
#---------------------------------------------------------------------#
instruments.LW <- function(x.panel, y.panel)
{
  
  N.t <- ncol(x.panel) #Number of time periods
  N.i <- nrow(x.panel) #Number of individuals
  
  x <- c(t(x.panel[,-c(1, N.t)]))
  y <- c(t(y.panel[,-c(N.t, N.t - 1)]))
  
  values <- c(y, x)

  row.index <- rep(seq_len(N.i * (N.t - 2)), times = 2)
  
  index <- seq_len(N.t - 2)
  col1 <- rep(index * 2 - 1, times = N.i)
  col2 <- rep(index * 2, times = N.i)
  col.index <- c(col1, col2)
  
  as.matrix(sparseMatrix(i = row.index, j = col.index, x = values))
  
  
}#END instruments.LW

#---------------------------------------------------------------------#
#                   END FUNCTION instruments.LW
#---------------------------------------------------------------------#






#---------------------------------------------------------------------#
#                     FUNCTION instruments.LS
#---------------------------------------------------------------------#
#DESCRIPTION:
#   
#     This function constructs the instrument matrix for the estimator b.LS described in the documentation above.
#---------------------------------------------------------------------#
#ARGUMENTS:
#
#             x.panel     matrix containing values of regressor x.
#                           Rows are individuals, columns time periods.
#             y.panel     matrix containing values of outcome y.
#                           Rows are individuals, columns time periods.
#---------------------------------------------------------------------#
#RETURNS:
#             sparse matrix of individual instrument matrices
#---------------------------------------------------------------------#
instruments.LS <- function(x.panel, y.panel)
{
  
  N.t <- ncol(x.panel) #Number of time periods
  N.i <- nrow(x.panel) #Number of individuals
  
  x1 <- c(t(x.panel[,-c(1, N.t)]))
  x2 <- c(t(x.panel[,-c(1, 2)]))
  y <- c(t(y.panel[,-c(N.t, N.t - 1)]))
  
  values <- c(y, x1, x2)
  
  row.index <- rep(seq_len(N.i * (N.t - 2)), times = 3)
  
  index <- seq_len(N.t - 2)
  col1 <- rep(index * 3 - 2, times = N.i)
  col2 <- rep(index * 3 - 1, times = N.i)
  col3 <- rep(index * 3, times = N.i)
  col.index <- c(col1, col2, col3)
  
  as.matrix(sparseMatrix(i = row.index, j = col.index, x = values))
  
  
}#END instruments.LS

#---------------------------------------------------------------------#
#                   END FUNCTION instruments.LS
#---------------------------------------------------------------------#




#Test these functions with a simple example
#Time <- 3
#N <- 5000
#A <- 1:N
#B <- (1:Time)/10
#x.panel <- outer(A, B, '+')
#y.panel <- outer(10*A, B, '+')

#system.time({
#Z.LS <- instruments.LS(x.panel, y.panel)
#Z.LW <- instruments.LW(x.panel, y.panel)
#Z.S <- instruments.S(x.panel)
#Z.W <- instruments.W(x.panel)
#})

