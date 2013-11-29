library(RcppArmadillo)
setwd("~/gfic")
sourceCpp("dgp.cpp")
source("GFIC_arellano_bond.R")



#This setup follows the "base design" of Arellano & Bond 1991 Table 1.
ABsim <- function(a){
  
  dgp_cpp(a1 = a, a2 = 0, g = 0.1, N_i = 100, N_t = 7, 
                 burn_in = 10, b = 1, r = 0.9, theta = 0, 
                 s_e = sqrt(0.9), s_eta = 1, s_v = 1)
}


set.seed(2871)
foo <- ABsim(0.5) #A draw from the DGP of Panel 1, Table 1, A&B 1991
x <- foo$x
y <- foo$y

N.i <- nrow(x)
N.t <- ncol(x)


#Label rows and columns for sanity-checking
rownames(y) <- rownames(x) <- paste('i=', 1:N.i, sep = '')
colnames(y) <- colnames(x) <- paste('t=', 1:N.t, sep = '')

#Take first differences but "pad" with NAs
x.diff <- t(apply(x, 1, function(v) diff(c(NA, v))))
y.diff <- t(apply(y, 1, function(v) diff(c(NA, v))))

#Convert to lists in which each element is an individual
x.diff <- lapply(seq_len(N.i), function(i) x.diff[i,])
y.diff <- lapply(seq_len(N.i), function(i) y.diff[i,])
y <- lapply(seq_len(N.i), function(i) y[i,])

#List of design matrices corresponding to each individual. Rows of each of the N.i dataframes correspond to time periods with NAs "padding" the missing time periods as a sanity check.
data.list <- mapply(function(y.diff, x.diff, y) data.frame(y.diff = y.diff, x.diff = x.diff, y = y), x.diff, y.diff, y, SIMPLIFY = FALSE)

#Instrument Matrix is block-diagonal. Construct using sparse matrices
library(Matrix)

#Instruments for person i
Z.i <- function(i){
  
  data.i <- data.list[[i]]
  Z.i <- lapply(3:nrow(data.i), function(j) with(data.i, c(y[1:(j - 2)], x.diff[j])))
  Z.i <- bdiag(Z.i)
  Z.i <- t(Z.i)
  
}

#Function for a given individual
XZ.i <- function(i){
  
  data.i <- data.list[[i]]
  Z.i <- lapply(3:nrow(data.i), function(j) with(data.i, c(y[1:(j - 2)], x.diff[j])))
  
  Z.i <- bdiag(Z.i)
  Z.i <- t(Z.i)
  
  X.tilde.i <- as.matrix(data.i[3:nrow(data.i), c("x.diff", "y.diff")])
  XZ.i <- crossprod(X.tilde.i, Z.i)
  return(XZ.i)
  
}


XZ <- lapply(seq_len(N.i), XZ.i)
XZ <- Reduce('+', XZ) #Sum over all individuals

H <- bandSparse(5, 5, k = c(0,-1, 1), list(rep(2, 5), rep(-1, 5), rep(-1, 5)))
Z <- lapply(seq_len(N.i), Z.i)
W.n <- lapply(seq_len(N.i), function(i) t(Z[[i]]) %*% H %*% Z[[i]])
W.n <- Reduce('+', W.n)
W.n <- W.n / N.i
W.n <- solve(W.n)

solve(XZ %*% W.n %*% t(XZ))

