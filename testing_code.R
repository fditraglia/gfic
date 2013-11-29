library(RcppArmadillo)
library(Matrix)
setwd("~/gfic")
sourceCpp("dgp.cpp")
source("GFIC_arellano_bond.R")



#This setup follows the "base design" of Arellano & Bond 1991 Table 1.
ABsim <- function(a){
  
  dgp_cpp(a1 = a, a2 = 0, g = 0.1, N_i = 100, N_t = 7, 
                 burn_in = 10, b = 1, r = 0.9, theta = 0, 
                 s_e = sqrt(0.9), s_eta = 1, s_v = 1)
}


#set.seed(871)
foo <- ABsim(0.5) #A draw from the DGP of Panel 1, Table 1, A&B 1991
x <- foo$x
y <- foo$y

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






