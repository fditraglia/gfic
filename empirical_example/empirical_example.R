#Source:  Baltagi and Levin (1992) and Baltagi, Griffin and Xiong (2000): To pool or not to pool:...
#Description:  Panel Data, 46 U.S. States over the period 1963-1992.
#Variables:
#(1) STATE = State abbreviation.
#(2) YR = YEAR.
#(3) Price per pack of cigarettes.
#(4) Population.
#(5) Population above the age of 16.
#(6) CPI = Consumer price index with (1983=100).
#(7) NDI = Per capita disposable income.
#(8) C = Cigarette sales in packs per capita.
#(9) PIMIN = Minimum price in adjoining states per pack of cigarettes.          

library(MASS)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("fast2SLS_empirical.cpp")

library("plm")
data(Cigar)

data = Cigar[Cigar$year < 86 & Cigar$year > 74, ]

year.f = factor(data$year)
year_dummies = model.matrix(~year.f)

state.f = factor(data$state)
state_dummies = model.matrix(~state.f)

olsY = log((data$sales)/(data$cpi)*(data$pop)/(data$pop16))  # consumption/sales
olsX = log((data$price)/(data$cpi))  # price (all in real terms)

N.t = length(unique(data$year))
N.i = length(unique(data$state))
N.sims = 1  # it's empirical example with the given data

y <- matrix(olsY, ncol = N.t, byrow = TRUE)   # data format is N*T panel data. tested with data$sales
x <- matrix(olsX, ncol = N.t, byrow = TRUE)   

olsW = cbind(log((data$pimin)/(data$cpi)), log((data$ndi)/(data$cpi)), year_dummies[,-1])


W_mat = matrix(0, N.i*(N.t-1), dim(olsW)[2])

for (cc in 1:dim(W_mat)[2]){
  
  w <- matrix(olsW[,cc], ncol = N.t, byrow = TRUE)   
  auxw.tilde <- c(t(w[, -1])) - c(t(w[, -N.t]))  
  W_mat[,cc] <- matrix(t(auxw.tilde), nrow = N.i*(N.t-1), byrow = TRUE)
  
}

Mw <- diag(dim(W_mat)[1]) - W_mat%*%solve(t(W_mat)%*%W_mat)%*%t(W_mat)


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
y.aux <- cbind(c(t(y[, -1])) - c(t(y[, -N.t])))
y.tilde <- c(Mw%*%y.aux)



#Construct the RHS variables for the model without a lag:
#----------------------#
#   x[i,2] - x[i,1]
#   x[i,3] - x[i,2]
#      .        .
#      .        .
#      .        .
#   x[i,T] - x[i,T-1]
#----------------------#
X.aux <- cbind(c(t(x[, -1])) - c(t(x[, -N.t])))
X.tilde <- c(Mw%*%X.aux)

#Construct the differenced LHS y vector for the model with a lag:
#----------------------#
#   y[i,3] - y[i,2]
#   y[i,4] - y[i,3]
#      .        .
#      .        .
#      .        .
#   y[i,T] - y[i,T-1]
#----------------------#
y.aux.L <- matrix(t(y.tilde), ncol = N.i, byrow = FALSE)
y.tilde.L <- c(y.aux.L[-1,])

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
X1.aux.L <- matrix(t(y.tilde), ncol = N.i, byrow = FALSE)
X1.L <- c(X1.aux.L[-(N.t-1),])

X2.aux.L <- matrix(t(X.tilde), ncol = N.i, byrow = FALSE)
X2.L <- c(X2.aux.L[-1,])
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
Z.x.minus.plus <- sparseMatrix(i = rows, j = cols, x = c(Mw%*%cbind(c(t(x[,-N.t])))))
Z.x.plus <- sparseMatrix(i = rows, j = cols, x = c(Mw%*%cbind(c(t(x[,-1])))))



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

Z.y.aux <- matrix(t(c(Mw%*%cbind(c(t(y[, -N.t]))))), ncol = N.i, byrow = FALSE)
Z.y.L <- c(Z.y.aux[-(N.t-1),])

Z.x.minus.aux <- matrix(t(c(Mw%*%cbind(c(t(x[,-N.t]))))), ncol = N.i, byrow = FALSE) # from Z.x.minus.plus above
Z.x.minus.L <- c(Z.x.minus.aux[-1,])

Z.x.aux <-  matrix(t(c(Mw%*%cbind(c(t(x[,-1]))))), ncol = N.i, byrow = FALSE) # from Z.x.plus above
Z.x.L <- c(Z.x.aux[-1,])


rows <- seq_len((N.t - 2) * N.i * N.sims)
cols <- rep(1:(N.t - 2), times = N.i * N.sims)
Z.y <- sparseMatrix(i = rows, j = cols, x = Z.y.L)
Z.x.minus <- sparseMatrix(i = rows, j = cols, x = Z.x.minus.L)
Z.x <- sparseMatrix(i = rows, j = cols, x = Z.x.L)



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

rm(rows, cols, x, y, Z.x, Z.x.minus, Z.y, Z.x.plus, Z.x.minus.plus)


#Now we'll break apart the instrument and model matrices for each replication
block.size.L <- length(y.tilde.L)/N.sims
block.size <- length(y.tilde)/N.sims

seq.block <- seq_len(block.size)
seq.block.L <- seq_len(block.size.L)

#Initialize dataframe to store estimates from each replication
theta.hat <- matrix(NA, nrow = N.sims, ncol = 4)
theta.hat <- data.frame(theta.hat)
names(theta.hat) <- c('LW', 'LS', 'W', 'S')

#Initialize dataframe to store GFIC and J-test from each replication
GFIC <- theta.hat
Biassq <- theta.hat
Var <- theta.hat

k = 1 # N.sims = 1

#Extract the block corresponding to the current replication for the unlagged model
block.k <- (k - 1) * block.size + seq.block
X.tilde.k <- X.tilde[block.k] 
y.tilde.k <- y.tilde[block.k] 
Z.W.k <- Z.W[block.k,]
Z.S.k <- Z.S[block.k,]

#Extract the block corresponding to the current replication for the lagged model
block.L.k <- (k - 1) * block.size.L + seq.block.L
X.tilde.L.k <- X.tilde.L[block.L.k,]
y.tilde.L.k <- y.tilde.L[block.L.k] 
Z.LW.k <- Z.LW[block.L.k,]
Z.LS.k <- Z.LS[block.L.k,]

#Fit 2SLS with centered covariance matrices for everything except LW
reg.LW <- fast2SLS(y.tilde.L.k, X.tilde.L.k, Z.LW.k, center = FALSE, n = N.i)
reg.LS <- fast2SLS(y.tilde.L.k, X.tilde.L.k, Z.LS.k, center = TRUE, n = N.i)
reg.W <- fast2SLS(y.tilde.k, as.matrix(X.tilde.k), as.matrix(Z.W.k), center = TRUE, n = N.i)
reg.S <- fast2SLS(y.tilde.k, as.matrix(X.tilde.k), Z.S.k, center = TRUE, n = N.i)

#Store the 2SLS estimates for each specification
theta.LW <- reg.LW$b[2]
theta.LS <- reg.LS$b[2]
theta.W <- reg.W$b[1]
theta.S <- reg.S$b[1]
theta.hat[k,] <- cbind(theta.LW, theta.LS, theta.W, theta.S)


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
# This is related to expression (60) in the paper
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
Biassq.LW <- 0
Var.LW <- reg.LW$V[2,2]
GFIC.LW <- Biassq.LW + Var.LW
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
Biassq.LS <- (reg.LS$K[2,] %*% inner) %*% t(reg.LS$K)[,2]
Var.LS <- reg.LS$V[2,2]
GFIC.LS <- Biassq.LS + Var.LS

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
xi.3 <- mean(rowSums(Z.x.plus.k) * y.tilde.k)

bias.mat <- matrix(0, nrow = 2, ncol = 2)
bias.mat[1,1] <- delta.squared * xi.3^2

if(tau.squared<0 | delta.squared<0) {bias.mat[2,2] <- delta.squared * xi.1^2 + 2 * tau.delta * xi.1 + tau.squared} else{bias.mat[2,2] <- (sqrt(delta.squared)*xi.1 + sqrt(tau.squared))^2}

bias.mat[1,2] <- bias.mat[2,1] <- delta.squared * xi.1 * xi.3 + tau.delta * xi.3  
i.Tminus1 <- matrix(1, nrow = (N.t - 1))
inner <- bias.mat %x% (i.Tminus1 %*% t(i.Tminus1))
Biassq.S <- reg.S$K %*% inner %*% t(reg.S$K)
Var.S <- reg.S$V
GFIC.S <- Biassq.S + Var.S


#-------------------------------------------------------------
# GFIC for specification W -- incorrect model, exogeneity
# assumption
#-------------------------------------------------------------
# The estimtor of theta in this specification inherits an 
# asymptotic bias from both delta and tau. It also involves
# the quantity xi.3 from above.
#-------------------------------------------------------------
inner <- xi.3^2 * delta.squared * i.Tminus1 %*% t(i.Tminus1)
Biassq.W <- reg.W$K %*% inner %*% t(reg.W$K)
Var.W <- reg.W$V
GFIC.W <- Biassq.W + Var.W

#Store all of the GFIC values
GFIC[k,] <- cbind(GFIC.LW, GFIC.LS, GFIC.W, GFIC.S)
Biassq[k,] <- cbind(Biassq.LW, Biassq.LS, Biassq.W, Biassq.S)
Var[k,] <- cbind(Var.LW, Var.LS, Var.W, Var.S)




#-------------------------------------------------------------
# Selection by GFIC
#-------------------------------------------------------------
which.GFIC <- apply(GFIC, 1, which.min)
theta.GFIC <- theta.hat[cbind(1:N.sims, which.GFIC)] 

theta.hat
theta.GFIC
GFIC
Biassq
Var
