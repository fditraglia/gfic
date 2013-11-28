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


set.seed(284)
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

#List of design matrices corresponding to each individual. Rows of each of the N.i dataframes correspond to time periods with NAs "padding" the missing time periods.
data.list <- mapply(function(x.diff, y.diff, y) data.frame(x.diff = x.diff, y.diff = y.diff, y = y), x.diff, y.diff, y, SIMPLIFY = FALSE)


#Function to handle a given individual. Construct as loop to avoid complicated matrix construction and storing really large sparse matrices
XZ.i <- function(data.i){

  out <- 0
  x.diff <- data.i$x.diff
  y.diff <- data.i$y.diff
  y <- data.i$y
  
  #Loop over time periods
  for(j in 3:N.t){
    
    x.tilde.j <- c(y.diff[j - 1], x.diff[j])
    z.j <- c(y[1:(j - 2)], x.diff[j])  
    out <- out + x.tilde.j %*% t(z.j)#THIS IS NOT CORRECT!
    
  }
  
}






