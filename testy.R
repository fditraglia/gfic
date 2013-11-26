library(Rcpp)
library(RcppArmadillo)
setwd("~/Dropbox/GFIC Paper/Arellano Bond Example")
sourceCpp("dgp_cpp.cpp")

set.seed(1234)

foo <- dgpR(a1 = 0.2, a2 = -0.1, g = 0.2, N.i = 100, N.t = 5, burn.in = 10, b = 1, r = 0, theta = 0.2, s.e = 1, s.eta = 1, s.v = 1)

set.seed(1234)

bar <- dgpR(a1 = 0.2, a2 = -0.1, g = 0.2, N.i = 100, N.t = 5, burn.in = 10, b = 1, r = 0, theta = 0.2, s.e = 1, s.eta = 1, s.v = 1)


foobar <- dgp_cpp(a1 = 0.2, a2 = -0.1, g = 0.2, N_i = 100, N_t = 5, burn_in = 10000, b = 1, r = 0, theta = 0.2, s_e = 1, s_eta = 1, s_v = 1)

sourceCpp("rnorm_cpp.cpp")

rnormR <- function(s, N){
  
  epsilon <- s * rnorm(N)
  return(list(e = epsilon))
  
}

set.seed(1234)
fooR <- rnormR(s = 5, N = 10)

set.seed(1234)
barR <- rnormR(s = 5, N = 10)

set.seed(1234)
fooCpp <- rnorm_cpp(s = 5, N = 10)

set.seed(1234)
barCpp <- rnorm_cpp(s = 5, N = 10)


#Now it works!!! Just need to adapt the DGP code for testing purposes...
fooR$e - barR$e 
barR$e - fooCpp$e
fooCpp$e - barCpp$e


#From Dirk - This works as expected! So I should adapt this code...
set.seed(42); rnorm(5)           ## Five N(0,1) draws in R
cppFunction('NumericVector foo(int n) { return rnorm(n); }')
set.seed(42); foo(5)       



