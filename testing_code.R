library(RcppArmadillo)
setwd("~/gfic")
sourceCpp("dgp.cpp")
source("GFIC_arellano_bond.R")


set.seed(837)
fooR <- dgpR(a1 = 0.3, a2 = -0.1, g = 0.1, N.i = 25, N.t = 4, 
             burn.in = 10, b = 1, r = 0.1, theta = 0.2, s.e = 1, 
             s.eta = 1, s.v = 1)

set.seed(837)
fooCpp <- dgp_cpp(a1 = 0.3, a2 = -0.1, g = 0.1, N_i = 25, N_t = 4, 
                  burn_in = 10, b = 1, r = 0.1, theta = 0.2, s_e = 1, 
                  s_eta = 1, s_v = 1)

all.equal(fooR$x, fooCpp$x)
all.equal(fooR$y, fooCpp$y)


all.equal(as.vector(fooCpp$eta), fooR$eta)
all.equal(fooCpp$v, fooR$v)
all.equal(fooCpp$epsilon, fooR$epsilon)



library(microbenchmark)
