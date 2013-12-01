library(RcppArmadillo)
library(Matrix)
setwd("~/gfic")
sourceCpp("dgp.cpp")
source("GFIC_arellano_bond.R")


#Function to replicate one draw from the simulation experiment reported in the first column of Table 1 in Arellano & Bond (1991)
#NOTE: this code is slow! Takes about 0.7 seconds to run...
test.sim <- function(a){
  
  sim.data <- ABsim(a)
  ABfit(x = sim.data$x, y = sim.data$y)
  
}


set.seed(871)

#Simulations corresponding to alpha = 0.2, 0.5, 0.8
#Each of these takes a little over a minute. I really need to speed things up...
sim.2 <- replicate(100, test.sim(0.2))
sim.5 <- replicate(100, test.sim(0.5))
sim.8 <- replicate(100, test.sim(0.8))

sim.2 <- t(sim.2)
sim.5 <- t(sim.5)
sim.8 <- t(sim.8)

sims <- list(sim.2, sim.5, sim.8)
names(sims) <- c("a=0.2", "a=0.5", "a=0.8")
g <- function(a.b.dataframe){
  
  MEAN <- apply(a.b.dataframe, 2, mean)
  STDEV <- apply(a.b.dataframe, 2, sd)
  
  return(rbind(MEAN, STDEV))
  
}

lapply(sims, g)
#$`a=0.2`
#a          b
#MEAN  0.17861145 1.00066024
#STDEV 0.07318933 0.06464598
# 
# $`a=0.5`
# a         b
# MEAN  0.42672098 1.0128821
# STDEV 0.09163048 0.0632998
# 
# $`a=0.8`
# a         b
# MEAN  0.76670162 0.9899416
# STDEV 0.05911683 0.0675210


#Test out the C++ version of ABfit
set.seed(3728)
testy <- ABsim(0.2)
x <- testy$x
y <- testy$y
sourceCpp("ABfit.cpp")
foo <- ABfit_cpp(x, y)
bCpp <- foo$b
bR <- ABfit(x,y)
all.equal(as.vector(bCpp), as.vector(bR))

microbenchmark(ABfit_cpp(x, y), ABfit(x, y))
# Unit: microseconds
# expr        min         lq      median
# ABfit_cpp(x, y)    414.549    422.486    444.6165
# ABfit(x, y) 681218.531 692722.582 698398.9735
# uq        max neval
# 455.7545    526.159   100
# 706567.4530 790436.044   100




