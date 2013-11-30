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

a.sim.2 <- sapply(seq_len(100), function(i) sim.2[[i]][1,1], simplify = "array")
mean(a.sim.2)
#[1] 0.219452
sd(a.sim.2)
#[1] 0.08032777

b.sim.2 <- sapply(seq_len(100), function(i) sim.2[[i]][2,1], simplify = "array")
mean(b.sim.2)
#[1] 0.9022707
sd(b.sim.2)
#[1] 0.07124266


a.sim.5 <- sapply(seq_len(100), function(i) sim.5[[i]][1,1], simplify = "array")
mean(a.sim.5)
#[1] 0.458137
sd(a.sim.5)
#[1] 0.09586758

b.sim.5 <- sapply(seq_len(100), function(i) sim.5[[i]][2,1], simplify = "array")
mean(b.sim.5)
#[1] 0.9258031
sd(b.sim.5)
#[1] 0.06598999


a.sim.8 <- sapply(seq_len(100), function(i) sim.8[[i]][1,1], simplify = "array")
mean(a.sim.8)
#[1] 0.7582259
sd(a.sim.8)
#[1] 0.05669517

b.sim.8 <- sapply(seq_len(100), function(i) sim.8[[i]][2,1], simplify = "array")
mean(b.sim.8)
#[1] 0.9045425
sd(b.sim.8)
#[1] 0.06523568





