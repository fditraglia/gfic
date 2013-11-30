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
b.sim.2 <- sapply(seq_len(100), function(i) sim.2[[i]][2,1], simplify = "array")
a.sim.5 <- sapply(seq_len(100), function(i) sim.5[[i]][1,1], simplify = "array")
b.sim.5 <- sapply(seq_len(100), function(i) sim.5[[i]][2,1], simplify = "array")
a.sim.8 <- sapply(seq_len(100), function(i) sim.8[[i]][1,1], simplify = "array")
b.sim.8 <- sapply(seq_len(100), function(i) sim.8[[i]][2,1], simplify = "array")


names(a.sim.2) <- names(b.sim.2) <- names(a.sim.5) <- names(b.sim.5) <- names(a.sim.8) <- names(b.sim.8) <- NULL

sim.2 <- data.frame(a = a.sim.2, b = b.sim.2)
sim.5 <- data.frame(a = a.sim.5, b = b.sim.5)
sim.8 <- data.frame(a = a.sim.8, b = b.sim.8)

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







