library(RcppArmadillo)
library(Matrix)
setwd("~/gfic")
sourceCpp("dgp.cpp")
sourceCpp("ABfit.cpp")
source("GFIC_arellano_bond.R")


#Function to replicate one draw from the simulation experiment reported in the first column of Table 1 in Arellano & Bond (1991)
#NOTE: this code is slow! Takes about 0.7 seconds to run...
testsimR <- function(a){
  
  sim.data <- ABsim(a)
  ABfit(x = sim.data$x, y = sim.data$y)
  
}

#C++ version uses ABfit_cpp which is over 1500 times faster than the
#R version of the same
testsimCpp <- function(a){
  
  sim.data <- ABsim(a)
  out <- ABfit_cpp(x = sim.data$x, y = sim.data$y)$b
  return(as.vector(out))
}


#Replicate Column 1 of Table 1 in Arellano & Bond (1991)
set.seed(871)
simR.2 <- replicate(100, testsimR(0.2))
simR.5 <- replicate(100, testsimR(0.5))
simR.8 <- replicate(100, testsimR(0.8))

simR.2 <- t(simR.2)
simR.5 <- t(simR.5)
simR.8 <- t(simR.8)

simsR <- list(simR.2, simR.5, simR.8)
names(simsR) <- c("a=0.2", "a=0.5", "a=0.8")



#C++ version of the simulation
set.seed(871)
simCpp.2 <- replicate(100, testsimCpp(0.2))
simCpp.5 <- replicate(100, testsimCpp(0.5))
simCpp.8 <- replicate(100, testsimCpp(0.8))

simCpp.2 <- t(simCpp.2)
simCpp.5 <- t(simCpp.5)
simCpp.8 <- t(simCpp.8)

simsCpp <- list(simCpp.2, simCpp.5, simCpp.8)
names(simsCpp) <- c("a=0.2", "a=0.5", "a=0.8")




g <- function(a.b.dataframe){
  
  MEAN <- apply(a.b.dataframe, 2, mean)
  STDEV <- apply(a.b.dataframe, 2, sd)
  
  return(rbind(MEAN, STDEV))
  
}



summaryR <- lapply(simsR, g)
summaryR
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


summaryCpp <- lapply(simsCpp, g)
summaryCpp
# $`a=0.2`
# [,1]       [,2]
# MEAN  0.17861145 1.00066024
# STDEV 0.07318933 0.06464598
# 
# $`a=0.5`
# [,1]      [,2]
# MEAN  0.42672098 1.0128821
# STDEV 0.09163048 0.0632998
# 
# $`a=0.8`
# [,1]      [,2]
# MEAN  0.76670162 0.9899416
# STDEV 0.05911683 0.0675210



