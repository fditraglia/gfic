library(RcppArmadillo)
library(Matrix)
setwd("~/gfic")
sourceCpp("dgp.cpp")
sourceCpp("ABfit.cpp")
source("GFIC_arellano_bond.R")


#Function to replicate one draw from the simulation experiment reported in the first column of Table 1 in Arellano & bond (1991)
testsim <- function(a){
  
  sim.data <- ABsim(a)
  out <- ABfit_cpp(x = sim.data$x, y = sim.data$y)$b
  out <- as.vector(out)
  names(out) <- c("a", "b")
  return(out)
}



#Replicate Column 1 of Table 1 in Arellano & Bond (1991) but with 10000 rather than 100 replications
set.seed(871)
simCpp.2 <- replicate(10000, testsim(0.2))
simCpp.5 <- replicate(10000, testsim(0.5))
simCpp.8 <- replicate(10000, testsim(0.8))

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

lapply(simsCpp, g)
# $`a=0.2`
# a          b
# MEAN  0.16694128 1.00629941
# STDEV 0.07441652 0.06315172
# 
# $`a=0.5`
# a          b
# MEAN  0.44370741 1.00504139
# STDEV 0.09239409 0.06239628
# 
# $`a=0.8`
# a          b
# MEAN  0.77083541 0.99246129
# STDEV 0.05900658 0.06452556


