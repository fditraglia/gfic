library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(tikzDevice)

setwd("/Users/minsuchang/Desktop/ReplicationFiles_CD/Appendix") 
sourceCpp("functionsREFE.cpp")

set.seed(1928)

#Set the number of cores for mclapply.
#If you're using Windows you can't
#make use of this function, so set 
#nCores to 1 and mclapply will revert 
#to mapply. If you're running Linux or 
#you can use mclapply. If you have k 
#cores, set nCores to k + 1 for best 
#performance.
nCores <- 4

#Calculate RMSE results for the simulation
#and store them in ./Results/rmse_results.Rdata
source("mse_calculationsREFE.R")

#Open ./Results/rmse_results.Rdata
#plot the results 
#and export in tikz format to ./Results/
source("rmse_plotsREFE.R")
