# This script runs the Random versus Fixed Effects simulation from the 
# described in the online appendix to the paper. Before running this script,
# create a subdirectory within the directory where it is stored called 'results'

library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(tikzDevice)

# Compile C++ functions for simulation
sourceCpp("functions.cpp")

set.seed(1928)

#Set the number of cores for mclapply.
#If you're using Windows you can't
#make use of this function, so set 
#nCores to 1 and mclapply will revert 
#to mapply. If you're running Linux or Mac
#you can use mclapply. If you have k 
#cores, set nCores to k + 1 for best 
#performance.
nCores <- 4

#Calculate RMSE results for the simulation
#and store them in ./results/rmse_results.Rdata
source("calculate_RMSE.R")

#Open ./results/rmse_results.Rdata
#plot the results 
#and export in tikz format to ./results/
source("make_plots.R")
