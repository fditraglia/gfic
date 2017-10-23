#-----------------------------------------------------------------#
#SCRIPT:  RunMe_SRLR.R
#AUTHOR:  Minsu Chang
#EMAIL:   minsuc@sas.upenn.edu
#DATE:    October 23rd, 2017
#-----------------------------------------------------------------#


N.SIMS <- 1000
GAMMA <- seq(from = 0.1, to = 0.2, by = 0.01)
R.X.V <- 0.1

setwd("/Users/minsuchang/Desktop/ReplicationFiles_CD/Section6/SRLR")
source("GFIC_simulation_SRLR.R")
set.seed(586133)
rmse.T5N250A4 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 5, N.i = 250)

rmse.T5N250A4
