#-----------------------------------------------------------------#
# Parameter values and sample size
#-----------------------------------------------------------------#
N.SIMS <- 2000
GAMMA <- seq(from = 0, to = 0.2, by = 0.005)
R.X.V <- seq(from = 0, to = 0.2, by = 0.005)

#-----------------------------------------------------------------#
# Load simulation functions
#-----------------------------------------------------------------#
source('functions.R')

#-----------------------------------------------------------------#
# Run simulation and write results to file.
#-----------------------------------------------------------------#

set.seed(1445)

setwd('./results/')

#----------- N = 250, T = 4
rmse.T4.N250 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, 
                               N.t = 4, N.i = 250)
write.csv(rmse.T4.N250, file = 'rmse_T4_N250.csv', row.names = FALSE)

#----------- N = 250, T = 5
rmse.T5.N250 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, 
                               N.t = 5, N.i = 250)
write.csv(rmse.T5.N250, file = 'rmse_T5_N250.csv', row.names = FALSE)


#----------- N = 500, T = 4
rmse.T4.N500 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, 
                               N.t = 4, N.i = 500)
write.csv(rmse.T4.N500, file = 'rmse_T4_N500.csv', row.names = FALSE)

#----------- N = 500, T = 5
rmse.T5.N500 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, 
                               N.t = 5, N.i = 500)
write.csv(rmse.T5.N500, file = 'rmse_T5_N500.csv', row.names = FALSE)

# Clean up
rm(list = ls())