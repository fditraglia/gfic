#-----------------------------------------------------------------#
#SCRIPT:  simulation_2013_02_15.R
#AUTHOR:  Frank DiTraglia
#EMAIL:   fditra@sas.upenn.edu
#DATE:    Marth 8th, 2013
#-----------------------------------------------------------------#


#-----------------------------------------------------------------#
# MACROS
#-----------------------------------------------------------------#

THIS.DIR <- 'simulation_2013_03_08'
N.SIMS <- 2000
GAMMA <- seq(from = 0, to = 0.2, by = 0.005)
R.X.V <- seq(from = 0, to = 0.2, by = 0.005)


#-----------------------------------------------------------------#
# END MACROS
#-----------------------------------------------------------------#





#-----------------------------------------------------------------#
# DEFINE DIRECTORY STRUCTURE
#-----------------------------------------------------------------#

#Path to Dropbox
ROOT.DIR <- '~/Dropbox'

#Set appropriate working directory
WORKING.DIR <- paste(ROOT.DIR, '/Phd_Dissertation/Generalized FIC Paper/simulations/', THIS.DIR, sep = '')
FUNCTIONS.DIR <- paste(ROOT.DIR, '/Phd_Dissertation/Generalized FIC Paper/simulations/functions', sep = '')

#Clean up
rm(ROOT.DIR)

#-----------------------------------------------------------------#
# END DEFINE DIRECTORY STRUCTURE
#-----------------------------------------------------------------#




#-----------------------------------------------------------------#
# LOAD FUNCTIONS 
#-----------------------------------------------------------------#

setwd(FUNCTIONS.DIR)
source('GFIC_simulation_v2.R')

#-----------------------------------------------------------------#
# END LOAD FUNCTIONS
#-----------------------------------------------------------------#




#-----------------------------------------------------------------#
# RUN SIMULATION AND WRITE RESULTS TO FILE
#-----------------------------------------------------------------#

set.seed(1445)

setwd(WORKING.DIR)

#----------- N = 250
rmse.T4.N250 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 4, N.i = 250)
write.csv(rmse.T4.N250, file = 'rmse_T4_N250.csv', row.names = FALSE)

rmse.T5.N250 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 5, N.i = 250)
write.csv(rmse.T5.N250, file = 'rmse_T5_N250.csv', row.names = FALSE)


#rmse.T6.N250 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 6, N.i = 250)
#write.csv(rmse.T6.N250, file = 'rmse_T6_N250.csv', row.names = FALSE)


#----------- N = 500
rmse.T4.N500 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 4, N.i = 500)
write.csv(rmse.T4.N500, file = 'rmse_T4_N500.csv', row.names = FALSE)

rmse.T5.N500 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 5, N.i = 500)
write.csv(rmse.T5.N500, file = 'rmse_T5_N500.csv', row.names = FALSE)

#rmse.T6.N500 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 6, N.i = 500)
#write.csv(rmse.T6.N500, file = 'rmse_T6_N500.csv', row.names = FALSE)


#----------- N = 1000
#rmse.T4.N1000 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 4, N.i = 1000)
#write.csv(rmse.T4.N1000, file = 'rmse_T4_N1000.csv', row.names = FALSE)

#rmse.T5.N1000 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 5, N.i = 1000)
#write.csv(rmse.T5.N1000, file = 'rmse_T5_N1000.csv', row.names = FALSE)

#rmse.T6.N1000 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, N.t = 6, N.i = 1000)
#write.csv(rmse.T6.N1000, file = 'rmse_T6_N1000.csv', row.names = FALSE)


#-----------------------------------------------------------------#
# END SCRIPT
#-----------------------------------------------------------------#
