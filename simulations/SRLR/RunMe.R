# This script runs the Short-run versus Long-run simulation exercise from
# Section 6.1 of the paper.

N.SIMS <- 1000
GAMMA <- seq(from = 0.1, to = 0.2, by = 0.01)
R.X.V <- 0.1

source("functions.R")
set.seed(586133)
rmse.T5N250A4 <- RMSE.grid.FAST(g = GAMMA, r = R.X.V, N.sims = N.SIMS, 
                                N.t = 5, N.i = 250)

# Display results
rmse.T5N250A4
