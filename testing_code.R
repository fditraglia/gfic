library(RcppArmadillo)
library(Matrix)
setwd("~/gfic")
sourceCpp("dgp.cpp")
source("GFIC_arellano_bond.R")


set.seed(871)
foo <- ABsim(0.5) #A draw from the DGP of Panel 1, Table 1, A&B 1991
x <- foo$x
y <- foo$y
ABfit(x, y)



