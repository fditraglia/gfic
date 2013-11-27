library(RcppArmadillo)
setwd("~/gfic")
sourceCpp("dgp.cpp")
source("GFIC_arellano_bond.R")



#This setup follows the "base design" of Arellano & Bond 1991 Table 1.
AB_table1_sim <- function(a){
  
  dgp_cpp(a1 = a, a2 = 0, g = 0.1, N_i = 100, N_t = 7, 
                 burn_in = 10, b = 1, r = 0.9, theta = 0, 
                 s_e = sqrt(0.9), s_eta = 1, s_v = 1)
}




