n.reps <- 10000
rho.fine <- seq(0.4, 0.9, 0.01)
va.fine <- seq(0.0, 0.7, 0.1)
rho.coarse <- c(0.4, 0.6, 0.9)
va.coarse <- c(0.1, 0.4, 0.7)
n.grid<- c(100, 200, 500)

params.va.coarse <- expand.grid(n = n.grid,
                                v = va.coarse, 
                                r = rho.fine)
params.rho.coarse <- expand.grid(n = n.grid, 
                                 v = va.fine, 
                                 r = rho.coarse)

results.rho.coarse <- mcmapply(mse_compare_cpp,
                               b = 0.5,
                               rho = params.rho.coarse$r,
                               gamma = params.rho.coarse$v,
                               N = params.rho.coarse$n,
                               T = 2,
                               n_reps = n.reps,
                               mc.cores = nCores)

results.rho.coarse <- cbind(params.rho.coarse, 
                            t(results.rho.coarse))

results.va.coarse <- mcmapply(mse_compare_cpp, 
                              b = 0.5,
                              rho = params.va.coarse$r,
                              gamma = params.va.coarse$v,
                              N = params.va.coarse$n,
                              T = 2,
                              n_reps = n.reps,
                              mc.cores = nCores)

results.va.coarse <- cbind(params.va.coarse, 
                           t(results.va.coarse))

#Output results as R object
results <- list(coarse.rho = results.rho.coarse,
                coarse.va = results.va.coarse)
save(results, file = "./Results/mse_results.Rdata")
#Clean up
rm(n.reps, n.grid)
rm(va.fine, va.coarse)
rm(rho.fine, rho.coarse)
rm(params.va.coarse, params.rho.coarse)
rm(results.rho.coarse, results.va.coarse)
rm(results)