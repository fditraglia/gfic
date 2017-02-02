#####LOAD DATASET
#Set appropriate working directory
WORKING.DIR <- paste('~/Dropbox/Phd_Dissertation/Generalized FIC Paper/simulations/simulation_2013_01_18', sep = '')

setwd(WORKING.DIR)

T4.N250 <- read.csv('rmse_T4_N250.csv')
T5.N250 <- read.csv('rmse_T5_N250.csv')
#T6.N250 <- read.csv('rmse_T6_N250.csv')

T4.N500 <- read.csv('rmse_T4_N500.csv')
T5.N500 <- read.csv('rmse_T5_N500.csv')
#T6.N500 <- read.csv('rmse_T6_N500.csv')

#T4.N1000 <- read.csv('rmse_T4_N1000.csv')
#T5.N1000 <- read.csv('rmse_T5_N1000.csv')
#T6.N1000 <- read.csv('rmse_T6_N1000.csv')



T4.N250 <- data.frame(T4.N250, N.t = rep(4, nrow(T4.N250)), N.i = rep(250, nrow(T4.N250)))
T5.N250 <- data.frame(T5.N250, N.t = rep(5, nrow(T5.N250)), N.i = rep(250, nrow(T5.N250)))
#T6.N250 <- data.frame(T6.N250, N.t = rep(6, nrow(T6.N250)), N.i = rep(250, nrow(T6.N250)))

T4.N500 <- data.frame(T4.N500, N.t = rep(4, nrow(T4.N500)), N.i = rep(500, nrow(T4.N500)))
T5.N500 <- data.frame(T5.N500, N.t = rep(5, nrow(T5.N500)), N.i = rep(500, nrow(T5.N500)))
#T6.N500 <- data.frame(T6.N500, N.t = rep(6, nrow(T6.N500)), N.i = rep(500, nrow(T6.N500)))

#T4.N1000 <- data.frame(T4.N1000, N.t = rep(4, nrow(T4.N1000)), N.i = rep(1000, nrow(T4.N1000)))
#T5.N1000 <- data.frame(T5.N1000, N.t = rep(5, nrow(T5.N1000)), N.i = rep(1000, nrow(T5.N1000)))
#T6.N1000 <- data.frame(T6.N1000, N.t = rep(6, nrow(T6.N1000)), N.i = rep(1000, nrow(T6.N1000)))



results <- rbind(T4.N250, T5.N250, T4.N500, T5.N500)
#results <- rbind(T4.N250, T5.N250, T6.N250, T4.N500, T5.N500, T6.N500, T4.N1000, T5.N1000, T6.N1000)


#Save the numeric versions for reordering the factors
order.N <- results$N.i
order.t <- results$N.t

results$N.i <- as.factor(paste('N = ', results$N.i, sep = ''))
results$N.t <- as.factor(paste('T = ', results$N.t, sep = ''))

#Set factor levels manually so that ggplot2 puts the facets in correct order
results$N.i <- reorder(results$N.i, order.N)
results$N.t <- reorder(results$N.t, order.t)

rmse <- results[,c('LW', 'LS', 'W', 'S')]
optimal <- as.factor(names(rmse)[apply(rmse, 1, which.min)])
min.rmse <- apply(rmse, 1, min)


f <- function(x){
  x.1 <- min(x)
  x.2 <- min(x[-c(which.min(x))])
  x.2 - x.1
}#END f

f2 <- function(x){
  x.1 <- min(x)
  x.2 <- min(x[-c(which.min(x))])
  (x.2 - x.1)/x.2 * 100
}#END f2


#Comparison to best specification at this point of paramter space
rmse.diff <- apply(rmse, 1, f)
rmse.percent <- apply(rmse, 1, f2)

#Comparisons to correct specification
rmse.diff.vs.LW <- apply(rmse - rmse$LW, 1, min)
rmse.percent.vs.LW <- apply((rmse - rmse$LW)/rmse$LW, 1, min)

results <- data.frame(results, optimal, min.rmse, rmse.diff, rmse.percent, rmse.diff.vs.LW, rmse.percent.vs.LW)

library(ggplot2)


#Compare RMSE of best specification against LW
v <- ggplot(results, aes(r.x.v, gamma))

#How much better is best specification than LW?
v + geom_raster(aes(fill = rmse.percent.vs.LW * -100)) + facet_grid(N.i ~ N.t) +  scale_fill_gradientn(colours = terrain.colors(4), guide = "colorbar") + ylab(expression(gamma)) + xlab(expression(rho[x*v])) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + labs(fill = '') #+ opts(panel.background = theme_blank(), panel_border = theme_blank(), strip.background = theme_blank()) 
                
#Color-blind palette
cbPalette <- c("#E69F00", "#009E73", "#999999", "#56B4E9")

#Which estimator is optimatal in each region of the parameter space?
v + geom_tile(aes(fill = optimal, alpha = rmse.percent))  + facet_grid(N.i ~ N.t) + ylab(expression(gamma)) + xlab(expression(rho[x*v])) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + labs(fill = '') + scale_fill_manual(values = cbPalette) + scale_alpha_continuous(range = c(0.2, 1), guide = 'none', trans = 'sqrt') #+ opts(panel.background = element_blank(), panel_border = element_blank(), strip.background = element_blank()) 



