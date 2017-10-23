setwd("~/gfic/simulations/DynamicPanel/results/")

T4_N250 <- read.csv("rmse_T4_N250.csv", stringsAsFactors = FALSE) 
T4_N250 <- data.frame(T4_N250, nt = rep(4, nrow(T4_N250)), 
                      ni = rep(250, nrow(T4_N250)))

T5_N250 <- read.csv("rmse_T5_N250.csv", stringsAsFactors = FALSE) 
T5_N250 <- data.frame(T5_N250, nt = rep(5, nrow(T5_N250)), 
                      ni = rep(250, nrow(T5_N250)))

T4_N500 <- read.csv("rmse_T4_N500.csv", stringsAsFactors = FALSE) 
T4_N500 <- data.frame(T4_N500, nt = rep(4, nrow(T4_N500)), 
                      ni = rep(500, nrow(T4_N500)))

T5_N500 <- read.csv("rmse_T5_N500.csv", stringsAsFactors = FALSE) 
T5_N500 <- data.frame(T5_N500, nt = rep(5, nrow(T5_N500)), 
                      ni = rep(500, nrow(T5_N500)))

results <- rbind(T4_N250, T5_N250, T4_N500, T5_N500)

# Change names of LW and W to match paper
names(results)[which(names(results) == 'LW')] <- 'LP'
names(results)[which(names(results) == 'W')] <- 'P'

rm(T4_N250, T4_N500, T5_N250, T5_N500)

#Save the numeric versions for reordering the factors
orderN <- results$ni
orderT <- results$nt

results$ni <- as.factor(paste('N = ', results$ni, sep = ''))
results$nt <- as.factor(paste('T = ', results$nt, sep = ''))

#Set factor levels manually so that ggplot2 puts the facets in correct order
results$ni <- reorder(results$ni, orderN)
results$nt <- reorder(results$nt, orderT)

rmse <- results[, c('LP', 'LS', 'P', 'S')]
optimal <- as.factor(names(rmse)[apply(rmse, 1, which.min)])
min_rmse <- apply(rmse, 1, min)


# Compare second best estimator to best at given point in parameter space
diff2minus1 <- function(x) { 
  x <- sort(x) 
  x[2] - x[1]
}

ratio2to1 <- function(x) { 
  x <- sort(x)  
  100 * (x[2] - x[1]) / x[2] 
}

rmse_diff <- apply(rmse, 1, diff2minus1)
rmse_percent <- apply(rmse, 1, ratio2to1)

#Comparisons to correct specification
rmse_diff_vs_LP <- apply(rmse - rmse$LP, 1, min)
rmse_percent_vs_LP <- apply((rmse - rmse$LP)/rmse$LP, 1, min)

results <- data.frame(results, optimal, min_rmse, rmse_diff, rmse_percent, 
                      rmse_diff_vs_LP, rmse_percent_vs_LP)

library(ggplot2)


#Compare RMSE of best specification against LP
v <- ggplot(results, aes(r.x.v, gamma))

#testy <- colorRampPalette(c("white", "green1", "green2", "green4"))(4)
#testy <- cm.colors(4)  
library(RColorBrewer)
testy <- brewer.pal(9,"YlGnBu")

#How much better is best specification than LP?
v + geom_raster(aes(fill = round(rmse_percent_vs_LP * -100))) + facet_grid(ni ~ nt) +  
  scale_fill_gradientn(colours = testy, guide = "colorbar") + 
  ylab(expression(gamma)) + xlab(expression(rho[x*v])) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  labs(fill = '') 
#+ opts(panel.background = theme_blank(), panel_border = theme_blank(), strip.background = theme_blank()) 
                
#Color-blind palette
cbPalette <- c("#E69F00", "#009E73", "#999999", "#56B4E9")

#Which estimator is optimal in each region of the parameter space?
pdf(file = 'Dpanel_oracle.pdf', width = 6, height = 5)
v + geom_tile(aes(fill = optimal, alpha = rmse_percent)) + 
  facet_grid(ni ~ nt) + ylab(expression(gamma)) + xlab(expression(rho[x*v])) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  labs(fill = '') + scale_fill_manual(values = cbPalette) + 
  scale_alpha_continuous(range = c(0.2, 1), guide = 'none', trans = 'sqrt') 
dev.off()
#+ opts(panel.background = element_blank(), panel_border = element_blank(), strip.background = element_blank()) 

# Compare GFIC to: Oracle, Conservative, Aggressive, competing criteria...
results$oracle <- apply(rmse, 1, min)

# Calculate the appropriate comparison (relative to oracle?), plot color "diverging"
# color scale that has white for zero, blue for GFIC better and red for alternative better
diverge <- brewer.pal(9, "RdBu")

# Relative RMSE of GFIC to the true DGP (LP)
pdf(file = 'Dpanel_GFIC_RMSE_rel_LP.pdf', width = 6, height = 5)
v + geom_raster(aes(fill = with(results, 100 * (GFIC - LP) / LP ))) + 
  facet_grid(ni ~ nt) +  
  scale_fill_gradientn(colours = diverge, guide = "colorbar") + 
  ylab(expression(gamma)) + xlab(expression(rho[x*v])) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  labs(fill = '') 
dev.off()

# Relative RMSE of GFIC to the oracle
pdf(file = 'Dpanel_GFIC_RMSE_rel_oracle.pdf', width = 6, height = 5)
v + geom_raster(aes(fill = with(results, 100 * (GFIC - oracle) / oracle))) + 
  facet_grid(ni ~ nt) +  
  scale_fill_gradientn(colours = testy, guide = "colorbar") + 
  ylab(expression(gamma)) + xlab(expression(rho[x*v])) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  labs(fill = '') 
dev.off()

# Relative RMSE of GFIC to GMM-BIC
pdf(file = 'Dpanel_GFIC_RMSE_rel_BIC.pdf', width = 6, height = 5)
v + geom_raster(aes(fill = with(results, 100 * (GFIC - BIC) / BIC))) + 
  facet_grid(ni ~ nt) +  
  scale_fill_gradientn(colours = diverge, guide = "colorbar") + 
  ylab(expression(gamma)) + xlab(expression(rho[x*v])) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  labs(fill = '') 
dev.off()


# Relative RMSE of GFIC to GMM-HQ
pdf(file = 'Dpanel_GFIC_RMSE_rel_HQ.pdf', width = 6, height = 5)
v + geom_raster(aes(fill = with(results, 100 * (GFIC - HQ) / HQ))) + 
  facet_grid(ni ~ nt) +  
  scale_fill_gradientn(colours = diverge, guide = "colorbar") + 
  ylab(expression(gamma)) + xlab(expression(rho[x*v])) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  labs(fill = '') 
dev.off()

# Relative RMSE of GFIC to GMM-AIC
pdf(file = 'Dpanel_GFIC_RMSE_rel_AIC.pdf', width = 6, height = 5)
v + geom_raster(aes(fill = with(results, 100 * (GFIC - AIC) / AIC))) + 
  facet_grid(ni ~ nt) +  
  scale_fill_gradientn(colours = diverge, guide = "colorbar") + 
  ylab(expression(gamma)) + xlab(expression(rho[x*v])) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  labs(fill = '') 
dev.off()

# Relative RMSE of GFIC to J10
pdf(file = 'Dpanel_GFIC_RMSE_rel_J10.pdf', width = 6, height = 5)
v + geom_raster(aes(fill = with(results, 100 * (GFIC - J10) / J10))) + 
  facet_grid(ni ~ nt) +  
  scale_fill_gradientn(colours = diverge, guide = "colorbar") + 
  ylab(expression(gamma)) + xlab(expression(rho[x*v])) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  labs(fill = '') 
dev.off()

# Relative RMSE of GFIC to J5
pdf(file = 'Dpanel_GFIC_RMSE_rel_J5.pdf', width = 6, height = 5)
v + geom_raster(aes(fill = with(results, 100 * (GFIC - J5) / J5))) + 
  facet_grid(ni ~ nt) +  
  scale_fill_gradientn(colours = diverge, guide = "colorbar") + 
  ylab(expression(gamma)) + xlab(expression(rho[x*v])) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  labs(fill = '') 
dev.off()

# Clean up
rm(list = ls())