# Probably should just ditch these tables and present the images instead...

setwd("~/Dropbox/GFIC_annex/simulations/simulation_2013_03_08/")

T4_N250 <- read.csv("rmse_T4_N250.csv", stringsAsFactors = FALSE) 
T4_N250 <- data.frame(T4_N250, nt = rep(4, nrow(T4_N250)), ni = rep(250, nrow(T4_N250)))

T5_N250 <- read.csv("rmse_T5_N250.csv", stringsAsFactors = FALSE) 
T5_N250 <- data.frame(T5_N250, nt = rep(5, nrow(T5_N250)), ni = rep(250, nrow(T5_N250)))

T4_N500 <- read.csv("rmse_T4_N500.csv", stringsAsFactors = FALSE) 
T4_N500 <- data.frame(T4_N500, nt = rep(4, nrow(T4_N500)), ni = rep(500, nrow(T4_N500)))

T5_N500 <- read.csv("rmse_T5_N500.csv", stringsAsFactors = FALSE) 
T5_N500 <- data.frame(T5_N500, nt = rep(5, nrow(T5_N500)), ni = rep(500, nrow(T5_N500)))

results <- rbind(T4_N250, T5_N250, T4_N500, T5_N500)
names(results)[1:2] <- c('g', 'r')
results$oracle <- apply(results[, c('LW', 'LS', 'W', 'S')], 1, min)

rm(T4_N250, T4_N500, T5_N250, T5_N500)

# Coarsen the simulation grid over g and r
g_coarse <- c(0, 0.05, 0.1, 0.15)
r_coarse <- c(0, 0.05, 0.1, 0.15)
coarse <- subset(results, (g %in% g_coarse) & (r %in% r_coarse))

# Use xtabs and ftable to make "flat" tables of RMSE in THOUSANDTHS
# i.e. 1000 * RMSE
make_tab <- function(col){
  tab <- xtabs(1000 * round(col, 3) ~ g + r + ni + nt, coarse)
  names(dimnames(tab)) <- c('gamma', 'rho', 'N', 'T')
  ftable(tab, row.vars = c('T', 'N', 'gamma'), col.vars = c('rho'))
}

LW <- make_tab(coarse$LW)
LS <- make_tab(coarse$LS)
W <- make_tab(coarse$W)
S <- make_tab(coarse$S)

# I think we just need to show GFIC, Oracle, and some competitors...
# Maybe show the best competitor?
oracle <- make_tab(coarse$oracle)
GFIC <- make_tab(coarse$GFIC)
J10 <- make_tab(coarse$J10)
J5 <- make_tab(coarse$J5)
AIC <- make_tab(coarse$AIC)
BIC <- make_tab(coarse$BIC)
HQ <- make_tab(coarse$HQ)

attr(oracle, "row.vars")
attr(oracle, "col.vars")

# Need to write code to put together the table in a reasonable way...

# Possibly easiest to output these as a bunch of separate tables?
library(xtable)
print(xtable(format(oracle)))


# Looks good! Now just output these in some kind of reasonable format and 
# decide which figures to include.

# May need to use either Hmisc as below or xtable...

# library(Hmisc)
# 
# 
# outfilename <- 'GFIC_summary_tables.tex'
#  
# 
# #Assign to variable table to keep it from trying to compile
# table <- latex(round(mean.vs), file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = FALSE)
# 
# #Append the remaining tables
# table <- latex(round(max.vs), file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE) 
# 
# table <- latex(round(mean.rmse,3), file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE) 
# 
# table <- latex(round(max.rmse,3), file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)
# 
# table <- latex(raw.T4.N250, file = 'GFIC_pointwise_tables_T4_N250.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)
# 
# table <- latex(raw.T5.N250, file = 'GFIC_pointwise_tables_T5_N250.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)
# 
# 
# table <- latex(raw.T4.N500, file = 'GFIC_pointwise_tables_T4_N500.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)
# 
# table <- latex(raw.T5.N500, file = 'GFIC_pointwise_tables_T5_N500.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)
# 
