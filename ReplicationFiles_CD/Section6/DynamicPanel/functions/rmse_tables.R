#####LOAD DATASET
#Set appropriate working directory
WORKING.DIR <- paste('~/Dropbox/Phd_Dissertation/Generalized FIC Paper/simulations/simulation_2013_03_08', sep = '')

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
names(results)[1:2] <- c('g', 'r')


rmse <- results[,c('LW', 'LS', 'W', 'S')]
optimal <- as.factor(names(rmse)[apply(rmse, 1, which.min)])
true.min <- apply(rmse, 1, min)


vs.min <- cbind(results[,c(1:2, 13:14)], 100 * (results[,c(3,7:12)] - true.min)/true.min)



params <- T4.N250[,1:2]
names(params) <- c('g', 'r')

f <- function(x, k = 0){tapply(x, params, round, digits = k)}

#Pointwise RMSE
raw.T4.N250 <- subset(results, N.t == 4 & N.i == 250)[,-c(1:2,13:14)]
raw.T4.N250 <- lapply(raw.T4.N250, f, 3)

raw.T5.N250 <- subset(results, N.t == 5 & N.i == 250)[,-c(1:2,13:14)]
raw.T5.N250 <- lapply(raw.T5.N250, f, 3)

raw.T4.N500 <- subset(results, N.t == 4 & N.i == 500)[,-c(1:2,13:14)]
raw.T4.N500 <- lapply(raw.T4.N500, f, 3)

raw.T5.N500 <- subset(results, N.t == 5 & N.i == 500)[,-c(1:2,13:14)]
raw.T5.N500 <- lapply(raw.T5.N500, f, 3)



#Pointwise RMSE versus pointwise minimum
vs.min.T4.N250 <- subset(vs.min, N.t == 4 & N.i == 250)[,-(1:4)]
vs.min.T4.N250 <- lapply(vs.min.T4.N250, f)

vs.min.T5.N250 <- subset(vs.min, N.t == 5 & N.i == 250)[,-(1:4)]
vs.min.T5.N250 <- lapply(vs.min.T5.N250, f)

vs.min.T4.N500 <- subset(vs.min, N.t == 4 & N.i == 500)[,-(1:4)]
vs.min.T4.N500 <- lapply(vs.min.T4.N500, f)

vs.min.T5.N500 <- subset(vs.min, N.t == 5 & N.i == 500)[,-(1:4)]
vs.min.T5.N500 <- lapply(vs.min.T5.N500, f)




#Construct Summary Tables
sample.size <- results[,13:14]

f2 <- function(x){tapply(x, sample.size, mean)}
mean.rmse <- lapply(results[,-c(1:2, 13:14)], f2)

f3 <- function(x){tapply(x, sample.size, max)}
max.rmse <- lapply(results[,-c(1:2, 13:14)], f3)

#This is really neat! Wow! How did I not know this!
#sapply(max.rmse, '[', 1, 2)
max.rmse <- as.data.frame(t(sapply(max.rmse, '['))) #By column!
names(max.rmse) <- c('T4.N250', 'T5.N250', 'T4.N500', 'T5.N500')

mean.rmse <- as.data.frame(t(sapply(mean.rmse, '['))) #By column!
names(mean.rmse) <- names(max.rmse)


c(tapply(true.min, sample.size, mean))

mean.true.min <- matrix(c(tapply(true.min, sample.size, mean)), nrow = nrow(max.rmse), ncol = ncol(max.rmse), byrow = TRUE)

mean.vs <- (mean.rmse - mean.true.min)/mean.true.min * 100


#To compare to the ``oracle'' for worst case, we need the worst-case oracle! That is, which fixed spec gives minimax risk for each sample size? 
minimax.oracle <- matrix(apply(max.rmse[1:4,], 2, min), nrow = nrow(max.rmse), ncol = ncol(max.rmse), byrow = TRUE)

max.vs <- (max.rmse - minimax.oracle)/minimax.oracle * 100

#------------Average/Worse-case Results
round(max.rmse,2)
round(mean.rmse,2)
round(max.rmse/0.5 * 100) #in %-points relative to theta
round(mean.rmse/0.5 * 100) #in %-points relative to theta
round(mean.vs) #%-points worse than pointwise oracle
round(max.vs)#%-points worse than minimax oracle

#------------#GFIC Pointwise results
vs.min.T4.N250$GFIC
vs.min.T5.N250$GFIC
vs.min.T4.N500$GFIC
vs.min.T5.N500$GFIC

library(Hmisc)


outfilename <- 'GFIC_summary_tables.tex'
 

#Assign to variable table to keep it from trying to compile
table <- latex(round(mean.vs), file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = FALSE)

#Append the remaining tables
table <- latex(round(max.vs), file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE) 

table <- latex(round(mean.rmse,3), file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE) 

table <- latex(round(max.rmse,3), file = outfilename, greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

table <- latex(raw.T4.N250, file = 'GFIC_pointwise_tables_T4_N250.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

table <- latex(raw.T5.N250, file = 'GFIC_pointwise_tables_T5_N250.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)


table <- latex(raw.T4.N500, file = 'GFIC_pointwise_tables_T4_N500.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

table <- latex(raw.T5.N500, file = 'GFIC_pointwise_tables_T5_N500.tex', greek = TRUE,  numeric.dollar = FALSE, na.blank = TRUE, landscape = FALSE, append = TRUE)

