# Tables of RMSE results for Dynamic Panel Simulation
setwd("~/gfic/simulations/DynamicPanel/results")

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
names(results)[1:2] <- c('g', 'r')
results$oracle <- apply(results[, c('LW', 'LS', 'W', 'S')], 1, min)

rm(T4_N250, T4_N500, T5_N250, T5_N500)

# Coarsen the simulation grid over g and r
g_coarse <- c(0, 0.05, 0.1, 0.15)
r_coarse <- c(0, 0.05, 0.1, 0.15)
coarse <- subset(results, (g %in% g_coarse) & (r %in% r_coarse))

# Use xtabs and ftable to make "flat" tables of RMSE in THOUSANDTHS
# i.e. 1000 * RMSE
make_R_table <- function(col) {
  tab <- xtabs(1000 * round(col, 3) ~ g + r + ni + nt, coarse)
  names(dimnames(tab)) <- c('$\\gamma$', '$\\rho$', '$N$', '$T$')
  tab <- ftable(tab, row.vars = c('$T$', '$N$', '$\\gamma$'), 
                col.vars = c('$\\rho$'))
  tab <- format(tab, justify = 'none', trim = FALSE)
  tab <- apply(tab, 1:2, gsub, pattern = '\"', replacement = '')
  tab <- apply(tab, 1:2, gsub, pattern = ' ', replacement = '')
  return(tab)
}

#--------------------- Make ftables of results for each estimator

LP <- make_R_table(coarse$LW) # Notice the new name for these moment conditions
LS <- make_R_table(coarse$LS)
P <- make_R_table(coarse$W) # Notice the new name for these moment conditions
S <- make_R_table(coarse$S)
oracle <- make_R_table(coarse$oracle)
GFIC <- make_R_table(coarse$GFIC)
J10 <- make_R_table(coarse$J10)
J5 <- make_R_table(coarse$J5)
AIC <- make_R_table(coarse$AIC)
BIC <- make_R_table(coarse$BIC)
HQ <- make_R_table(coarse$HQ)

# -------------------- Convenience functions for converting to TeX
make_tex_row <- function(char_vec) {
  out <- paste0(char_vec, collapse = ' & ')
  out <- paste(out, '\\\\')
  return(out)
}

make_tex_tabular <- function(char_mat) {
  rows_vec <- apply(char_mat, 1, make_tex_row)
  paste0(rows_vec, collapse = '\n')
}

make_table_header <- function(names_vec) {
  first_part <- paste0('\\begin{tabular}{', 
                       paste0(rep('cccc', 6), collapse = '|'),  
                       '} \n \\hline \\hline \n')
  second_part <- paste0('\\multicolumn{4}{c}{', c('', names_vec), '}',
                        collapse = '&')
  third_part <-  '\\\\ \n \\hline'
  return(paste0(first_part, second_part, third_part))
}

make_tex_table <- function(header_names, body_mat){
  table_header <- make_table_header(header_names)
  table_body <- make_tex_tabular(body_mat)
  paste0(c(table_header, table_body, '\\hline', 
           '\\end{tabular}'), collapse = '\n')
}

#---------------------------- Assemble ftables into three nice LaTeX tables
table1 <- make_tex_table(c('GFIC', 'LP', 'LS', 'P', 'S'),
                         cbind(GFIC, LP[,-c(1:4)], LS[,-c(1:4)], P[,-c(1:4)],
                               S[,-c(1:4)]))

table2 <- make_tex_table(c('Oracle', 'GFIC', 'J-test 5\\%', 'GMM-BIC', 'GMM-AIC'),
                         cbind(oracle, GFIC[,-c(1:4)], J5[,-c(1:4)], 
                               BIC[,-c(1:4)], AIC[,-c(1:4)]))

table3 <- make_tex_table(c('Oracle', 'GFIC', 'J-test 10\\%', 'GMM-HQ'),
                         cbind(oracle, GFIC[,-c(1:4)], J10[,-c(1:4)], 
                               HQ[,-c(1:4)]))

#-------------------------- Clean up
cat(table1, file = 'Dpanel_RMSE_GFIC_vs_fixed_spec.tex')
cat(table2, file = 'Dpanel_RMSE_GFIC_vs_alternatives_main.tex')
cat(table3, file = 'Dpanel_RMSE_GFIC_vs_alternatives_append.tex')
rm(list = ls())