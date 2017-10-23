Minsu Chang and Francis J. DiTraglia, "A Generalized Focused Information Criterion for GMM 
with Applications to Panel Data Models", Journal of Applied Econometrics, ...

There are three folders, each with a set of programs and data. 

1. Section6: 

This folder contains two folders including the R codes to replicate the results in Section 6.

1) Folder "SRLR" is for replicating Table 1: For long-run versus short-run effects, run the code RunMe_SRLR.R. 
This code uses the functions in GFIC_simulation_SRLR.R which computes the mean absolute deviation (MAD) for 
estimators of the short-run and long-run effects. Before running the codes, change the working directory 
in RunMe_SRLR.R and GFIC_simulation_SRLR.R.


2) Folder "DynamicPanel" is for replicating Table 2: For model/moment selection for the short run effect, 
run the code rmse_tables.R. 
(Explanation...? I am not sure which files can be removed. For example, funtions/prelim folder can be removed?)
 

2. Section7: 

This folder includes the R codes and the data to replicate Table 3 in section 7. 
Cigar.csv is the data file used in the paper Baltagi, Griffin and Xiong (2000). To replicate Table 3-(b), run the code RunMe_empirical.R. 
To replicate Table 3-(a), you can change the time period of interest to be 1975-1980.

3. Appendix: 

This folder includes the R codes to replicate the results in section B of the online appendix. This section is about
the fixed vs. random effects example. 

To replicate the Figure B.1., run the code RunMe_REFE.R. 
This code is based on mse_calculationsREFE.R (which computes RMSE results for simulation), 
rmse_plotsREFE.R (which plots the results and export in tikz format to Results folder. You need to create a folder 
named Results so that the outputs are saved there.), and functionsREFE.cpp (which includes all the relevant functions).