gfic
====
This repository contains all .tex files and source code for the following paper:

* Minsu Chang & Francis J. DiTraglia, "A Generalized Focused Information Criterion for GMM", *Journal of Applied Econometrics*, forthcoming. 

The file "main.pdf" contains a copy of the paper itself, including the online-only appendix at the end of the document.


Required Software and Packages
------------------------------
Calculations for the simulations and empirical examples were carried out in a mixture of `R` and `C++` using the package `Rcpp` as an interface to `C++` from `R`. A full list of the packages needed to replicate our results is as follows:

  * `MASS`
  * `Matrix`
  * `Rcpp`
  * `RcppArmadillo`
  * `parallel`
  * `ggplot`
  * `RColorBrewer`
  * `tikzDevice`
  * `plm`

Note that `plm` is used only to load the data for our empirical example (see below for details) and `ggplot`, `RcolorBrewer` and `tikzDevice` are used only to generate the Figures (B1-B10) that appear in the online-only appendix. The package `parallel` is used only the generate the Random versus Fixed Effects results from the online-only appendix (Figures B1 and B2). For more information, see the comments beginning in line 15 of `simulations/REvsFE/RUN_ME.R`.


Simulation Studies
-------------------
Code for our simulation studies appears in the directory `simulations` which
contains three subdirectories:

  * `DynamicPanel` contains code to replicate Table 2 from the body of the paper 
  and figures B3 through B10 from the online Appendix.
  * `SRvsLR` contains code to replicate Table 1 from the body of the paper.
  * `REvsFE` contains code to replicate Figures B1 and B2 from the online
  Appendix.
  
To replicate any of these examples, navigate to the corresponding directory and
create a subdirectory called `results` if it does not already exist. Then run 
the script `RUN_ME.R`. 

Empirical Example
-----------------

Our empirical example uses the dataset of Baltagi and Levin (1992) and Griffin and Xiong (2000), as provided by the `plm` package in R. 
The dataset is a panel of 1309 observations for 46 U.S. states over the period 1963-1992 containing the following variables:

* (1) `state` = State abbreviation.
* (2) `year` = Year.
* (3) `price` = Price per pack of cigarettes.
* (4) `pop` = Population.
* (5) `pop16` = Population above the age of 16.
* (6) `cpi` = Consumer price index with (1983=100).
* (7) `ndi` = Per capita disposable income.
* (8) `sales` = Cigarette sales in packs per capita.
* (9) `pimin` = Minimum price in adjoining states per pack of cigarettes. 

All code for our empirical example is contained in the directory `empirical_example`. 
Our code loads the dataset directly from the "plm" package.
In case this dataset is excluded from future versions of the package, we include a .csv file `Cigar.csv` in the `empirical_example` directory.
To replicate the empirical example, navigate to the appropriate directory and run the script `RUN_ME.R`
