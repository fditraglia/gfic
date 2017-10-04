gfic
====
This repository contains all .tex files and source code for the following paper:

* Minsu Chang & Francis J. DiTraglia, "A Generalized Focused Information Criterion for GMM", *Journal of Applied Econometrics*, forthcoming. 

The file "main.pdf" contains a copy of the paper itself, including the online-only appendix at the end of the document.

Code for our simulation studies appears in the directory "simulations." 
The subdirectories "DynamicPanel" and "REvsFE" give code for the dynamic panel example and random versus fixed effects example, respectively.
(Note that the random versus fixed effects example appears in the online only appendix to the paper.)

Our empirical example uses the dataset of Baltagi and Levin (1992) and Griffin and Xiong (2000), as provided by the "plm" package in R. 
The dataset is a panel of 46 U.S. states over the period 1963-1992 containing the following variables:

* (1) STATE = State abbreviation.
* (2) YR = YEAR.
* (3) Price per pack of cigarettes.
* (4) Population.
* (5) Population above the age of 16.
* (6) CPI = Consumer price index with (1983=100).
* (7) NDI = Per capita disposable income.
* (8) C = Cigarette sales in packs per capita.
* (9) PIMIN = Minimum price in adjoining states per pack of cigarettes. 

All code for our empirical example is contained in the directory "empirical_example." 
Our code loads the dataset directly from the "plm" package.
In case this dataset is excluded from future versions of teh package, we include a .csv file "Cigar.csv" in the "empirical_example" directory.

