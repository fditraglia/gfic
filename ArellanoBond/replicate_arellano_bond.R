#------------------------------------------------------------
#Filename:        replicate_arellano_bond.R
#Author:          Frank DiTraglia
#First Version:   2013-27-11
#This Version:    2013-27-11
#------------------------------------------------------------
#This script is an attempt to replicate some of the results
#from Arellano & Bond (1991): first using the plm package,
#then using my own R code and finally using my own C++ code.
#Once I'm sure everything is working correctly, I will then
#incorporate this code into my GFIC simulation example. For
#relevant details see Arellano & Bond (1991), the plm Package Vignette (Croissant & Millo), and Cameron & Trivedi (2005)
#section 22.5.3.
#------------------------------------------------------------

library(plm) #R Package plm: Linear Models for Panel Data
data(EmplUK) #UK Unemployment data from Arellano & Bond (1992)

#Example from Section 5.4 of the plm package vignette
#Replicates Table 4b of Arellano & Bond (1991)
#Employment is explained by past values of employment (two lags), current and first lag of wages and output and current value of capital

#effect = "twoways" estimates in first differences to eliminate individual effects and includes time dummies
emp.gmm <- pgmm(log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) + lag(log(output), 0:1) | lag(log(emp), 2:99), data = EmplUK, effect = "twoways", model = "twosteps")

#Estimates in the vignette agree with the paper, but standard errors are much bigger than in the paper. 
summary(emp.gmm)



