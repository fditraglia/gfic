/*------------------------------------------------------------
Filename:        dgp.cpp
Author:          Frank DiTraglia
First Version:   2013-26-11
This Version:    2013-26-11
--------------------------------------------------------------
Generates simulated data for the Arellano-Bond example from my
GFIC paper.
------------------------------------------------------------*/
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List dgp_cpp(double a1, double a2, int g, int N_i, 
                      int N_t, int burn_in, double b, 
                      double r, double theta, double s_e, 
                      double s_eta, double s_v){
          
    //Generate exogenous errors
    arma::mat epsilon = s_e * arma::randn(N_i, N_t + burn_in);
    arma::mat v = s_v * arma::randn(N_i, N_t + burn_in);
    arma::colvec eta = s_eta * arma::randn(N_i);

    //Initialize matrices to store x, y, xi 
    arma::mat x(N_i, N_t + burn_in);
    arma::mat y(N_i, N_t + burn_in);
    arma::mat xi(N_i, N_t + burn_in);
    
    //First Initialization Step - Presample Obs = 0
    xi.col(0) = epsilon.col(0); //Remember: Arma zero-indexes!
    x.col(0) = theta * eta + xi.col(0); 
    y.col(0) = b * x.col(0) + eta + v.col(0);
  
    //Second Initialization Step - Presample Obs = 0
    xi.col(1) = r * xi.col(0) + epsilon.col(1);
    x.col(1) = theta * eta + g * v.col(0) + xi.col(1);
    y.col(1) = a1 * y.col(0) + b * x.col(1) + eta + v.col(1);
    
    //Generate Endogenous Variables by looping over columns
    for(int col = 2; col < (N_t + burn_in); col++){
  
      xi.col(col) = r * xi.col(col - 1) + epsilon.col(col);
      x.col(col) = theta * eta + g * v.col(col - 1) + xi.col(col);
      y.col(col) = a1 * y.col(col - 1) + a2 * y.col(col - 2) + 
                    b * x.col(col) + eta + v.col(col);
      
    }
    
    //Discard Burn-in Samples
    x = x.cols(burn_in, burn_in + N_t - 1);
    y = y.cols(burn_in, burn_in + N_t - 1);
    
    return List::create(Named("x") = x, Named("y") = y);
                        
}