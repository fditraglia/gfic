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
List dgp_cpp(double a1, double a2, double g, int N_i, 
                      int N_t, int burn_in, double b, 
                      double r, double theta, double s_e, 
                      double s_eta, double s_v){
    
    RNGScope scope;
    
    //Generate exogenous errors
    arma::colvec epsilon_sims = rnorm(N_i * (N_t + burn_in), 0 , s_e);
    arma::mat epsilon = reshape(epsilon_sims, N_i, N_t + burn_in);
    arma::colvec v_sims = rnorm(N_i * (N_t + burn_in), 0, s_v);
    arma::mat v = reshape(v_sims, N_i, N_t + burn_in);
    arma::colvec eta = rnorm(N_i, 0, s_eta);

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
    for(int i = 2; i < (N_t + burn_in); i++){
  
      xi.col(i) = r * xi.col(i - 1) + epsilon.col(i);
      x.col(i) = theta * eta + g * v.col(i - 1) + xi.col(i);
      y.col(i) = a1 * y.col(i - 1) + a2 * y.col(i - 2) + 
                    b * x.col(i) + eta + v.col(i);
      
    }
    
    //Discard Burn-in Samples
    x = x.cols(burn_in, burn_in + N_t - 1);
    y = y.cols(burn_in, burn_in + N_t - 1);
    
    return List::create(Named("x") = x, Named("y") = y);
                  
}