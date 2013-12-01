/*------------------------------------------------------------
Filename:        ABfit.cpp
Author:          Frank DiTraglia
First Version:   2013-30-11
This Version:    2013-30-11
--------------------------------------------------------------
This function calculates the Arellano-Bond estimator in a 
simple setting: one lag of y and a single exogenous regressor.
This corresponds to the simulation study in Table 1 of
Arellano & Bond (1991).
------------------------------------------------------------*/
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ABfit_cpp(NumericMatrix x_r, NumericMatrix y_r) {
   
   int N_i = x_r.nrow();
   int N_t = x_r.ncol();
   
   //Number of regressors in this example
   int N_reg = 2;
   
   //Initialize armadillo matrices for R input matrices
   //and reuse original memory
   arma::mat x(x_r.begin(), N_i, N_t, false);
   arma::mat y(y_r.begin(), N_i, N_t, false);
   
   //First-differences: resulting columns correspond to t = 2,...,T
   arma::mat xdiff = x.cols(1, N_t - 1) - x.cols(0, N_t - 2);
   arma::mat ydiff = y.cols(1, N_t - 1) - y.cols(0, N_t - 2);
   
   //Initialize Covariance matrix of differenced errors
   //with 2 on main diagonal, -1 on the sub and super
   //diagonals
   arma::mat H = 2 * arma::mat(N_t - 2, N_t - 2, arma::fill::eye);
   H.diag(-1) = -1 * arma::vec(N_t - 3, arma::fill::ones);
   H.diag(1) = -1 * arma::vec(N_t - 3, arma::fill::ones);

   //Number of instruments in each time period
   arma::vec N_y_instruments = arma::linspace(1, N_t - 2, N_t - 2);
   arma::vec N_x_instruments = arma::ones(N_t - 2);
   arma::vec N_instruments = N_y_instruments + N_x_instruments;
   
   //Total number of moment conditions
   int m = arma::sum(N_instruments);
   
   //Starting and ending columns for embedding instruments in Z_i
   //Remember: zero indexing!
   arma::vec end_cols = arma::cumsum(N_instruments) - 1;
   arma::vec start_cols = end_cols - N_instruments + 1;
   
   //Initialize instrument matrix Z_i and fill with zeros
   arma::mat Z_i(N_t - 2, m, arma::fill::zeros);
   
   //Initialize all other matrices for individual i that
   //appear in the loop over i and fill with zeros
   arma::colvec ydiff_lag_i(N_t - 2, arma::fill::zeros);
   arma::colvec xdiff_i(N_t - 2, arma::fill::zeros);
   arma::mat X_tilde_i(N_t - 2, N_reg, arma::fill::zeros);
   arma::colvec y_tilde_i(N_t - 2, 1, arma::fill::zeros);
   arma::mat XZ_i(N_reg, m, arma::fill::zeros);
   arma::colvec Zy_i(m, 1, arma::fill::zeros);
   arma::mat ZHZ_i(m, m, arma::fill::zeros);
   
   //Initialize Matrices that accumulate over individuals in i loop
   arma::mat XZ = XZ_i;
   arma::mat ZHZ = ZHZ_i;
   arma::colvec Zy = Zy_i;
   
   //Loop over individuals
   for(int i = 0; i < N_i; i++){
     
     //Loop over time periods to construct matrix products
     //For a given individual: j indexes rows of Z_i as well as
     //elements of end_cols and start_cols
     for(int j = 0; j < N_t - 2; j++){
       
       //Careful: y starts at t = 1, xdiff starts at t = 2
       Z_i(arma::span(j,j), arma::span(start_cols(j), end_cols(j))) = 
           arma::join_rows(  y( arma::span(i, i), arma::span(0,j) ),  
           xdiff( arma::span(i, i), arma::span(j + 1, j + 1) ));
     }
     
     //Note: both ydiff and xdiff start at t = 2 but ydiff enters
     //as a lagged RHS variable.
     ydiff_lag_i = ydiff( arma::span(i, i), 
                          arma::span(0, N_t - 3) ).t();
     xdiff_i = xdiff( arma::span(i, i), 
                      arma::span(1, N_t - 2) ).t();
     X_tilde_i = arma::join_rows(ydiff_lag_i, xdiff_i);
   
     y_tilde_i = ydiff( arma::span(i, i), 
                        arma::span(1, N_t - 2) ).t();
     
     //Construct XZ, ZHZ, and Zy for person i 
     XZ_i = X_tilde_i.t() * Z_i;
     Zy_i = Z_i.t() * y_tilde_i;
     ZHZ_i = Z_i.t() * H * Z_i;     
   
     //Accumulate over i
     XZ += XZ_i;
     Zy += Zy_i;
     ZHZ += ZHZ_i;
     
   }
   

   //Note that chol(M) returns upper-triangular matrix R
   //M = L * L.t() = R.t() * R
   arma::mat L = chol(ZHZ / N_i).t();
   arma::mat Xw = solve(trimatl(L), XZ.t());
   arma::colvec yw = solve(trimatl(L), Zy);
   
   //Now we simply have a least squares problem in Xw, yw
   arma::colvec b = solve(Xw, yw);
   
   arma::mat W_sqrt = chol(ZHZ / N_i);
   
   return List::create(Named("b") = b,
                       Named("W_inv") = ZHZ/N_i,
                       Named("XZ") = XZ,
                       Named("Zy") = Zy);
   
}
