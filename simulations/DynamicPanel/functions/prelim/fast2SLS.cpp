// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List fast2SLS(NumericVector yr, NumericMatrix Xr, 
    NumericMatrix Zr, bool center) {

   int n = Xr.nrow(), k = Xr.ncol(), p = Zr.ncol();
   
   //Initialize armadillo matrices corresponding to R input matrices
   arma::mat X(Xr.begin(), n, k, false); 
   arma::mat Z(Zr.begin(), n, p, false); 
   arma::colvec y(yr.begin(), yr.size(), false);
   
   //Calculate the 2SLS estimator and residuals
   arma::mat first_stage = arma::solve(Z,X);
   arma::colvec coef = arma::solve(Z * first_stage, y); 
   arma::colvec resid = y - X * coef; 
   
   //Calculate the matrix K needed for estimated variance matrix
   arma::mat Qz, Rz;
   arma::qr_econ(Qz, Rz, Z);
   arma::mat Rzinv = arma::inv(arma::trimatu(Rz));
   arma::mat ZZinv = Rzinv * arma::trans(Rzinv);
   arma::mat Xstar = arma::trans(Qz) * X;
   arma::mat Qstar, Rstar;
   arma::qr_econ(Qstar, Rstar, Xstar);
   arma::mat Rstarinv = arma::inv(arma::trimatu(Rstar));
   arma::mat XXstarinv = Rstarinv * arma::trans(Rstarinv);
   arma::mat K = n * XXstarinv * arma::trans(X) * Z * ZZinv;
   
   //Estimate variance matrix of moment conditions
   arma::mat ZU = arma::trans(Z) * arma::diagmat(resid);
   arma::mat Omega = ZU * arma::trans(ZU)/n;
   
   //Center estimated variance matrix of moment conditions 
   if(center == true){
     
     arma::mat Zu = arma::trans(Z) * resid;
     Omega = Omega - (Zu/n * arma::trans(Zu)/n);
     
   }//END if
   
   //Estimated covariance matrix for 2SLS estimator
   arma::mat V = K * Omega * arma::trans(K)/n;
   
   //Calculate 2-step GMM estimator needed for J-test
   arma::mat Ro = arma::chol(Omega); 
   arma::mat Roinv = arma::inv(arma::trimatu(Ro));
   arma::mat W = Z * Roinv;
   arma::colvec GMMcoef = arma::solve(arma::trans(W)*X, arma::trans(W)*y);
   
   //J-test Statistic
   arma::colvec GMMresid = y - X * GMMcoef; 
   arma::colvec Wu = arma::trans(W) * GMMresid;
   double J = arma::dot(Wu, Wu)/n;
   
   return List::create(Named("b") = coef, Named("K") = K, 
                        Named("Omega") = Omega, Named("V") = V, 
                        Named("bGMM") = GMMcoef, 
                        Named("Jtest") = J);

}