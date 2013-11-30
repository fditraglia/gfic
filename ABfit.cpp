//C++ version of ABfit

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ABfit_cpp(NumericMatrix x_r, NumericMatrix y_r) {
   
   int N_i = x_r.nrow();
   int N_t = x_r.ncol();
   
   //Initialize armadillo matrices for R input matrices
   //and reuse original memory
   arma::mat x(x_r.begin(), N_i, N_t, false);
   arma::mat y(y_r.begin(), N_i, N_t, false);
   
   //First-differences: resulting columns correspond to t = 2,...,T
   arma::mat xdiff = x.cols(1, N_t - 1) - x.cols(0, N_t - 2);
   arma::mat ydiff = y.cols(1, N_t - 1) - y.cols(0, N_t - 2);
   
   //Initialize Covariance matrix of differenced errors
   //with 2 on main diagonal, -1 on the sib and super
   //diagonals
   arma::mat H = 2 * arma::mat(N_t - 2, N_t - 2, arma::fill::eye);
   H.diag(-1) = -1 * arma::vec(N_t - 3, arma::fill::ones);
   H.diag(1) = -1 * arma::vec(N_t - 3, arma::fill::ones);
   
   //Construct XZ, Zy for each i
   
   return List::create(Named("xdiff") = xdiff, 
        Named("ydiff") = ydiff, Named("H") = H);
    
  
//  #Instruments for each individual (list of lists)
//  Z.i <- function(i){
//    
//    lapply(3:N.t, function(j) c(y[[i]][1:(j - 2)], x.diff[[i]][j]))
//    
//  }
//  
//  Z <- lapply(individuals, Z.i)
//  
//  #Convert to list of sparse matrices
//  Z <- lapply(individuals, function(i) t(bdiag(Z[[i]])))
//  
//  
//  #Regressors for each individual (list of matrices)
//  X.tilde.i <- function(i){
//    
//    cbind(y.diff[[i]][2:(N.t - 1)], x.diff[[i]][3:N.t])
//    
//  }
//  
//  X.tilde <- lapply(individuals, X.tilde.i)
//  
//  #Outcomes for each individual (list of vectors)
//  y.tilde.i <- function(i){
//    
//    y.diff[[i]][3:N.t]
//    
//  }
//  y.tilde <- lapply(individuals, y.tilde.i)
//  
//  
//  XZ <- lapply(individuals, function(i) crossprod(X.tilde[[i]], Z[[i]]))
//  XZ <- Reduce('+', XZ) #Sum over all individuals
//  
//  Zy <- lapply(individuals, function(i) crossprod(Z[[i]], y.tilde[[i]]))
//  Zy <- Reduce('+', Zy)
//  
//  
//  #Probably better not to use "solve" here but I'm not sure how to handle a qr decomposition for sparse matrices...
//  H <- bandSparse(5, 5, k = c(0,-1, 1), list(rep(2, 5), rep(-1, 5), rep(-1, 5)))
//  W.inv <- lapply(individuals, function(i) t(Z[[i]]) %*% H %*% Z[[i]])
//  W.inv <- Reduce('+', W.inv) / N.i
//  W <- solve(W.inv)
//  
//  XZW <- XZ %*% W
//  K.inv <- XZW %*% t(XZ)
//  K <- chol2inv(qr.R(qr(K.inv)))
//  b <- solve(K.inv) %*% XZW %*% Zy
//  #Direct calculation is: solve(XZ %*% W %*% t(XZ)) %*% XZ %*% W %*% Zy
//  b <- as.vector(b)
//  names(b) <- c("a", "b")
//  return(b)
   
   
}
