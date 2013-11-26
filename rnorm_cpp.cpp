// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rnorm_cpp(double s, int N){
   
    arma::colvec epsilon = s * rnorm(N); 
    return List::create(Named("e") = epsilon);
  
}
  
  
  
