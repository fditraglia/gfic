/*-------------------------------------------------------
  # Filename:        functionsREFE.cpp  
  # Author:          Minsu Chang
  # First Version:   2015-08-25
  # Last Updated:    2017-09-20
  #--------------------------------------------------------
# C++ functions for Fixed Effect vs. Random Effect simulation experiment.
#-------------------------------------------------------*/

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat mvrnorm(int n, vec mu, mat Sigma){
  /*-------------------------------------------------------
    # Generate draws from a multivariate normal distribution
    #--------------------------------------------------------
  #  n        number of samples
  #  mu       mean vector
  #  Sigma    covariance matrix
  #--------------------------------------------------------
  # Details:
  #           This is essentially a stripped-down version
  #           of the mvrnorm function from the MASS library
  #           in R. Through the magic of Rcpp we're 
  #           transforming the *same* standard normal draws
  #           as the R version. However, since Armadillo
  #           follows a different convention from R in its
  #           definition of the eign-decomposition, the 
  #           output of this function will *not* be the
  #           same as that of its R counterpart. Since we
  #           access R's function for generating normal
  #           draws, we can set the seed from R.
  #-------------------------------------------------------*/
  RNGScope scope;
  int p = Sigma.n_cols;
  mat X = reshape(vec(rnorm(p * n)), p, n);
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Sigma);
  X = eigvec * diagmat(sqrt(eigval)) * X;
  X.each_col() += mu;
  return(X.t());
}



double sample_quantile(colvec x, double p){
  /*-------------------------------------------------------
    # Calculates a sample quantile
    #--------------------------------------------------------
  #  x        vector of data
  #  p        probability for desired quantile (e.g. 0.25
  #             gives the first quartile, 0.5 the median)
  #--------------------------------------------------------
  # Details:
  #           There are many competing definitions of
  #           sample quantiles (see Hyndman & Fan, 1996).
  #           Here we simply use the R default definition,
  #           which corresponds to Definition 7 in Hyndman
  #           & Fan. See ?quantile in R for more details.
  #-------------------------------------------------------*/
  int n = x.n_elem;
  double m = 1 - p;
  int j = floor(n * p + m);
  double g = n * p + m - j;
  colvec x_sort = sort(x);
  return((1 - g) * x_sort(j - 1) + g * x_sort(j));
}

double MSE_trim(colvec x, double truth, double trim){
  /*-------------------------------------------------------
    # Calculates trimmed mean-squared error.
    #--------------------------------------------------------
  #  x        vector of estimates
  #  truth    true value of the parameter
  #  trim     fraction of estimates to discard (half from
  #             each tail) before calculating MSE
  #-------------------------------------------------------*/
  int k = x.n_elem;
  int tail_drop = ceil(k * trim / 2);
  
  colvec x_trimmed = sort(x);
  x_trimmed = x_trimmed(span(tail_drop, k - tail_drop - 1));
  
  colvec truth_vec = truth * ones(x_trimmed.n_elem);
  colvec errors = x_trimmed - truth_vec;
  double MSE = dot(errors, errors) / errors.n_elem;
  
  return(MSE);  
}

double MAD(colvec x, double truth){
  /*-------------------------------------------------------
    # Calculates median absolute deviation.
    #--------------------------------------------------------
  #  x        vector of estimates
  #  truth    true value of the parameter
  #-------------------------------------------------------*/
  colvec truth_vec = truth * ones(x.n_rows);
  colvec abs_dev = abs(x - truth_vec);
  return(median(abs_dev)); 
}


double coverage_prob(mat conf_intervals, double truth){
  /*-------------------------------------------------------
    # Calculates the coverage probability of a matrix of
    # confidence intervals.
    #--------------------------------------------------------
  #  conf_intervals   matrix of confidence intervals in 
  #                     which each row is a CI, the 1st
  #                     column is the lower limit, and the
  #                     2nd column is the upper limit 
  #                           
  #  truth            true value of the parameter for which
  #                       the CIs were constructed
  #-------------------------------------------------------*/
  colvec truth_vec = truth * ones(conf_intervals.n_rows);
  colvec cover_lower = conv_to<colvec>
    ::from(conf_intervals.col(0) < truth_vec);
  colvec cover_upper = conv_to<colvec>
    ::from(conf_intervals.col(1) > truth_vec);
  colvec cover = cover_lower % cover_upper;
  return(sum(cover) / cover.n_elem);
}


double median_width(mat conf_intervals){
  /*-------------------------------------------------------
    # Calculates the median width of a matrix of confidence
    # intervals.
    #--------------------------------------------------------
  #  conf_intervals   matrix of confidence intervals in 
  #                     which each row is a CI, the 1st
  #                     column is the lower limit, and the
  #                     2nd column is the upper limit 
  #-------------------------------------------------------*/
  colvec width = conf_intervals.col(1) - conf_intervals.col(0);
  return(median(width));
}

double sumAll(mat NTmatrix){
  /*------------------------------------------------------
    # Calculates the sum of all the elements in NxT matrix
    #-------------------------------------------------------
  # NTmatrix   matrix consisting of N*T elements for panel data
  #------------------------------------------------------*/
  
  int T = NTmatrix.n_cols;
  int N = NTmatrix.n_rows;
  double total=0.0;
  
  for (int ii=0; ii<N; ii++){
    for (int tt=0; tt<T; tt++){
      total = total + NTmatrix(ii,tt);
    }
  }
  return(total);
}

class dgp_Panel {
  /*--------------------------------------------------
    # Class for simulating from dgp for panel data
    #-------------------------------------------------*/
    public:
    mat x;  //sims for observed regressor
  mat y;  //sims for outcome
  //Class constructor
  dgp_Panel(double, double, double, int, int);
};

dgp_Panel::dgp_Panel(double b, double rho, double gamma, 
                     int N, int T){
  /*--------------------------------------------------
    # Constructor: generates simulations from dgp
    #---------------------------------------------------
  # Arguments:
  #  b            coefficient on x 
  #  rho          persistence of x_{it}
  #  gamma        correlation between x_{it} and \alpha_i (unobserved heterogeneity)
  #  N            sample size for i
  #  T            sample size for t
  #---------------------------------------------------
  # Initializes:
  #     x, y
  #-------------------------------------------------*/
  mat epsilons = reshape(vec(rnorm(N*T)),N,T);
  epsilons = epsilons*2.5;
  
  arma::mat Mu = zeros<mat>(T+1,1);
  arma::mat CovVar = zeros<mat>(T,T);
  
  for (int ii=0; ii<T; ii++){
    for (int jj=0; jj<T; jj++){
    if(ii == jj){
    CovVar(ii,jj) = 1;
  } else {
    CovVar(ii,jj) = rho;
  }
    }
  }

arma::mat Gammas = ones<mat>(T,1);
Gammas = Gammas*gamma;

CovVar.insert_cols(T, Gammas);
arma::mat va = ones<mat>(1,1); 
Gammas.insert_rows(T, va);

CovVar.insert_rows(T, Gammas.t());
mat CSigma = CovVar;
  
mat DataDraws = mvrnorm(N, Mu, CSigma);

vec alphas = DataDraws.col(T);
DataDraws.shed_col(T);
x = DataDraws;

  arma::mat A=zeros<mat>(N,T);
  y=A;
    
  for (int ii=0; ii<N; ii++){
    for (int tt=0; tt<T; tt++){
      y(ii,tt) = b*x(ii,tt)+alphas(ii,0)+epsilons(ii,tt);
    }
  }
}

class fmsc_FE_RE {
  /*--------------------------------------------------
    # Class for FMSC calculations in FE vs RE example
    #-------------------------------------------------*/
    public:
    //class constructor
  fmsc_FE_RE(const mat&, const mat&);
  //member functions
  //    double Tfmsc();       //return FMSC "test statistic"
  double b_fe();       //return FE estimate
  colvec elements_all();
  colvec vec_test();
  mat resid_fe();
  double b_re();      //return RE estimate
  double b_fmsc();      //return FMSC-selected estimate
  double b_AVG();       //return feasible averaging estimate
  
  private:
    double fe_estimate, re_estimate, ols_estimate, v_eps_hat, v_u_hat, v_alpha_hat, var_fe, var_re, cov_fr;
  double tau_hat, v_tau_hat, K_hat, K_tilde, cc, sigma_sq, eta_sq, Tfmsc;
  int T, N;
  colvec xbar, ybar;
  mat fe_resid_sq, re_resid_sq, ols_resid_sq;
};


fmsc_FE_RE::fmsc_FE_RE(const mat& x, const mat& y){
  /*--------------------------------------------------
    # Constructor: calculates all "basic quantities"
    #---------------------------------------------------
  # Arguments:
  #  x          vector of obs. for regressor
  #  y          vector of obs. for outcome
  #---------------------------------------------------
  # Initializes:
  #    n_z, n, xx, g_sq, s_e_sq_ols, s_e_tsls, 
  #    s_x_sq, s_v_sq, tau, tau_var, ols_estimate, 
  #    ols_resid, tsls_estimate, tsls_resid, 
  #    first_stage_coefs, zx, zz_inv, zz
  #-------------------------------------------------*/
  T = y.n_cols;
  N = y.n_rows;
  
  arma::mat A=zeros<mat>(N,T);
  arma::mat B=eye(T,T);
  arma::mat C=ones<mat>(T,1);
  arma::mat D=zeros<mat>(N,1);
  
  mat Pt=C*C.t()/T;
  xbar=D;
  ybar=D;
  
  double term1=0.0;
  double term2=0.0;
  var_fe=0.0;
  var_re=0.0;
  cov_fr=0.0;
  
  for (int ii=0; ii<N; ii++){
    
    rowvec xxrow=x.row(ii);
    rowvec yyrow=y.row(ii);
    xbar(ii,0)=sum(xxrow)/T;
    ybar(ii,0)=sum(yyrow)/T;
    
    colvec xx_vfe=xxrow.t()-Pt*xxrow.t();
    colvec yy_vfe=yyrow.t()-Pt*yyrow.t();
    term1=as_scalar(xx_vfe.t()*xx_vfe)+term1;
    term2=as_scalar(xx_vfe.t()*yy_vfe)+term2;
    
  }
  
  fe_estimate = term2/term1;
  
  fe_resid_sq=A;
  
  for (int ii=0; ii<N; ii++){
    for (int tt=0; tt<T; tt++){
      fe_resid_sq(ii,tt) = pow((y(ii,tt)-ybar(ii,0))-(x(ii,tt)-xbar(ii,0))*fe_estimate,2);
    }
  }
  
  v_eps_hat= sumAll(fe_resid_sq)/(N*(T-1)-1); 
    
  colvec yvec=reshape(y, N*T, 1);
  colvec xvec=reshape(x, N*T, 1);
  
  double xx = dot(xvec, xvec);
  ols_estimate = dot(xvec, yvec) / xx;
  ols_resid_sq = reshape(pow(yvec - xvec * ols_estimate,2), N,T);
  
  v_u_hat = sumAll(ols_resid_sq)/(N*T-1); 
  
  v_alpha_hat = v_u_hat-v_eps_hat; 
  
  mat Omega_inv = 1/v_eps_hat*(B-v_alpha_hat/(T*v_alpha_hat+v_eps_hat)*C*C.t());
  
  
  double lambda = 1-sqrt(v_eps_hat/(T*v_alpha_hat+v_eps_hat));
 
  double term3=0.0;
  double term4=0.0;
  double term5=0.0;
  double term6=0.0;
  double term7=0.0;
  double term8=0.0;
  double term9=0.0;
  
  for (int ii=0; ii<N; ii++){
    
    rowvec xxrow=x.row(ii);
    rowvec yyrow=y.row(ii);
    
    term3=as_scalar(xxrow*Omega_inv*xxrow.t())+term3;
    term4=as_scalar(xxrow*Omega_inv*yyrow.t())+term4;
    term5=as_scalar(xxrow*(B-Pt)*Omega_inv*xxrow.t())+term5;
    term6=as_scalar(xxrow*Omega_inv*(yyrow.t()-xxrow.t()*fe_estimate))+term6;
    
    colvec xx_vre=xxrow.t()-lambda*Pt*xxrow.t();
    colvec yy_vre=yyrow.t()-lambda*Pt*yyrow.t();
    colvec xx_vfe=xxrow.t()-Pt*xxrow.t();
    
    term7=as_scalar(xx_vre.t()*xx_vre)+term7;
    term8=as_scalar(xx_vfe.t()*xx_vfe)+term8;
    term9=as_scalar(xx_vre.t()*yy_vre)+term9;
    
  }
  
  re_estimate = term4/term3;
  
  var_re= v_eps_hat*N/term7;
  var_fe= v_eps_hat*N/term8; 
  cov_fr = (v_eps_hat/term1)*term5*N/term3;
  
  tau_hat = (1/sqrt(N))*(T*v_alpha_hat+v_eps_hat)*term6;
  
  K_hat = term1/N;
  K_tilde = N/term3;
  
  v_tau_hat = pow(T*v_alpha_hat+v_eps_hat,2)*(var_fe/var_re-1)/K_tilde; // positive guaranteed due to var_fe>=var_re
  
  cc = K_tilde/(T*v_alpha_hat+v_eps_hat);
  sigma_sq = v_tau_hat;
  eta_sq = K_tilde;
  Tfmsc = pow(tau_hat,2)/sigma_sq; // FMSC test statistic
  
}



double fmsc_FE_RE::b_fe(){
  //Member function of class fmsc_FE_RE
  //Extracts FE estimate
  return(fe_estimate);
}

double fmsc_FE_RE::b_re(){
  //Member function of class fmsc_FE_RE
  //Extracts RE estimate
  return(re_estimate);
}

double fmsc_FE_RE::b_fmsc(){
  //Member function of class fmsc_FE_RE
  //Calculates estimator selected by FMS
  double out;
  if(Tfmsc < 2){
    out = re_estimate;
  } else {
    out = fe_estimate;
  }
  return(out);
}


double fmsc_FE_RE::b_AVG(){
  //Member function of class fmsc_FE_RE
  //Calculates feasible version of AMSE-optimal averaging estimator
  
  double nume = var_fe - cov_fr;
  double tau_sq = pow(tau_hat,2) - v_tau_hat;
  
  double sq_bias_est;
  //If the squared bias estimate is negative, set it to zero
  //so weight lies in [0,1]
  if(tau_sq >= 0){
    sq_bias_est = tau_sq;
  } else {
    sq_bias_est = 0;
  }
  
  double denom = (sq_bias_est/pow(T*v_alpha_hat+v_eps_hat,2))*pow(var_re,2) + var_re - 2*cov_fr + var_fe;
  
  double omega = nume/denom;  
  double out = omega * re_estimate + (1 - omega) * fe_estimate;
  return(out);
}


colvec fmsc_FE_RE::elements_all(){
   NumericVector vout = NumericVector::create(Tfmsc, K_tilde, K_hat, fe_estimate, re_estimate);
 
  return(vout);
}

mat fmsc_FE_RE::resid_fe(){
  return(fe_resid_sq);
}



// [[Rcpp::export]]
double PanelTry(double b, double rho, double gamma, 
                int N, int T){
  
  dgp_Panel sim_data(b, rho, gamma, N, T);
  fmsc_FE_RE sim_results(sim_data.x, sim_data.y);
  
  double out = sim_results.b_re();
  return(out);    
}

// [[Rcpp::export]]
colvec CI_elements(double b, double rho, double gamma,
                   int N, int T){
  dgp_Panel sim_data(b, rho, gamma, N, T);
  fmsc_FE_RE sim_results(sim_data.x, sim_data.y);
  
  colvec out = sim_results.elements_all();                        
  return(out);
}

// [[Rcpp::export]]
NumericVector mse_compare_cpp(double b, double rho, double gamma, 
                              int N, int T, int n_reps){
  //Function to run n_reps of the simulation study and calculate the MSE
  //of various estimators
  
  colvec fe(n_reps);
  colvec re(n_reps);
  colvec fmsc(n_reps);
  colvec AVG(n_reps);
  
  for(int ii = 0; ii < n_reps; ii++){
    
    dgp_Panel sim_data(b, rho, gamma, N, T);
    fmsc_FE_RE sim_results(sim_data.x, sim_data.y);
    
    fe(ii) = sim_results.b_fe();
    re(ii) = sim_results.b_re();
    fmsc(ii) = sim_results.b_fmsc();
    AVG(ii) = sim_results.b_AVG();
    
  }
  
  double const trim_frac = 0; //Change this if you want trimmed MSE
  
  double MSE_fe = MSE_trim(fe, b, trim_frac);
  double MSE_re = MSE_trim(re, b, trim_frac);
  double MSE_fmsc = MSE_trim(fmsc, b, trim_frac);
  double MSE_star = MSE_trim(AVG, b, trim_frac);
  
  //Create and return vector of results
  NumericVector out = NumericVector::create(MSE_fe, MSE_re, MSE_fmsc, MSE_star);
  out.names() = CharacterVector::create("FE", "RE", "FMSC", "AVG"); 
  return(out);    
}


// [[Rcpp::export]]
mat Panelmat(double b, double rho, double gamma, 
             int N, int T){
  
  dgp_Panel sim_data(b, rho, gamma, N, T);
  fmsc_FE_RE sim_results(sim_data.x, sim_data.y);
  mat outMat = sim_results.resid_fe();
  
  return(outMat);    
}



// [[Rcpp::export]]
mat CI_elements_mat(double b , double rho, double gamma,
                    int N, int T){                        
  int n_reps = 1;
  
  mat Output(4, n_reps, fill::zeros);
  
  for(int ii = 0; ii < n_reps; ii++){
    
    //Simulate data and calculate estimators
    dgp_Panel data(b, rho, gamma, N, T);
    fmsc_FE_RE results(data.x, data.y);
    
    //CIs that do not require draw_CI_sims()
    Output.col(ii) = results.elements_all();
  }
  
  return(Output);
}


