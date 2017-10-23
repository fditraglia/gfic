#-----------------------------------------------------------------#
#SCRIPT:  rmse_simulation.R
#AUTHOR:  Frank DiTraglia
#EMAIL:   fditraglia@gmail.com 
#DATE:    June 4th, 2012
#-----------------------------------------------------------------#
#DESCRIPTION:   
#
#     This script defines functions to carry out a simulation studying the effect of mis-specification on the MSE of a target parameter in a linear, dynamic panel model. This script assumes that the following scripts have been loaded into memory:
#
# Dpanel_sim.R              Function to generate a simulated panel dataset
# reg2SLS.R                 Function to get 2SLS estimate
# instrument_matrices.R     Functions to set up instrument matrices
# model_matrices.R          Functions to set up ``X'' and ``y''
#
#     The correct model is:
#
#     y[i,t] = gamma * y[i,t-1] + theta * x[i,t] + eta[i] + v[i,t]
#
#The random variables x[i,t], eta[i] and v[i,t] are normal with mean zero and constant variances: s.x.sq, s.eta.sq, s.v.sq. The covariance between x[i,t] and eta[i] is s.x.eta; the covariance between x[i,t] and v[i,t-1] is s.x.v while that between x[i,t] and v[i,k] for (k != t-1) is zero. That is, x is weakly exogenous (predetermined) but is correlated with the one-period lagged error.
#     The purpose of this simulation is to calculate the RMSE of four estimators of the paramters for varying values of gamma and s.x.v. The four estimators are as follows:
#
#     Name        Lagged y?     Assume x Strictly Exogenous?
#     --------    ---------     ----------------------------
#     b.LW        YES           NO (assume weakly exogenous)
#     b.LS        YES           YES
#     b.W         NO            NO (assume weakly exogenous)
#     b.S         NO            YES
#
#The first two estimators, b.LW and b.LS, use the correct model and hence estimate (gamma, theta). The second two estimators omit the lagged y, estimating only theta.
#     The idea here is to compare bias versus variance. By omitting the lagged effect, we gain a time period for estimation and have one less paramter to estimate; by assuming strict exogeneity, we gain additional instruments. Both of these decrease the variance of the estimator. However, if the assumptions are wrong they introduce bias. The question is, how small must gamma be for it to make sense to omit the dynamic effect from an RMSE perspective? Similarly, how small must s.x.v be for assuming strict exogneity to make sense?
#     To answer this question, this script compares the RMSE of the 2SLS estimator for b.LW, b.LS, b.W and b.S. The simulation setup described above is similar to that of Andrews and Lu (2001) with two main differences. First, I do not include a constant term. Second, I do not create extra pre-sample observations to accomodate the lagged y. That is, in my simulation omitting the lag gives an extra time period. This is important because in a short panel additional time periods could considerably reduce the variance of an estimator. Third, I set the initial value of y, y.0, equal to the mean of its stationary distribution (zero) as in Arellano and Bond (1991), unlike Andrews and Lu. Fourth, I use only strong instruments rather than all instruments available under the assumptions, essentially the standard Arellano-Bond Moment Conditions.
#-----------------------------------------------------------------#





#-----------------------------------------------------------------#
# FUNCTION FOR SINGLE SIMULATION REP AT FIXED PARAMETERS
#-----------------------------------------------------------------#
sim.rep <- function(gamma, r.x.v, N.t = 3, N.i  = 250)
{  
  #Generate a simulation from the true model
  sim <- Dpanel.sim(N = N.i, Time = N.t, a = gamma, s.x.v = r.x.v)#Since both x and v default to unit variance, the correlation r.x.v equals the covariance
  x.panel <- sim$x
  y.panel <- sim$y
  
  #Construct instrument matrices
  Z.LS <- instruments.LS(x.panel, y.panel)
  Z.LW <- instruments.LW(x.panel, y.panel)
  Z.S <- instruments.S(x.panel)
  Z.W <- instruments.W(x.panel)
  
  #Construct model matrices ``X'' and vectors ``y'' for the model with no lag
  model <- model.matrices(x.panel, y.panel)
  X.tilde <- model$X.tilde
  y.tilde <- model$y.tilde
  
  #Construct model matrices ``X'' and vectors ``y'' for the model with a lag
  model.L <- model.matrices.L(x.panel, y.panel)
  X.tilde.L <- model.L$X.tilde.L
  y.tilde.L <- model.L$y.tilde.L
  
  #Carry out 2SLS for each combination of model and instrument set
  reg.LW <- reg2SLS(X = X.tilde.L, y = y.tilde.L, Z = Z.LW)
  reg.LS <- reg2SLS(X = X.tilde.L, y = y.tilde.L, Z = Z.LS)
  reg.W <- reg2SLS(X = X.tilde, y = y.tilde, Z = Z.W)
  reg.S <- reg2SLS(X = X.tilde, y = y.tilde, Z = Z.S)
  
  #Extract the target parameter theta for each
  theta.LW <- reg.LW$b[2]
  theta.LS <- reg.LS$b[2]
  theta.W <- reg.W$b
  theta.S <- reg.S$b
  
  theta <- cbind(theta.LW, theta.LS, theta.W, theta.S)
  return(theta)
    
}#END sim.rep
#-----------------------------------------------------------------#
# END FUNCTION FOR SINGLE SIMULATION REP AT FIXED PARAMETERS
#-----------------------------------------------------------------#




#-----------------------------------------------------------------#
# FUNCTION TO CALCULATE RMSE OF THETA AT GIVEN PARAMETER VALUES
#-----------------------------------------------------------------#
RMSE.sim <- function(gamma, r.x.v, b = 0.5, N.sims = 1000, N.t = 3, N.i = 250){

  
  sim.reps <- replicate(N.sims, sim.rep(gamma, r.x.v, N.t, N.i), simplify = 'matrix')
  
  bias <- apply(sim.reps - b, 1, mean)
  variance <- apply(sim.reps, 1, var)
  RMSE <- sqrt(bias^2 + variance)
  
  names(RMSE) <- c('LW', 'LS', 'W', 'S')
  return(RMSE)
  
}#END RMSE.sim
#-----------------------------------------------------------------#
# END FUNCTION TO CALCULATE RMSE OF THETA AT GIVEN PARAMETER VALUES
#-----------------------------------------------------------------#




#-----------------------------------------------------------------#
# FUNCTION TO CALCULATE RMSE OF THETA OVER A GRID OF PARAMETERS
#-----------------------------------------------------------------#
RMSE.grid <- function(g, r, N.sims, N.i = 250, N.t = 3)
{
  
  g.vector <- rep(g, each = length(r))
  r.vector <- rep(r, times = length(g))
  
  result <- mapply(RMSE.sim, gamma = g.vector, r.x.v = r.vector, N.sims = N.sims, N.t = N.t, N.i = N.i, SIMPLIFY = 'matrix')
  
  parameters <- cbind(g.vector, r.vector)
  colnames(parameters) <- c('gamma', 'r.x.v')
  cbind(parameters, t(result))
  
  
}#END RMSE.grid
#-----------------------------------------------------------------#
# END FUNCTION TO CALCULATE RMSE OF THETA OVER A GRID
#-----------------------------------------------------------------#



#-----------------------------------------------------------------#
# END SCRIPT
#-----------------------------------------------------------------#
