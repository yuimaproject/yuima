// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
arma::cube calc_filter_vcov(arma::cube un_dr_sl, arma::cube un_diff, arma::cube ob_dr_sl, arma::cube inv_sq_ob_diff, arma::mat init, double delta) {
    int d_un = un_dr_sl.n_rows; // the number of observed variables
    int n = un_diff.n_slices;  // the number of observations
    
    // initialize vcov with suitable size, no value
    arma::cube vcov = arma::cube(d_un, d_un, n, arma::fill::none);
    vcov.slice(0) = init;

    for(int i = 1; i < n; i++){
        arma::mat vcov_prev = vcov.slice(i - 1);
        arma::mat vcov_next = vcov_prev + delta * (
            un_diff.slice(i-1) * un_diff.slice(i-1).t()
             + un_dr_sl.slice(i-1) * vcov_prev
             + vcov_prev * un_dr_sl.slice(i-1).t()
             - vcov_prev * ob_dr_sl.slice(i-1).t() * inv_sq_ob_diff.slice(i-1) * ob_dr_sl.slice(i-1) * vcov_prev);
        vcov.slice(i) = vcov_next;
    }

    //return vcov and inv_sq_ob_diff;
    return vcov;
}

// [[Rcpp::export]]
arma::cube calc_filter_vcov_time_homogeneous(arma::mat un_dr_sl, arma::mat un_diff, arma::mat ob_dr_sl, arma::mat inv_sq_ob_diff, arma::mat init, double delta, int n) {
  int d_un = un_dr_sl.n_rows; // the number of observed variables

    // initialize vcov with suitable size, no value
  arma::cube vcov(d_un, d_un, n, arma::fill::none);
  vcov.slice(0) = init;
  
  arma::mat intercept = delta * un_diff * un_diff.t();
  arma::mat vcov_vcov_coef = delta * ob_dr_sl.t() * inv_sq_ob_diff * ob_dr_sl;
  arma::mat vcov_coef = delta * un_dr_sl;
  for(int i = 1; i < n; i++){
    arma::mat vcov_prev(&vcov(0, 0, i-1), d_un, d_un, false, true);
    vcov.slice(i) = vcov_prev
                  + vcov_coef * vcov_prev + vcov_prev * vcov_coef.t()
                  - vcov_prev * vcov_vcov_coef * vcov_prev + intercept;
  }
  
  return vcov;
}

// [[Rcpp::export]]
arma::mat calc_filter_vcov_are(arma::mat un_dr_sl, arma::mat un_diff, arma::mat ob_dr_sl, arma::mat inv_sq_ob_diff) {
  arma::mat H(un_dr_sl.n_cols + un_dr_sl.n_rows, un_dr_sl.n_cols + un_dr_sl.n_rows);
  /*
   H = 
   \left(
   \begin{array}{cc}
   a(\theta_2)^\top & c(\theta_2)^\top\{\sigma(\theta_1)\sigma(\theta_1)^\top\}^{-1}c(\theta_2) \\
   b(\theta_2)b(\theta_2)^\top & -a(\theta_2) \                                                  \
   \end{array}
   \right), 
   a = \mathrm{un\_dr\_sl}, 
   b = \mathrm{un\_diff}, 
   c = \mathrm{ob\_dr\_sl}, 
   \sigma = \mathrm{ob\_diff}
   */
  for(unsigned int i = 0; i < un_dr_sl.n_rows; i++) {
    for(unsigned int j = 0; j < un_dr_sl.n_cols; j++) {
      H(j,i) = -un_dr_sl(i,j);
      H(un_dr_sl.n_cols + i, un_dr_sl.n_rows + j) = un_dr_sl(i,j);
    }
  }
  arma::mat lower_left_matrix = un_diff * un_diff.t();
  for(unsigned int i = 0; i < un_dr_sl.n_rows; i++) {
    for(unsigned int j = 0; j < un_dr_sl.n_rows; j++) {
      H(un_dr_sl.n_cols + i, j) = lower_left_matrix(i,j);
    }
  }
  arma::mat upper_right_matrix = ob_dr_sl.t() * inv_sq_ob_diff * ob_dr_sl;
  for(unsigned int i = 0; i < un_dr_sl.n_cols; i++) {
    for(unsigned int j = 0; j < un_dr_sl.n_cols; j++) {
      H(i, un_dr_sl.n_rows + j) = upper_right_matrix(i,j);
    }
  }
  /////////////////////////
  // `QZ` inplementation //
  /////////////////////////
  // dummy matrix for QZ
  arma::mat B(arma::size(H), arma::fill::eye);
  arma::mat AA;
  arma::mat BB;
  arma::mat Q;
  arma::mat Z;
  
  bool qz_res = qz(AA, BB, Q, Z, H, B, "rhp");
  if(qz_res == false) {
    stop("Failed in QZ decomposition in vcov calculation."); 
  }
  
  arma::mat generalized_eigenvec_mat = Q.t();
  
  arma::mat upper_right_q = generalized_eigenvec_mat.submat(0              , 0, un_dr_sl.n_rows - 1    , un_dr_sl.n_rows - 1);
  arma::mat lower_right_q = generalized_eigenvec_mat.submat(un_dr_sl.n_cols, 0, un_dr_sl.n_rows * 2 - 1, un_dr_sl.n_rows - 1);
  arma::mat gamma = lower_right_q * arma::inv(upper_right_q);
  return gamma;
}

// [[Rcpp::export]]
arma::mat calc_filter_mean(arma::cube un_dr_sl, arma::cube un_dr_in, arma::cube ob_dr_sl, arma::cube ob_dr_in, arma::cube inv_sq_ob_diff, arma::cube vcov, arma::vec init, double delta, arma::mat deltaY) {
    int d_un = un_dr_sl.n_rows; // the number of observed variables
    int d_ob = ob_dr_sl.n_rows; // the number of unobserved variables
    int n = deltaY.n_cols + 1;  // the number of observations

    // initialize mean with suitable size, no value
    arma::mat mean(d_un, n, arma::fill::none);
    mean.col(0) = init;
    for(int i = 1; i < n; i++){
        arma::vec mean_prev(&mean(0, i-1), d_un, false, true);
        arma::mat ob_dr_sl_slice(&ob_dr_sl(0, 0, i-1), d_ob, d_ob, false, true);
        arma::mat ob_dr_sl_slice_t = ob_dr_sl_slice.t();
        arma::vec ob_dr_in_slice(&ob_dr_in(0, 0, i-1), d_ob, false, true);
        arma::mat un_dr_sl_slice(&un_dr_sl(0, 0, i-1), d_un, d_ob, false, true);
        arma::vec un_dr_in_slice(&un_dr_in(0, 0, i-1), d_un, false, true);
        arma::mat vcov_slice(&vcov(0, 0, i-1), d_un, d_un, false, true);
        arma::mat inv_sq_ob_diff_slice(&inv_sq_ob_diff(0, 0, i-1), d_ob, d_ob, false, true);
        arma::vec deltaY_col(&deltaY(0, i-1), d_ob, false, true);
        arma::mat coef = vcov_slice * ob_dr_sl_slice_t * inv_sq_ob_diff_slice;
        mean.col(i) = mean_prev
                    + (un_dr_sl_slice * mean_prev + un_dr_in_slice - coef * (ob_dr_sl_slice * mean_prev + ob_dr_in_slice)) * delta
                    + coef * deltaY_col;
    }

    return mean;
}

// [[Rcpp::export]]
arma::mat calc_filter_mean_time_homogeneous_with_vcov_are(arma::mat un_dr_sl, arma::vec un_dr_in, arma::mat ob_dr_sl, arma::vec ob_dr_in, arma::mat inv_sq_ob_diff, arma::mat vcov, arma::vec init, double delta, arma::mat deltaY) {
  int d_un = un_dr_sl.n_rows; // the number of observed variables
  int d_ob = ob_dr_sl.n_rows; // the number of unobserved variables
  int n = deltaY.n_cols + 1;  // the number of observations
  
  // initialize mean with suitable size, no value
  arma::mat mean(d_un, n, arma::fill::none);
  mean.col(0) = init;
  
  arma::mat deltaY_coef = vcov * ob_dr_sl.t() * inv_sq_ob_diff;
  arma::mat mean_prev_coef = arma::eye(d_un, d_un) + (un_dr_sl - deltaY_coef * ob_dr_sl) * delta;
  arma::vec intercept = (un_dr_in - deltaY_coef * ob_dr_in) * delta;
  for(int i = 1; i < n; i++){
    arma::vec deltaY_col(&deltaY(0, i-1), d_ob, false, true);
    arma::vec mean_prev(&mean(0, i-1), d_un, false, true);
    mean.col(i) = mean_prev_coef * mean_prev + deltaY_coef * deltaY_col + intercept;
  }
  
  return mean;
}

// [[Rcpp::export]]
arma::mat calc_filter_mean_time_homogeneous(arma::mat un_dr_sl, arma::vec un_dr_in, arma::mat ob_dr_sl, arma::vec ob_dr_in, arma::mat inv_sq_ob_diff, arma::cube vcov, arma::vec init, double delta, arma::mat deltaY) {
  // initialize mean with suitable size, no value
  int d_un = un_dr_sl.n_rows; // the number of observed variables
  int d_ob = ob_dr_sl.n_rows; // the number of unobserved variables
  int n = deltaY.n_cols + 1;  // the number of observations
  arma::mat mean(d_un, n, arma::fill::none);
  mean.col(0) = init;
  
  arma::mat vcov_deltaY_coef = ob_dr_sl.t() * inv_sq_ob_diff;
  arma::mat mean_prev_coef = arma::eye(d_un, d_un) + un_dr_sl * delta;
  arma::mat vcov_mean_prev_coef = vcov_deltaY_coef * ob_dr_sl * delta;
  arma::vec vcov_coef = vcov_deltaY_coef * ob_dr_in * delta;
  arma::vec intercept = un_dr_in * delta;
  for(int i = 1; i < n; i++){
    arma::mat vcov_slice(&vcov(0, 0, i-1), d_un, d_un, false, true);
    arma::vec deltaY_col(&deltaY(0, i-1), d_ob, false, true);
    arma::vec mean_prev(&mean(0, i-1), d_un, false, true);
    mean.col(i) = mean_prev_coef * mean_prev + vcov_slice * (vcov_mean_prev_coef * mean_prev + vcov_mean_prev_coef * deltaY_col + vcov_coef) + intercept;
  }
  
  return mean;
}

//[[Rcpp::export]]
arma::mat calc_filter_mean_explicit(arma::mat un_dr_sl, arma::mat un_dr_in, arma::mat ob_dr_sl, arma::mat ob_dr_in, arma::mat inv_sq_ob_diff, arma::mat vcov, arma::vec init, double delta, arma::mat deltaY) {
    /*
    calculate mean explicitly if coefficients are time-independent.
    use when estimated vcov with Algebric Riccati Equation.
    
    m_{i+1} = e*m_{i} + e*a_2*h + e*\gamma(\theta1, \theta2)c(\theta2)^T\Sigma(\theta1)^{-1}\Delta_{i+1}Y + e*\gamma(\theta1, \theta2)c(\theta2)^T\Sigma(\theta1)^{-1}*c_2*h
    where 
    e = \exp(-\alpha(\theta1, \theta2)h)
    \alpha(\theta1, \theta2) = a(\theta1) + \gamma(\theta1, \theta2)c(\theta2)^T\Sigma(\theta1)c(\theta2)
    \Sigma = \sigma(\theta1)\sigma^T(\theta1)
    */
    int d_un = un_dr_sl.n_rows; // the number of observed variables
    int d_ob = ob_dr_sl.n_rows; // the number of unobserved variables
    int n_deltaY = deltaY.n_cols;  // the number of observations - 1

    // initialize mean with suitable size, no value
    arma::mat mean = arma::mat(d_un, n_deltaY + 1, arma::fill::none);
    mean.col(0) = init;

    // Compute exp_alpha_h, deltaY_coeff, and intercept within an inner scope
    arma::mat exp_alpha_h;
    arma::mat deltaY_coeff;
    arma::mat intercept;

    {
      // Intermediate calculations are scoped within this block
      arma::mat intermed = vcov * ob_dr_sl.t() * inv_sq_ob_diff;
      arma::mat alpha = un_dr_sl - intermed * ob_dr_sl;
      exp_alpha_h = arma::expmat(alpha * delta);
      deltaY_coeff = exp_alpha_h * intermed;
      intercept = (exp_alpha_h * un_dr_in + deltaY_coeff * ob_dr_in) * delta;
    }

    for(int i = 0; i < n_deltaY; i++){
        arma::vec mean_prev(&mean(0, i), d_un, false, true);
        arma::vec deltaY_col(&deltaY(0, i), d_ob, false, true);
        mean.col(i + 1) = exp_alpha_h * mean_prev + deltaY_coeff * deltaY_col + intercept;
    }

    return mean;
}
/*** R
un_dr_sl <- array(1:9, dim = c(3, 3))
un_dr_in <- array(1:3, dim = c(3, 1))
un_diff <- array(1:9, dim = c(3, 3))
ob_dr_sl <- array(1:9, dim = c(2, 3))
ob_dr_in <- array(1:3, dim = c(2, 1))
ob_diff <- array(c(1,0,0,1), dim = c(2, 2))
vcov = calc_filter_vcov_are(un_dr_sl, un_diff, ob_dr_sl, ob_diff)
mean_init <- c(1,1,1)
delta <- 0.1
deltaY <- matrix(rnorm(6), nrow=2)
calc_filter_mean_explicit(un_dr_sl, un_dr_in, ob_dr_sl, ob_dr_in, ob_diff, vcov, mean_init, delta, deltaY)
*/
/*
un_dr_sl <- array(rep(1:9, 4), dim = c(3, 3, 4))
un_dr_in <- array(rep(1:3, 4), dim = c(3, 1, 4))
un_diff <- array(rep(1:6, 4), dim = c(3, 3, 4))
ob_dr_sl <- array(rep(1:9, 4), dim = c(2, 3, 4))
ob_dr_in <- array(rep(1:3, 4), dim = c(2, 1, 4))
ob_diff <- array(rep(c(1,0,0,1), 4), dim = c(2, 2, 4))
vcov_init <- diag(3)
mean_init <- c(1,1,1)
delta <- 0.1
deltaY <- matrix(rnorm(6), nrow=2)
res1 <- calc_filter_vcov(un_dr_sl, un_diff, ob_dr_sl, ob_diff, vcov_init, delta)
vcov <- res1$vcov
inv_sq_ob_diff <- res1$inv_sq_ob_diff
calc_filter_mean(un_dr_sl, un_dr_in, un_diff, ob_dr_sl, ob_dr_in, ob_diff, vcov, inv_sq_ob_diff, mean_init, delta, deltaY)
*/
