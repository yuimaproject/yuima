// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]


List calc_filter_vcov(arma::cube un_dr_sl, arma::cube un_diff, arma::cube ob_dr_sl, arma::cube ob_diff, arma::mat init, double delta) {
    // initialize vcov with suitablesize, no value
    arma::cube vcov = arma::cube(un_dr_sl.n_rows, un_dr_sl.n_cols, un_dr_sl.n_slices, arma::fill::none);
    vcov.slice(0) = init;

    arma::cube inv_sq_ob_diff = arma::cube(ob_diff.n_rows, ob_diff.n_rows, ob_diff.n_slices, arma::fill::none);
    int n_data = un_dr_sl.n_slices;
    for(int i = 1; i < n_data; i++){
        inv_sq_ob_diff.slice(i-1) = arma::inv_sympd(ob_diff.slice(i-1) * ob_diff.slice(i-1).t());
        arma::mat vcov_prev = vcov.slice(i - 1);
        arma::mat vcov_next = vcov_prev + delta * (
            un_diff.slice(i-1) * un_diff.slice(i-1).t()
             + un_dr_sl.slice(i-1) * vcov_prev
             + vcov_prev * un_dr_sl.slice(i-1).t()
             - vcov_prev * ob_dr_sl.slice(i-1).t() * inv_sq_ob_diff.slice(i-1) * ob_dr_sl.slice(i-1) * vcov_prev);
        vcov.slice(i) = vcov_next;
    }

    //return vcov and inv_sq_ob_diff;
    List res = List::create(Named("vcov") = vcov , Named("inv_sq_ob_diff") = inv_sq_ob_diff);
    return res;
}

// [[Rcpp::export]]
arma::mat calc_filter_mean(arma::cube un_dr_sl, arma::cube un_dr_in, arma::cube un_diff, arma::cube ob_dr_sl, arma::cube ob_dr_in, arma::cube ob_diff, arma::cube vcov, arma::cube inv_sq_ob_diff, arma::vec init, double delta, arma::mat deltaY) {
    // initialize vcov with suitable size, no value
    arma::mat mean = arma::mat(un_dr_sl.n_rows, un_dr_sl.n_slices, arma::fill::none);
    mean.col(0) = init;
    int n_data = un_dr_sl.n_slices;
    for(int i = 1; i < n_data; i++){
        arma::vec mean_prev = mean.col(i - 1);
        arma::vec mean_next = mean_prev 
                           + (un_dr_sl.slice(i-1) * mean_prev + un_dr_in.slice(i-1) - vcov.slice(i-1) * ob_dr_sl.slice(i-1).t() * inv_sq_ob_diff.slice(i-1) * (ob_dr_sl.slice(i-1) * mean_prev + ob_dr_in.slice(i-1))) * delta
                           + vcov.slice(i-1) * ob_dr_sl.slice(i-1).t() * inv_sq_ob_diff.slice(i-1) * deltaY.col(i-1);
        mean.col(i) = mean_next;
    }

    return mean;
}

// [[Rcpp::export]]
arma::mat calc_filter_vcov_are(arma::mat un_dr_sl, arma::mat un_diff, arma::mat ob_dr_sl, arma::mat ob_diff) {
    arma::mat H(un_dr_sl.n_cols + un_dr_sl.n_rows, un_dr_sl.n_cols + un_dr_sl.n_rows);
    /*
    H = 
    \left(
        \begin{array}{cc}
        a(\theta_2)^\top & c(\theta_2)^\top\{\sigma(\theta_1)\sigma(\theta_1)^\top\}^{-1}c(\theta_2) \\
        b(\theta_2)b(\theta_2)^\top & -a(\theta_2) \\
        \end{array}
    \right), 
    a = \mathrm{un\_dr\_sl}, 
    b = \mathrm{un\_diff}, 
    c = \mathrm{ob\_dr\_sl}, 
    \sigma = \mathrm{ob\_diff}
    */
    for(int i = 0; i < un_dr_sl.n_rows; i++) {
        for(int j = 0; j < un_dr_sl.n_cols; j++) {
            H(j,i) = -un_dr_sl(i,j);
            H(un_dr_sl.n_cols + i, un_dr_sl.n_rows + j) = un_dr_sl(i,j);
        }
    }
    arma::mat lower_left_matrix = un_diff * un_diff.t();
    for(int i = 0; i < un_dr_sl.n_rows; i++) {
        for(int j = 0; j < un_dr_sl.n_rows; j++) {
            H(un_dr_sl.n_cols + i, j) = lower_left_matrix(i,j);
        }
    }
    arma::mat upper_right_matrix = ob_dr_sl.t() * arma::inv_sympd(ob_diff * ob_diff.t()) * ob_dr_sl;
    for(int i = 0; i < un_dr_sl.n_cols; i++) {
        for(int j = 0; j < un_dr_sl.n_cols; j++) {
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

//[[Rcpp::export]]
arma::mat calc_filter_mean_explicit(arma::mat un_dr_sl, arma::mat un_dr_in, arma::mat ob_dr_sl, arma::mat ob_dr_in, arma::mat ob_diff, arma::mat vcov, arma::vec init, double delta, arma::mat deltaY) {
    /*
    calculate mean explicitly if coefficients are time-independent.
    use when estimated vcov with Algebric Riccati Equation.
    
    m_{i+1} = e*m_{i} + e*a_2*h + e*\gamma(\theta1, \theta2)c(\theta2)^T\Sigma(\theta1)^{-1}\Delta_{i+1}Y + e*\gamma(\theta1, \theta2)c(\theta2)^T\Sigma(\theta1)^{-1}*c_2*h
    where 
    e = \exp(-\alpha(\theta1, \theta2)h)
    \alpha(\theta1, \theta2) = a(\theta1) + \gamma(\theta1, \theta2)c(\theta2)^T\Sigma(\theta1)c(\theta2)
    \Sigma = \sigma(\theta1)\sigma^T(\theta1)
    */
    // initialize mean with suitable size, no value
    int n_data = deltaY.n_cols + 1;
    arma::mat mean = arma::mat(un_dr_sl.n_rows, n_data, arma::fill::none);
    mean.col(0) = init;
    
    arma::mat Sigma = ob_diff * ob_diff.t();
    arma::mat inv_Sigma = arma::inv_sympd(Sigma);
    arma::mat alpha = un_dr_sl - vcov * ob_dr_sl.t() * inv_Sigma * ob_dr_sl;
    arma::mat exp_alpha_h = arma::expmat(alpha * delta);
    arma::mat deltaY_coeff = exp_alpha_h * vcov * ob_dr_sl.t() * inv_Sigma;

    for(int i = 1; i < n_data; i++){
        arma::vec mean_prev = mean.col(i - 1);
        arma::vec mean_next = exp_alpha_h * mean_prev + exp_alpha_h * un_dr_in * delta + deltaY_coeff * deltaY.col(i-1) + deltaY_coeff * ob_dr_in * delta;
        mean.col(i) = mean_next;
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
