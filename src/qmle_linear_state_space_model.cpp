// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
double minusloglcpp_linear_state_space_theta1(
    double logdet_sq_observed_diffusion, arma::mat inv_sq_observed_diffusion,
    arma::mat dx, double h, int drop_terms) {
  // observed_diffusion : 2-dim array of evaluated observed diffusion
  // coefficients dx : 2-dim matrix of delta of observerd variables in each time
  // point h : interval of observations
  int d_ob = dx.n_rows;
  int n_dx = dx.n_cols;

  // calculate quasi-log-likelihood
  double QL1 = 0;
  double QL2 = 0;
  for (int i = drop_terms; i < n_dx; i++) {
    arma::vec dx_col(&dx(0, i), d_ob, false, true);
    QL1 += logdet_sq_observed_diffusion;
    QL2 += arma::as_scalar(dx_col.t() * inv_sq_observed_diffusion * dx_col);
  }
  return -(QL1 + QL2 / h) / 2;
}

// [[Rcpp::export]]
double minusloglcpp_linear_state_space_theta2(
    arma::mat observed_drift, arma::cube inv_sq_observed_diffusion, arma::mat dx,
    double h, int drop_terms) {
  // observed_drift : 2-dim matrix of evaluated observed drift coefficients
  // inv_sq_observed_diffusion : (\sigma\sigma^T)^{-1} \sigma is 2-dim array of
  // evaluated observed diffusion coefficients). it is cube of n_slices = 1 for consistency.
  // dx : 2-dim matrix of delta of
  // observerd variables in each time point h : interval of observations
  int d_ob = observed_drift.n_rows;
  int n_dx = dx.n_cols;
  arma::mat inv_sq_observed_diffusion_slice(&inv_sq_observed_diffusion(0,0,0), inv_sq_observed_diffusion.n_rows, inv_sq_observed_diffusion.n_cols, false, true);
  
  // calculate quasi-log-likelihood
  double QL = 0;
  for (int i = drop_terms; i < n_dx; i++) {
    arma::vec observed_drift_col(&observed_drift(0, i), d_ob, false, true);
    arma::vec dx_col(&dx(0, i), d_ob, false, true);
    arma::vec tmp = observed_drift_col * h - dx_col;
    QL += arma::as_scalar(tmp.t() * inv_sq_observed_diffusion_slice * tmp);
  }
  return -QL / 2 / h;
}

// [[Rcpp::export]]
arma::cube calc_inverce_square(arma::cube cube) {
  // calculate (AA^T)^(-1) for each slice.
  arma::cube res(cube.n_rows, cube.n_rows, cube.n_slices);
  int slices = cube.n_slices;
  for (int i = 0; i < slices; i++) {
    arma::mat cube_slice(&cube(0, 0, i), cube.n_rows, cube.n_rows, false, true);
    res.slice(i) = arma::inv_sympd(cube_slice * cube_slice.t());
  }

  return res;
}