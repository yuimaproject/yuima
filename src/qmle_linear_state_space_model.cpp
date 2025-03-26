// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
double minusloglcpp_linear_state_space_theta1(double logdet_sq_ob_diff,
                                              arma::mat& inv_sq_ob_diff,
                                              arma::mat& deltaY, double delta) {
  // (sq_ob_diff : 2-dim array of evaluated observed diffusion matrix)
  // logdet_sq_ob_diff : log determinant of sq_ob_diff
  // inv_sq_ob_diff : 2-dim array of inverse of sq_ob_diff
  // deltaY : 2-dim matrix of delta of observed variables in each time point
  // delta : interval of observations

  int n_deltaY = deltaY.n_cols;  // the number of observations - 1

  // calculate quasi-log-likelihood
  double ql = logdet_sq_ob_diff * n_deltaY;
  ql += arma::trace(inv_sq_ob_diff * deltaY * deltaY.t()) / delta;
  return -ql * 0.5;
}

// [[Rcpp::export]]
arma::cube calc_inverse_square(arma::cube cube) {
  // calculate (AA^T)^(-1) for each slice.
  arma::cube res(cube.n_rows, cube.n_rows, cube.n_slices);
  int slices = cube.n_slices;
  for (int i = 0; i < slices; i++) {
    arma::mat cube_slice(&cube(0, 0, i), cube.n_rows, cube.n_rows, false, true);
    res.slice(i) = arma::inv_sympd(cube_slice * cube_slice.t());
  }

  return res;
}