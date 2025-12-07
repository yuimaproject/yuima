// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

double calc_minuslogl_time_homogeneous(const arma::mat& ob_dr_sl,
                                       const arma::mat& ob_dr_in,
                                       const arma::mat& inv_sq_ob_diff,
                                       const arma::mat& deltaY,
                                       const arma::mat& mean, double delta,
                                       unsigned int drop_terms) {
  unsigned int d_un = ob_dr_sl.n_cols;    // the number of unobserved variables
  unsigned int d_ob = ob_dr_sl.n_rows;    // the number of observed variables
  unsigned int n_deltaY = deltaY.n_cols;  // the number of observations - 1
  arma::mat tmp =
      (ob_dr_sl * mean.submat(0, drop_terms, d_un - 1, n_deltaY - 1) +
       ob_dr_in * arma::ones(1, n_deltaY - drop_terms)) *
          delta -
      deltaY.submat(0, drop_terms, d_ob - 1, n_deltaY - 1);
  return arma::trace(inv_sq_ob_diff * tmp * tmp.t()) * 0.5 / delta;
}

void calc_filter_vcov_are(arma::mat& vcov, const arma::mat& un_dr_sl,
                          const arma::mat& un_diff, const arma::mat& ob_dr_sl,
                          const arma::mat& inv_sq_ob_diff) {
  unsigned int d_un = un_dr_sl.n_cols;  // the number of unobserved variables
  unsigned int d_qz = 2 * d_un;         // dimension for QZ decomposition

  arma::mat H = arma::join_vert(
      arma::join_horiz(-un_dr_sl.t(), ob_dr_sl.t() * inv_sq_ob_diff * ob_dr_sl),
      arma::join_horiz(un_diff * un_diff.t(), un_dr_sl));

  // QZ decomposition
  arma::mat B = arma::eye(d_qz, d_qz);
  arma::mat AA(d_qz, d_qz), BB(d_qz, d_qz), Q(d_qz, d_qz), Z(d_qz, d_qz);

  bool qz_res = arma::qz(AA, BB, Q, Z, H, B, "rhp");
  if (!qz_res) {
    Rcpp::stop("Failed in QZ decomposition in vcov calculation.");
  }

  arma::mat generalized_eigenvec_mat = Q.t();

  arma::mat upper_right_q =
      generalized_eigenvec_mat.submat(0, 0, d_un - 1, d_un - 1);
  arma::mat lower_right_q =
      generalized_eigenvec_mat.submat(d_un, 0, d_qz - 1, d_un - 1);

  vcov = arma::solve(upper_right_q.t(), lower_right_q.t()).t();
}

void calc_filter_mean_time_homogeneous_with_vcov_are(
    arma::mat& mean, const arma::mat& un_dr_sl, const arma::vec& un_dr_in,
    const arma::mat& ob_dr_sl, const arma::vec& ob_dr_in,
    const arma::mat& inv_sq_ob_diff, const arma::mat& vcov,
    const arma::vec& init, double delta, arma::mat& deltaY) {
  unsigned int d_un = un_dr_sl.n_rows;    // the number of unobserved variables
  unsigned int n_deltaY = deltaY.n_cols;  // the number of observations - 1

  mean.col(0) = init;

  arma::mat deltaY_coef = vcov * ob_dr_sl.t() * inv_sq_ob_diff;
  arma::mat deltaY_term = deltaY_coef * deltaY;
  arma::mat mean_prev_coef =
      arma::eye(d_un, d_un) + (un_dr_sl - deltaY_coef * ob_dr_sl) * delta;
  arma::vec intercept = (un_dr_in - deltaY_coef * ob_dr_in) * delta;
  arma::mat deltaY_term_plus_intercept =
      deltaY_term.each_col() + (un_dr_in - deltaY_coef * ob_dr_in) * delta;
  ;

  for (int i = 0; i < n_deltaY; i++) {
    arma::vec deltaY_term_plus_intercept_col(&deltaY_term_plus_intercept(0, i),
                                             d_un, false, true);
    arma::vec mean_prev(&mean(0, i), d_un, false, true);
    mean.col(i + 1) =
        mean_prev_coef * mean_prev + deltaY_term_plus_intercept_col;
  }
}

void calc_filter_mean_explicit(arma::mat& mean, const arma::mat& un_dr_sl,
                               const arma::vec& un_dr_in,
                               const arma::mat& ob_dr_sl,
                               const arma::vec& ob_dr_in,
                               const arma::mat& inv_sq_ob_diff,
                               const arma::mat& vcov, const arma::vec& init,
                               double delta, arma::mat& deltaY) {
  /*
  calculate mean explicitly if coefficients are time-independent.
  use when estimated vcov with Algebric Riccati Equation.

  m_{i+1} = e*m_{i} + e*a_2*h + e*\gamma(\theta1,
  \theta2)c(\theta2)^T\Sigma(\theta1)^{-1}\Delta_{i+1}Y + e*\gamma(\theta1,
  \theta2)c(\theta2)^T\Sigma(\theta1)^{-1}*c_2*h where e = \exp(-\alpha(\theta1,
  \theta2)h)
  \alpha(\theta1, \theta2) = a(\theta1) + \gamma(\theta1,
  \theta2)c(\theta2)^T\Sigma(\theta1)c(\theta2)
  \Sigma = \sigma(\theta1)\sigma^T(\theta1)
  */
  unsigned int d_un = un_dr_sl.n_rows;    // the number of unobserved variables
  unsigned int n_deltaY = deltaY.n_cols;  // the number of observations - 1

  mean.col(0) = init;

  arma::mat intermed = vcov * ob_dr_sl.t() * inv_sq_ob_diff;
  arma::mat exp_alpha_h =
      arma::expmat((un_dr_sl - intermed * ob_dr_sl) * delta);
  arma::mat deltaY_coeff = exp_alpha_h * intermed;
  arma::mat deltaY_term = deltaY_coeff * deltaY;
  arma::mat deltaY_term_plus_intercept =
      deltaY_term.each_col() +
      (exp_alpha_h * un_dr_in + deltaY_coeff * ob_dr_in) * delta;

  for (int i = 0; i < n_deltaY; i++) {
    arma::vec mean_prev(&mean(0, i), d_un, false, true);
    arma::vec deltaY_term_plus_intercept_col(&deltaY_term_plus_intercept(0, i),
                                             d_un, false, true);
    mean.col(i + 1) = exp_alpha_h * mean_prev + deltaY_term_plus_intercept_col;
  }
}

Rcpp::List calc_kalman_bucy_filter_no_are_no_time_homogeneous(
    arma::cube& un_dr_sl, arma::mat& un_dr_in, arma::cube& un_diff,
    arma::cube& ob_dr_sl, arma::mat& ob_dr_in, arma::cube& inv_sq_ob_diff,
    arma::mat& vcov_init, arma::vec& mean_init, double delta, arma::mat& deltaY,
    unsigned int upsamp_rate = 1) {
  unsigned int d_un = un_dr_sl.n_rows;    // the number of unobserved variables
  unsigned int d_ob = ob_dr_sl.n_rows;    // the number of observed variables
  unsigned int n_deltaY = deltaY.n_cols;  // the number of observations - 1

  double upsamp_delta =
      delta / upsamp_rate;  // the interval of upsampled observations

  // initialize vcov and mean with suitable size, and set initial values.
  arma::cube vcov(d_un, d_un, n_deltaY + 1, arma::fill::none);
  vcov.slice(0) = vcov_init;
  arma::mat mean(d_un, n_deltaY + 1, arma::fill::none);
  mean.col(0) = mean_init;

  // temporary object for interpolated vcovs.
  arma::mat tmp_vcov = vcov_init;

  for (int i = 0; i < n_deltaY; i++) {
    // calc mean
    int prev_i = i * upsamp_rate;
    arma::vec mean_prev(&mean(0, i), d_un, false, true);
    arma::mat ob_dr_sl_slice(&ob_dr_sl(0, 0, prev_i), d_ob, d_ob, false, true);
    arma::vec ob_dr_in_col(&ob_dr_in(0, prev_i), d_ob, false, true);
    arma::mat un_dr_sl_slice(&un_dr_sl(0, 0, prev_i), d_un, d_ob, false, true);
    arma::vec un_dr_in_col(&un_dr_in(0, prev_i), d_un, false, true);
    arma::mat inv_sq_ob_diff_slice(&inv_sq_ob_diff(0, 0, prev_i), d_ob, d_ob,
                                   false, true);
    arma::vec deltaY_col(&deltaY(0, i), d_ob, false, true);
    arma::mat coef = tmp_vcov * ob_dr_sl_slice.t() * inv_sq_ob_diff_slice;
    mean.col(i + 1) = mean_prev +
                      (un_dr_sl_slice * mean_prev + un_dr_in_col -
                       coef * (ob_dr_sl_slice * mean_prev + ob_dr_in_col)) *
                          delta +
                      coef * deltaY_col;

    // calc vcov
    for (int j = 0; j < upsamp_rate; j++) {
      int upsamp_i = i * upsamp_rate + j;

      arma::mat un_diff_slice(&un_diff(0, 0, upsamp_i), d_un, un_diff.n_cols,
                              false, true);
      arma::mat un_dr_sl_slice(&un_dr_sl(0, 0, upsamp_i), d_un, d_un, false,
                               true);
      arma::mat ob_dr_sl_slice(&ob_dr_sl(0, 0, upsamp_i), d_ob, d_un, false,
                               true);
      arma::mat inv_sq_ob_diff_slice(&inv_sq_ob_diff(0, 0, upsamp_i), d_ob,
                                     d_ob, false, true);
      arma::mat first_order_term = un_dr_sl_slice * tmp_vcov;
      tmp_vcov =
          tmp_vcov +
          upsamp_delta * (un_diff_slice * un_diff_slice.t() + first_order_term +
                          first_order_term.t() -
                          tmp_vcov * ob_dr_sl_slice.t() * inv_sq_ob_diff_slice *
                              ob_dr_sl_slice * tmp_vcov);
    }
    vcov.slice(i + 1) = tmp_vcov;
  }

  return Rcpp::List::create(Rcpp::Named("vcov") = vcov,
                            Rcpp::Named("mean") = mean,
                            Rcpp::Named("minuslogl") = 0.0);
}

Rcpp::List calc_kalman_bucy_filter_time_homogeneous(
    arma::mat& un_dr_sl, arma::vec& un_dr_in, arma::mat& un_diff,
    arma::mat& ob_dr_sl, arma::vec& ob_dr_in, arma::mat& inv_sq_ob_diff,
    arma::mat& vcov_init, arma::vec& mean_init, double delta, arma::mat& deltaY,
    unsigned int upsamp_rate = 1) {
  unsigned int d_un = un_dr_sl.n_rows;    // the number of unobserved variables
  unsigned int d_ob = ob_dr_sl.n_rows;    // the number of observed variables
  unsigned int n_deltaY = deltaY.n_cols;  // the number of observations - 1

  double upsamp_delta =
      delta / upsamp_rate;  // the interval of upsampled observations

  // initialize vcov with suitable size, no value
  arma::cube vcov(d_un, d_un, n_deltaY + 1);
  vcov.slice(0) = vcov_init;
  arma::mat mean(d_un, n_deltaY + 1);
  mean.col(0) = mean_init;

  // coefficients for mean calc
  arma::mat mean_vcov_deltaY_coef = ob_dr_sl.t() * inv_sq_ob_diff;
  arma::mat mean_mean_coef = arma::eye(d_un, d_un) + un_dr_sl * delta;
  arma::mat mean_vcov_mean_coef = -mean_vcov_deltaY_coef * ob_dr_sl * delta;
  arma::vec mean_vcov_coef = -mean_vcov_deltaY_coef * ob_dr_in * delta;
  arma::vec mean_intercept = un_dr_in * delta;

  // coefficients for vcov calc
  arma::mat vcov_intercept = upsamp_delta * un_diff * un_diff.t();
  arma::mat vcov_second_order_coef =
      upsamp_delta * ob_dr_sl.t() * inv_sq_ob_diff * ob_dr_sl;
  arma::mat vcov_first_order_coef = upsamp_delta * un_dr_sl;

  // temporary object for interpolated vcovs.
  arma::mat tmp_vcov = vcov_init;

  for (int i = 0; i < n_deltaY; i++) {
    // calc mean
    arma::vec deltaY_col(&deltaY(0, i), d_ob, false, true);
    arma::vec mean_prev(&mean(0, i), d_un, false, true);
    mean.col(i + 1) =
        mean_mean_coef * mean_prev +
        tmp_vcov * (mean_vcov_mean_coef * mean_prev +
                    mean_vcov_deltaY_coef * deltaY_col + mean_vcov_coef) +
        mean_intercept;

    // calc vcov
    for (int j = 0; j < upsamp_rate; j++) {
      arma::mat first_order_term = vcov_first_order_coef * tmp_vcov;
      tmp_vcov = tmp_vcov + first_order_term + first_order_term.t() -
                 tmp_vcov * vcov_second_order_coef * tmp_vcov + vcov_intercept;
    }
    vcov.slice(i + 1) = tmp_vcov;
  }

  return Rcpp::List::create(Rcpp::Named("vcov") = vcov,
                            Rcpp::Named("mean") = mean,
                            Rcpp::Named("minuslogl") = 0.0);
}

// [[Rcpp::export]]
Rcpp::List calc_kalman_bucy_filter_cpp(
    arma::cube& un_dr_sl, arma::mat& un_dr_in, arma::cube& un_diff,
    arma::cube& ob_dr_sl, arma::mat& ob_dr_in, arma::cube& inv_sq_ob_diff,
    arma::mat& vcov_init, arma::vec& mean_init, double delta, arma::mat& deltaY,
    bool use_are, bool is_explicit, bool is_time_homogeneous,
    bool calc_minuslogl, unsigned int drop_terms,
    unsigned int upsump_rate = 1) {
  unsigned int d_un = un_dr_sl.n_rows;  // the number of observed variables
  unsigned int d_ob = ob_dr_sl.n_rows;  // the number of unobserved variables
  if (use_are) {
    // coefficients of SDE are independent of time.
    // So n_slices of coefficients should be 1.
    arma::mat un_dr_sl_slice(&un_dr_sl(0, 0, 0), d_un, d_un, false, true);
    arma::vec un_dr_in_col(&un_dr_in(0, 0), d_un, false, true);
    arma::mat un_diff_slice(&un_diff(0, 0, 0), d_un, un_diff.n_cols, false,
                            true);
    arma::mat ob_dr_sl_slice(&ob_dr_sl(0, 0, 0), d_ob, d_un, false, true);
    arma::vec ob_dr_in_col(&ob_dr_in(0, 0), d_ob, false, true);
    arma::mat inv_sq_ob_diff_slice(&inv_sq_ob_diff(0, 0, 0), d_ob, d_ob, false,
                                   true);

    // calc vcov
    arma::cube vcov(d_un, d_un, 1);
    arma::mat vcov_slice(d_un, d_un);
    calc_filter_vcov_are(vcov_slice, un_dr_sl_slice, un_diff_slice,
                         ob_dr_sl_slice, inv_sq_ob_diff_slice);
    vcov.slice(0) = vcov_slice;

    // calc mean
    arma::mat mean(d_un, deltaY.n_cols + 1);
    if (is_explicit) {
      calc_filter_mean_explicit(
          mean, un_dr_sl_slice, un_dr_in_col, ob_dr_sl_slice, ob_dr_in_col,
          inv_sq_ob_diff_slice, vcov_slice, mean_init, delta, deltaY);
    } else {
      calc_filter_mean_time_homogeneous_with_vcov_are(
          mean, un_dr_sl_slice, un_dr_in_col, ob_dr_sl_slice, ob_dr_in_col,
          inv_sq_ob_diff_slice, vcov_slice, mean_init, delta, deltaY);
    }

    // calc minuslogl
    double minuslogl = 0.0;
    if (calc_minuslogl) {
      minuslogl = calc_minuslogl_time_homogeneous(ob_dr_sl_slice, ob_dr_in_col,
                                                  inv_sq_ob_diff_slice, deltaY,
                                                  mean, delta, drop_terms);
    }

    // make a cube of vcov for consistency
    return Rcpp::List::create(Rcpp::Named("vcov") = vcov,
                              Rcpp::Named("mean") = mean,
                              Rcpp::Named("minuslogl") = minuslogl);
  } else {
    if (calc_minuslogl) {
      Rcpp::warning("minuslogl is not yet implemented for are=FALSE.");
    }
    if (is_time_homogeneous) {
      // coefficients of SDE are independent of time.
      // So n_slices of coefficients should be 1.
      arma::mat un_dr_sl_slice(&un_dr_sl(0, 0, 0), d_un, d_un, false, true);
      arma::vec un_dr_in_col(&un_dr_in(0, 0), d_un, false, true);
      arma::mat un_diff_slice(&un_diff(0, 0, 0), d_un, un_diff.n_cols, false,
                              true);
      arma::mat ob_dr_sl_slice(&ob_dr_sl(0, 0, 0), d_ob, d_un, false, true);
      arma::vec ob_dr_in_col(&ob_dr_in(0, 0), d_ob, false, true);
      arma::mat inv_sq_ob_diff_slice(&inv_sq_ob_diff(0, 0, 0), d_ob, d_ob,
                                     false, true);

      return calc_kalman_bucy_filter_time_homogeneous(
          un_dr_sl_slice, un_dr_in_col, un_diff_slice, ob_dr_sl_slice,
          ob_dr_in_col, inv_sq_ob_diff_slice, vcov_init, mean_init, delta,
          deltaY, upsump_rate);
    } else {
      return calc_kalman_bucy_filter_no_are_no_time_homogeneous(
          un_dr_sl, un_dr_in, un_diff, ob_dr_sl, ob_dr_in, inv_sq_ob_diff,
          vcov_init, mean_init, delta, deltaY, upsump_rate);
    }
  }
}
