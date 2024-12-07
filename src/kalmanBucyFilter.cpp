// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

arma::mat calc_filter_vcov_are(const arma::mat& un_dr_sl,
                               const arma::mat& un_diff,
                               const arma::mat& ob_dr_sl,
                               const arma::mat& inv_sq_ob_diff) {
  unsigned int m = un_dr_sl.n_cols;
  unsigned int n = un_dr_sl.n_rows;

  arma::mat H = arma::join_vert(
      arma::join_horiz(-un_dr_sl.t(), ob_dr_sl.t() * inv_sq_ob_diff * ob_dr_sl),
      arma::join_horiz(un_diff * un_diff.t(), un_dr_sl));

  // QZ decomposition
  arma::mat B(H.n_rows, H.n_cols, arma::fill::eye);
  arma::mat AA, BB, Q, Z;

  bool qz_res = arma::qz(AA, BB, Q, Z, H, B, "rhp");
  if (!qz_res) {
    Rcpp::stop("Failed in QZ decomposition in vcov calculation.");
  }

  arma::mat generalized_eigenvec_mat = Q.t();

  arma::mat upper_right_q = generalized_eigenvec_mat.submat(0, 0, n - 1, n - 1);
  arma::mat lower_right_q =
      generalized_eigenvec_mat.submat(m, 0, m + n - 1, n - 1);

  arma::mat gamma = lower_right_q * arma::inv(upper_right_q);

  return gamma;
}

std::tuple<arma::mat, double> calc_filter_mean_time_homogeneous_with_vcov_are(
    arma::mat& un_dr_sl, arma::vec& un_dr_in, arma::mat& ob_dr_sl,
    arma::vec& ob_dr_in, arma::mat& inv_sq_ob_diff, arma::mat& vcov,
    arma::vec& init, double delta, arma::mat& deltaY, bool calc_minuslogl,
    int drop_terms) {
  int d_un = un_dr_sl.n_rows;  // the number of unobserved variables
  int d_ob = ob_dr_sl.n_rows;  // the number of observed variables
  int n = deltaY.n_cols + 1;   // the number of observations

  // initialize mean with suitable size, no value
  arma::mat mean(d_un, n, arma::fill::none);
  mean.col(0) = init;

  double minuslogl = 0;

  arma::mat deltaY_coef = vcov * ob_dr_sl.t() * inv_sq_ob_diff;
  arma::mat mean_prev_coef =
      arma::eye(d_un, d_un) + (un_dr_sl - deltaY_coef * ob_dr_sl) * delta;
  arma::vec intercept = (un_dr_in - deltaY_coef * ob_dr_in) * delta;
  for (int i = 1; i < n; i++) {
    arma::vec deltaY_col(&deltaY(0, i - 1), d_ob, false, true);
    arma::vec mean_prev(&mean(0, i - 1), d_un, false, true);
    mean.col(i) =
        mean_prev_coef * mean_prev + deltaY_coef * deltaY_col + intercept;
    if (calc_minuslogl) {
      if (drop_terms < i) {
        arma::vec tmp = (ob_dr_sl * mean_prev + ob_dr_in) * delta - deltaY_col;
        minuslogl += arma::as_scalar(tmp.t() * inv_sq_ob_diff * tmp);
      }
    }
  }
  minuslogl = minuslogl / 2 / delta;

  return std::forward_as_tuple(mean, minuslogl);
}

std::tuple<arma::mat, double> calc_filter_mean_explicit(
    const arma::mat& un_dr_sl, const arma::vec& un_dr_in,
    const arma::mat& ob_dr_sl, const arma::vec& ob_dr_in,
    const arma::mat& inv_sq_ob_diff, const arma::mat& vcov,
    const arma::vec& init, double delta, arma::mat deltaY, bool calc_minuslogl,
    int drop_terms) {
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
  int d_un = un_dr_sl.n_rows;    // the number of unobserved variables
  int d_ob = ob_dr_sl.n_rows;    // the number of observed variables
  int n_deltaY = deltaY.n_cols;  // the number of observations - 1

  // initialize mean with suitable size, no value
  arma::mat mean = arma::mat(d_un, n_deltaY + 1, arma::fill::none);
  mean.col(0) = init;
  double minuslogl;

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

  for (int i = 0; i < n_deltaY; i++) {
    arma::vec mean_prev(&mean(0, i), d_un, false, true);
    arma::vec deltaY_col(&deltaY(0, i), d_ob, false, true);

    mean.col(i + 1) =
        exp_alpha_h * mean_prev + deltaY_coeff * deltaY_col + intercept;
    if (calc_minuslogl) {
      if (drop_terms < i) {
        arma::vec tmp = (ob_dr_sl * mean_prev + ob_dr_in) * delta - deltaY_col;
        minuslogl += arma::as_scalar(tmp.t() * inv_sq_ob_diff * tmp);
      }
    }
  }
  minuslogl = minuslogl / 2 / delta;

  return std::forward_as_tuple(mean, minuslogl);
}

Rcpp::List calc_kalman_bucy_filter_no_are_no_time_homogeneous(
    arma::cube un_dr_sl, arma::mat un_dr_in, arma::cube un_diff,
    arma::cube ob_dr_sl, arma::mat ob_dr_in, arma::cube inv_sq_ob_diff,
    arma::mat vcov_init, arma::vec mean_init, double delta, arma::mat deltaY,
    int upsamp_rate = 1) {
  int d_un = un_dr_sl.n_rows;  // the number of unobserved variables
  int d_ob = ob_dr_sl.n_rows;  // the number of observed variables
  int n = deltaY.n_cols + 1;   // the number of observations

  double upsamp_delta =
      delta / upsamp_rate;  // the interval of upsampled observations

  // initialize vcov and mean with suitable size, no value
  arma::cube vcov(d_un, d_un, n, arma::fill::none);
  vcov.slice(0) = vcov_init;
  arma::mat mean(d_un, n, arma::fill::none);
  mean.col(0) = mean_init;

  // temporary object for interpolated vcovs.
  arma::mat tmp_vcov = vcov_init;

  for (int i = 1; i < n; i++) {
    // calc mean
    int prev_i = (i - 1) * upsamp_rate;
    arma::vec mean_prev(&mean(0, i - 1), d_un, false, true);
    arma::mat ob_dr_sl_slice(&ob_dr_sl(0, 0, prev_i), d_ob, d_ob, false, true);
    arma::vec ob_dr_in_col(&ob_dr_in(0, prev_i), d_ob, false, true);
    arma::mat un_dr_sl_slice(&un_dr_sl(0, 0, prev_i), d_un, d_ob, false, true);
    arma::vec un_dr_in_col(&un_dr_in(0, prev_i), d_un, false, true);
    arma::mat inv_sq_ob_diff_slice(&inv_sq_ob_diff(0, 0, prev_i), d_ob, d_ob,
                                   false, true);
    arma::vec deltaY_col(&deltaY(0, i - 1), d_ob, false, true);
    arma::mat coef = tmp_vcov * ob_dr_sl_slice.t() * inv_sq_ob_diff_slice;
    mean.col(i) = mean_prev +
                  (un_dr_sl_slice * mean_prev + un_dr_in_col -
                   coef * (ob_dr_sl_slice * mean_prev + ob_dr_in_col)) *
                      delta +
                  coef * deltaY_col;

    // calc vcov
    for (int j = 0; j < upsamp_rate; j++) {
      int upsamp_i = (i - 1) * upsamp_rate + j;

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
    vcov.slice(i) = tmp_vcov;
  }

  return Rcpp::List::create(Rcpp::Named("vcov") = vcov,
                            Rcpp::Named("mean") = mean,
                            Rcpp::Named("minuslogl") = 0.0);
}

Rcpp::List calc_kalman_bucy_filter_time_homogeneous(
    arma::mat un_dr_sl, arma::vec un_dr_in, arma::mat un_diff,
    arma::mat ob_dr_sl, arma::vec ob_dr_in, arma::mat inv_sq_ob_diff,
    arma::mat vcov_init, arma::vec mean_init, double delta, arma::mat deltaY,
    int upsamp_rate = 1) {
  int d_un = un_dr_sl.n_rows;  // the number of unobserved variables
  int d_ob = ob_dr_sl.n_rows;  // the number of observed variables
  int n = deltaY.n_cols + 1;   // the number of observations

  double upsamp_delta =
      delta / upsamp_rate;  // the interval of upsampled observations

  // initialize vcov with suitable size, no value
  arma::cube vcov(d_un, d_un, n, arma::fill::none);
  vcov.slice(0) = vcov_init;
  arma::mat mean(d_un, n, arma::fill::none);
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

  for (int i = 1; i < n; i++) {
    // calc mean
    arma::vec deltaY_col(&deltaY(0, i - 1), d_ob, false, true);
    arma::vec mean_prev(&mean(0, i - 1), d_un, false, true);
    mean.col(i) =
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
    vcov.slice(i) = tmp_vcov;
  }

  return Rcpp::List::create(Rcpp::Named("vcov") = vcov,
                            Rcpp::Named("mean") = mean,
                            Rcpp::Named("minuslogl") = 0.0);
}

// [[Rcpp::export]]
Rcpp::List calc_kalman_bucy_filter_cpp(
    arma::cube un_dr_sl, arma::mat un_dr_in, arma::cube un_diff,
    arma::cube ob_dr_sl, arma::mat ob_dr_in, arma::cube inv_sq_ob_diff,
    arma::mat vcov_init, arma::vec mean_init, double delta, arma::mat deltaY,
    bool use_are, bool is_explicit, bool is_time_homogeneous,
    bool calc_minuslogl, int drop_terms, int upsump_rate = 1) {
  int d_un = un_dr_sl.n_rows;  // the number of observed variables
  int d_ob = ob_dr_sl.n_rows;  // the number of unobserved variables
  // int drop_terms = 0;
  arma::cube vcov;
  arma::mat mean;
  double minuslogl;
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
    arma::mat vcov_slice = calc_filter_vcov_are(
        un_dr_sl_slice, un_diff_slice, ob_dr_sl_slice, inv_sq_ob_diff_slice);

    // calc mean
    if (is_explicit) {
      std::tie(mean, minuslogl) = calc_filter_mean_explicit(
          un_dr_sl_slice, un_dr_in_col, ob_dr_sl_slice, ob_dr_in_col,
          inv_sq_ob_diff_slice, vcov_slice, mean_init, delta, deltaY,
          calc_minuslogl, drop_terms);
    } else {
      std::tie(mean, minuslogl) =
          calc_filter_mean_time_homogeneous_with_vcov_are(
              un_dr_sl_slice, un_dr_in_col, ob_dr_sl_slice, ob_dr_in_col,
              inv_sq_ob_diff_slice, vcov_slice, mean_init, delta, deltaY,
              calc_minuslogl, drop_terms);
    }

    // make a cube of vcov for consistency
    vcov = arma::cube(vcov_slice.n_rows, vcov_slice.n_cols, 1);
    vcov.slice(0) = vcov_slice;

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
