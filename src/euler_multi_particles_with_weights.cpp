#include <RcppArmadillo.h>
// #include <time.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List euler_multi_particles_with_weights(
    arma::mat x0s, arma::vec weight_init, double t0, int r, double dt, int n,
    arma::vec dW, std::string modeltime, CharacterVector modelstate,
    ExpressionVector observed_drift, ExpressionVector unobserved_drift,
    ExpressionVector observed_diffusion, ExpressionVector unobserved_diffusion,
    arma::mat deltaY, Environment env, Environment rho) {
  int nsim = x0s.n_rows;              // number of particles
  int d_obs = observed_drift.size();  // number of observed dimensions
  int d_unobs = x0s.n_cols;           // number of unobserved dimensions

  // Create a matrix of (nsim, d, n+1) to hold the solution trajectory
  arma::cube X(nsim, d_unobs, n + 1);
  X.slice(0) = x0s;

  arma::mat weights(nsim, n + 1);
  weights.col(0) = weight_init;

  double t = t0;

  // Time stepping loop
  for (int i = 0; i < n; i++) {
    // assign the current time variable to the environment
    rho.assign(modeltime, t);
    for (int k = 0; k < nsim; k++) {
      // assign the current state for each dimension to the environment
      for (int j = 0; j < d_unobs; j++) {
        rho.assign(as<std::string>(modelstate[j]), X(k, j, i));
      }

      // Evaluate the drift and diffusion expressions
      arma::vec observed_drift_values =
          as<arma::vec>(Rf_eval(observed_drift[0], rho));
      arma::vec unobserved_drift_values =
          as<arma::vec>(Rf_eval(unobserved_drift[0], rho));
      NumericVector vec_observed_diffusion_values_sexp =
          Rf_eval(observed_diffusion[0], rho);
      arma::mat t_observed_diffusion_values(
          vec_observed_diffusion_values_sexp.begin(), r, d_obs, false);
      NumericVector vec_unobserved_diffusion_values_sexp =
          Rf_eval(unobserved_diffusion[0], rho);
      arma::mat t_unobserved_diffusion_values(
          vec_unobserved_diffusion_values_sexp.begin(), r, d_unobs, false);

      // Euler update: x(t+dt) = x(t) + drift*dt + diffusion*dW
      for (int j = 0; j < d_unobs; j++) {
        double temp = X(k, j, i) + unobserved_drift_values(j) * dt;
        for (int l = 0; l < r; l++) {
          temp +=
              t_unobserved_diffusion_values(l, j) * dW[l + i * r + k * n * r];
        }
        X(k, j, i + 1) = temp;
      }

      // Update weights
      weights(k, i + 1) =
          weights(k, i) +
          weights(k, i) *
              arma::as_scalar(observed_drift_values.t() *
                              arma::inv_sympd(t_observed_diffusion_values.t() *
                                              t_observed_diffusion_values) *
                              deltaY.col(i));
    }
    // Update time parameter
    t += dt;
  }
  return Rcpp::List::create(Rcpp::Named("values") = X,
                            Rcpp::Named("weights") = weights);
}
