#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube euler_multi_particles_with_weights(
    arma::mat x0s, arma::vec weight_init, double t0, int r, double dt, int n,
    arma::vec dW, std::string modeltime, CharacterVector modelstate,
    ExpressionVector observed_drift, ExpressionVector unobserved_drift,
    ExpressionVector observed_diffusion, ExpressionVector unobserved_diffusion,
    Environment env, Environment rho) {
  int nsim = x0s.n_rows;                  // number of particles
  int d_obs = observed_drift.size();      // number of observed dimensions
  int d_unobs = unobserved_drift.size();  // number of unobserved dimensions
  int d = x0s.n_cols;                     // number of dimensions

  // Create a matrix of (nsim, d, n+1) to hold the solution trajectory
  arma::cube X(nsim, d, n + 1);
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
      for (int j = 0; j < d; j++) {
        rho.assign(as<std::string>(modelstate[j]), X(k, j, i));
      }

      // Evaluate the drift and diffusion expressions
      NumericVector drift_values(d);
      NumericVector diffusion_values(d * r);
      for (int j = 0; j < d_obs; j++) {
        drift_values[j] = as<double>(Rf_eval(observed_drift[j], rho));
        for (int l = 0; l < r; l++) {
          diffusion_values[l + j * r] =
              as<double>(Rf_eval(observed_diffusion[l + j * r], rho));
        }
      }
      for (int j = 0; j < d_unobs; j++) {
        drift_values[j + d_obs] = as<double>(Rf_eval(unobserved_drift[j], rho));
        for (int l = 0; l < r; l++) {
          diffusion_values[l + (j + d_obs) * r] =
              as<double>(Rf_eval(unobserved_diffusion[l + j * r], rho));
        }
      }

      // Euler update: x(t+dt) = x(t) + drift*dt + diffusion*dW
      for (int j = 0; j < d; j++) {
        double temp = X(k, j, i) + drift_values[j] * dt;
        for (int l = 0; l < r; l++) {
          temp += diffusion_values[l + j * r] * dW[l + i * r + k * n * r];
        }
        X(k, j, i + 1) = temp;
      }
    }

    // Update weights

    // Update time parameter
    t += dt;
  }

  return X;
}
