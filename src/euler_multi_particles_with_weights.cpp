#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::cube euler_multi_particles_with_weights(arma::mat x0s, double t0, int r,
                                             double dt, int n, arma::vec dW,
                                             std::string modeltime,
                                             CharacterVector modelstate,
                                             SEXP drift, SEXP diffusion,
                                             Environment env, Environment rho) {
  int nsim = x0s.n_rows;  // number of particles
  int d = x0s.n_cols;     // number of dimensions

  // Create a d x (n+1) matrix to hold the solution trajectory
  arma::cube X(nsim, d, n + 1);
  // debug: print the dimensions of the cube
  X.slice(0) = x0s;
  // Initialize time parameter vector
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
      NumericVector b0 = as<NumericVector>(Rf_eval(VECTOR_ELT(drift, 0), rho));
      NumericVector sigma0 = as<NumericVector>(Rf_eval(VECTOR_ELT(diffusion, 0), rho));
      // Euler update: x(t+dt) = x(t) + drift*dt + diffusion*dW
      for (int j = 0; j < d; j++) {
        double temp = X(k, j, i) + b0[j] * dt;
        for (int l = 0; l < r; l++) {
          temp += sigma0[l + j * r] * dW[l + i * r];
        }
        X(k, j, i + 1) = temp;
      }
    }

    // Update time parameter
    t += dt;
  }

  return X;
}
