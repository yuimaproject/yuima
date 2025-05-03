#include <RcppArmadillo.h>
// #include <time.h>

using namespace Rcpp;

arma::cube euler_multi_particles(
    arma::mat xinits, double t0, double dt, int step0, int steps, 
    int d_random, arma::vec dW,
    std::string time_var, CharacterVector unobserved_vars,
    ExpressionVector unobserved_drift,
    ExpressionVector unobserved_diffusion,
    Environment eval_env
){
  int num_particles = xinits.n_rows;
  int d_unob = xinits.n_cols;
  int random_index_offset = (step0 - 1) * num_particles * d_random;
  
  arma::cube X(num_particles, d_unob, steps + 1);
  X.slice(0) = xinits;

  for (int i = 0; i < steps; i++) { 
    eval_env.assign(time_var, t0 + (i + 1) * dt);
    for (int p = 0; p< num_particles; p++) {
      // assign
      for (int v = 0; v < d_unob; v++) {
        eval_env.assign(as<std::string>(unobserved_vars[v]), X(p, v, i));
      }

      // evaluate drift and diffusion
      arma::vec unobserved_drift_values =
        as<arma::vec>(Rf_eval(unobserved_drift[0], eval_env));
      NumericVector vec_unobserved_diffusion_values_sexp =
        Rf_eval(unobserved_diffusion[0], eval_env);
      arma::mat t_unobserved_diffusion_values(
          vec_unobserved_diffusion_values_sexp.begin(), d_random, d_unob, false); // TODO: consider transpose
      
      // simulate next step
      arma::uword start_index = static_cast<arma::uword>(random_index_offset + (i * num_particles + p) * d_random);
      arma::vec random_vec(dW.memptr() + start_index, d_random, false);
      X.slice(i+1).row(p) = X.slice(i).row(p) + unobserved_drift_values * dt + t_unobserved_diffusion_values.t() * random_vec;
    }
  }
  return X;
}


// TODO: call simulate, update_weights, resample in a signle gateway function.
// TODO: standardize the function names and arguments.
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
  
  arma::mat weights(nsim, n + 1);
  weights.col(0) = weight_init;

  double t = t0;
  
  // simulate the unobserved process
  arma::cube X = euler_multi_particles(
    x0s, t0, dt, 1, n, r, dW,
    modeltime, modelstate,
    unobserved_drift, unobserved_diffusion,
    env);

  // update weights through time
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
      NumericVector vec_observed_diffusion_values_sexp =
          Rf_eval(observed_diffusion[0], rho);
      arma::mat t_observed_diffusion_values(
          vec_observed_diffusion_values_sexp.begin(), r, d_obs, false);

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
