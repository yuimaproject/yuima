#include <RcppArmadillo.h>
// #include <time.h>

using namespace Rcpp;

arma::cube euler_multi_particles(arma::mat xinits, double t0, double dt,
                                 int steps, int d_random, std::string time_var,
                                 CharacterVector unobserved_vars,
                                 ExpressionVector unobserved_drift,
                                 ExpressionVector unobserved_diffusion,
                                 Environment eval_env) {
  int num_particles = xinits.n_rows;
  int d_unob = xinits.n_cols;

  arma::cube X(num_particles, d_unob, steps + 1);
  X.slice(0) = xinits;

  for (int i = 0; i < steps; i++) {
    eval_env.assign(time_var, t0 + i * dt);
    for (int p = 0; p < num_particles; p++) {
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
          vec_unobserved_diffusion_values_sexp.begin(), d_random, d_unob,
          false);  // TODO: consider transpose

      // simulate next step
      arma::vec random_vec = Rcpp::rnorm(d_random, 0, sqrt(dt));
      X.slice(i + 1).row(p) = X.slice(i).row(p) + unobserved_drift_values * dt +
                              t_unobserved_diffusion_values.t() * random_vec;
    }
  }
  return X;
}

arma::mat update_weights(arma::cube X, int r, arma::vec weights_init,
                         std::string time_var, double t0, double dt,
                         CharacterVector unobserved_vars,
                         int simulations_per_weight_update, int steps,
                         ExpressionVector observed_drift,
                         ExpressionVector observed_diffusion, arma::mat deltaY,
                         Environment eval_env) {
  int d_ob = observed_drift.size();  // number of observed dimensions
  int d_unob = X.n_cols;
  arma::mat weights(weights_init.n_rows, steps + 1);
  weights.col(0) = weights_init;

  for (int i = 0; i < steps; i++) {
    eval_env.assign(time_var, t0 + i * dt);
    for (int p = 0; p < weights_init.n_rows; p++) {
      // if the weight is zero, set the next weight to zero
      if (weights(p, i) == 0) {
        weights(p, i + 1) = 0;
        continue;
      }

      // assign the current state for each dimension to the environment
      for (int v = 0; v < d_unob; v++) {
        eval_env.assign(as<std::string>(unobserved_vars[v]),
                        X(p, v, i * simulations_per_weight_update));
      }

      // Evaluate the drift and diffusion expressions
      arma::vec observed_drift_values =
          as<arma::vec>(Rf_eval(observed_drift[0], eval_env));
      NumericVector vec_observed_diffusion_values_sexp =
          Rf_eval(observed_diffusion[0], eval_env);
      arma::mat t_observed_diffusion_values(
          vec_observed_diffusion_values_sexp.begin(), r, d_ob, false);

      arma::mat inv_sigma = arma::inv_sympd(t_observed_diffusion_values.t() *
                                            t_observed_diffusion_values);
      // Update weights
      double coeff = std::exp(arma::as_scalar(
          observed_drift_values.t() * inv_sigma * deltaY.col(i) -
          0.5 * observed_drift_values.t() * inv_sigma *
              observed_drift_values.t() * dt));
      weights(p, i + 1) = weights(p, i) * coeff;

      // clip negative weights to zero
      if (weights(p, i + 1) < 0) {
        weights(p, i + 1) = 0;
      }
    }
  }

  return weights;
}

// TODO: call simulate, update_weights, branch in a signle gateway function.
// TODO: standardize the function names and arguments.
// [[Rcpp::export]]
Rcpp::List euler_multi_particles_with_weights(
    arma::mat xinits, arma::vec weight_init, double t0, int r, double dt,
    int steps, std::string time_var, CharacterVector unobserved_vars,
    int simulations_per_weight_update, ExpressionVector observed_drift,
    ExpressionVector unobserved_drift, ExpressionVector observed_diffusion,
    ExpressionVector unobserved_diffusion, arma::mat deltaY,
    Environment eval_env) {
  // simulate the unobserved process
  arma::cube X = euler_multi_particles(
      xinits, t0, dt / simulations_per_weight_update,
      steps * simulations_per_weight_update, r, time_var, unobserved_vars,
      unobserved_drift, unobserved_diffusion, eval_env);

  // update weights through time
  arma::mat weights =
      update_weights(X, r, weight_init, time_var, t0, dt, unobserved_vars,
                     simulations_per_weight_update, steps, observed_drift,
                     observed_diffusion, deltaY, eval_env);

  return Rcpp::List::create(Rcpp::Named("values") = X,
                            Rcpp::Named("weights") = weights);
}

// [[Rcpp::export]]
arma::vec branch_particles(arma::vec weights) {
  const int num_particles = weights.n_elem;
  if (num_particles <= 1) {
    warning("number of particles should be greater than 1");
    return arma::ones(num_particles);
  }

  // normalize weights to make the sum of weights equal to 1
  double sum_weights = arma::sum(weights);
  if (sum_weights == 0) {
    warning("sum of weights is zero");
    return arma::ones(num_particles);
  }
  weights /= sum_weights;

  arma::vec us = Rcpp::runif(num_particles);
  double g = static_cast<double>(num_particles);
  int h = num_particles;
  arma::vec numbers(num_particles);

  for (int i = 0; i < num_particles - 1; ++i) {
    double g_int = std::floor(g);
    double g_dec = g - g_int;
    double nw = weights(i) * num_particles;
    double nw_int = std::floor(nw);
    double nw_dec = nw - nw_int;
    double u = us(i);
    int o;
    if (nw_dec <= g_dec) {
      if (g_dec > 0.0 && u <= nw_dec / g_dec)
        o = static_cast<int>(nw_int + h - g_int);
      else
        o = static_cast<int>(nw_int);
    } else {
      if ((1.0 - g_dec) > 0.0 && u <= (1.0 - nw_dec) / (1.0 - g_dec))
        o = static_cast<int>(nw_int + h - g_int);
      else
        o = static_cast<int>(nw_int + 1);
    }
    numbers(i) = o;
    g -= nw;
    h -= o;
  }
  numbers(num_particles - 1) = h;
  return numbers;
}

// [[Rcpp::export]]
Rcpp::List euler_multi_particles_with_weights_and_branching(
    arma::mat xinits, arma::vec weight_init, double t0, int r, double dt,
    int steps, std::string time_var, CharacterVector unobserved_vars,
    int simulations_per_weight_update, int weight_update_per_branching,
    ExpressionVector observed_drift, ExpressionVector unobserved_drift,
    ExpressionVector observed_diffusion, ExpressionVector unobserved_diffusion,
    arma::mat deltaY, Environment eval_env) {
  arma::cube all_values(xinits.n_rows, xinits.n_cols,
                        steps * simulations_per_weight_update + 1);
  all_values.slice(0) = xinits;
  arma::mat all_weights(weight_init.n_rows, steps + 1);
  all_weights.col(0) = weight_init;

  int num_particles = all_weights.n_rows;

  for (int i = 0; i < steps; i++) {
    // simulate the unobserved process
    arma::cube X = euler_multi_particles(
        all_values.slice(i * simulations_per_weight_update), t0 + i * dt,
        dt / simulations_per_weight_update, simulations_per_weight_update, r,
        time_var, unobserved_vars, unobserved_drift, unobserved_diffusion,
        eval_env);
    // std::cout << "debug: " << X.n_slices << "==" <<
    // simulations_per_weight_update + 1 << "?" << std::endl;
    all_values.slices(i * simulations_per_weight_update,
                      (i + 1) * simulations_per_weight_update) = X;

    // update weights through time
    arma::mat weights = update_weights(
        X, r, all_weights.col(i), time_var, t0 + i * dt, dt, unobserved_vars,
        simulations_per_weight_update, 1, observed_drift, observed_diffusion,
        deltaY.cols(i, i), eval_env);
    all_weights.col(i + 1) = weights.col(1);
    // branch particles
    if (i % weight_update_per_branching == weight_update_per_branching - 1) {
      arma::vec numbers = branch_particles(all_weights.col(i + 1));
      int current_num_particles = 0;
      for (int j = 0; j < num_particles; ++j) {
        int o = static_cast<int>(numbers(j));
        if (o == 0) {
          continue;
        }
        for (int k = 0; k < o; ++k) {
          all_values.slice((i + 1) * simulations_per_weight_update)
              .row(current_num_particles) =
              all_values.slice(i * simulations_per_weight_update).row(j);
          all_weights.col(i + 1).row(current_num_particles) =
              1.0 / num_particles;
          current_num_particles++;
        }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("values") = all_values,
                            Rcpp::Named("weights") = all_weights);
}