#include <RcppArmadillo.h>

#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

struct WeightUpdate {
  arma::vec weights;
  double log_normalizer;
};

void assign_state(Environment eval_env,
                  const CharacterVector& unobserved_vars,
                  const arma::vec& unobserved_values,
                  const CharacterVector& observed_vars,
                  const arma::vec& observed_values) {
  for (int i = 0; i < unobserved_vars.size(); ++i) {
    eval_env.assign(as<std::string>(unobserved_vars[i]), unobserved_values(i));
  }
  for (int i = 0; i < observed_vars.size(); ++i) {
    eval_env.assign(as<std::string>(observed_vars[i]), observed_values(i));
  }
}

arma::vec evaluate_vector(const ExpressionVector& expression,
                          Environment eval_env,
                          int expected_size,
                          const char* coefficient_name) {
  NumericVector values = Rf_eval(expression[0], eval_env);
  if (values.size() != expected_size) {
    stop("%s evaluated to length %d; expected length %d.", coefficient_name,
         values.size(), expected_size);
  }
  arma::vec result(values.begin(), values.size(), true);
  if (!result.is_finite()) {
    stop("%s evaluated to non-finite values.", coefficient_name);
  }
  return result;
}

arma::mat evaluate_diffusion(const ExpressionVector& expression,
                             Environment eval_env,
                             int state_dimension,
                             int noise_dimension,
                             const char* coefficient_name) {
  NumericVector values = Rf_eval(expression[0], eval_env);
  const int expected_size = state_dimension * noise_dimension;
  if (values.size() != expected_size) {
    stop("%s evaluated to length %d; expected length %d.", coefficient_name,
         values.size(), expected_size);
  }

  // R supplies coefficients one equation-row at a time. Constructing the
  // transpose in column-major order recovers a state-by-noise matrix.
  arma::mat transposed(values.begin(), noise_dimension, state_dimension, false);
  arma::mat result = transposed.t();
  if (!result.is_finite()) {
    stop("%s evaluated to non-finite values.", coefficient_name);
  }
  return result;
}

arma::cube euler_multi_particles(
    const arma::mat& initial_particles,
    double t0,
    double dt,
    int steps,
    int noise_dimension,
    const std::string& time_var,
    const CharacterVector& unobserved_vars,
    const CharacterVector& observed_vars,
    const arma::vec& observed_values,
    const ExpressionVector& unobserved_drift,
    const ExpressionVector& unobserved_diffusion,
    Environment eval_env,
    const arma::cube& noises) {
  const int num_particles = initial_particles.n_rows;
  const int unobserved_dimension = initial_particles.n_cols;

  arma::cube values(num_particles, unobserved_dimension, steps + 1);
  values.slice(0) = initial_particles;

  for (int i = 0; i < steps; ++i) {
    eval_env.assign(time_var, t0 + i * dt);
    for (int particle = 0; particle < num_particles; ++particle) {
      arma::vec current_state = values.slice(i).row(particle).t();
      assign_state(eval_env, unobserved_vars, current_state, observed_vars,
                   observed_values);

      arma::vec drift = evaluate_vector(
          unobserved_drift, eval_env, unobserved_dimension,
          "Unobserved drift");
      arma::mat diffusion = evaluate_diffusion(
          unobserved_diffusion, eval_env, unobserved_dimension,
          noise_dimension, "Unobserved diffusion");

      arma::vec next_state = current_state + drift * dt +
                             diffusion * noises.slice(i).col(particle);
      if (!next_state.is_finite()) {
        stop("Particle propagation produced non-finite values.");
      }
      values.slice(i + 1).row(particle) = next_state.t();
    }
  }
  return values;
}

WeightUpdate update_weights(
    const arma::cube& values,
    int noise_dimension,
    int observed_dimension,
    const arma::vec& initial_weights,
    const std::string& time_var,
    double time,
    double dt,
    const CharacterVector& unobserved_vars,
    const CharacterVector& observed_vars,
    const arma::vec& observed_values,
    const ExpressionVector& observed_drift,
    const ExpressionVector& observed_diffusion,
    const arma::vec& delta_y,
    Environment eval_env) {
  const int num_particles = initial_weights.n_rows;
  arma::vec log_weights(num_particles);
  log_weights.fill(-std::numeric_limits<double>::infinity());

  eval_env.assign(time_var, time);
  for (int particle = 0; particle < num_particles; ++particle) {
    if (initial_weights(particle) <= 0.0) {
      continue;
    }

    // Section 9.6 uses the left endpoint for the Ito/Euler weight update.
    arma::vec current_state = values.slice(0).row(particle).t();
    assign_state(eval_env, unobserved_vars, current_state, observed_vars,
                 observed_values);

    arma::vec drift = evaluate_vector(
        observed_drift, eval_env, observed_dimension, "Observed drift");
    arma::mat diffusion = evaluate_diffusion(
        observed_diffusion, eval_env, observed_dimension, noise_dimension,
        "Observed diffusion");
    arma::mat covariance = diffusion * diffusion.t();
    arma::mat inverse_covariance;
    if (!covariance.is_finite() ||
        !arma::inv_sympd(inverse_covariance, covariance)) {
      stop("The observation covariance must be finite and positive definite.");
    }

    const arma::vec scaled_delta = inverse_covariance * delta_y;
    const arma::vec scaled_drift = inverse_covariance * drift;
    const double log_increment =
        arma::dot(drift, scaled_delta) -
        0.5 * arma::dot(drift, scaled_drift) * dt;
    log_weights(particle) = std::log(initial_weights(particle)) + log_increment;
  }

  double maximum = -std::numeric_limits<double>::infinity();
  for (int particle = 0; particle < num_particles; ++particle) {
    if (std::isfinite(log_weights(particle)) && log_weights(particle) > maximum) {
      maximum = log_weights(particle);
    }
  }
  if (!std::isfinite(maximum)) {
    stop("All particle weights became non-finite.");
  }

  arma::vec weights(num_particles, arma::fill::zeros);
  for (int particle = 0; particle < num_particles; ++particle) {
    if (std::isfinite(log_weights(particle))) {
      weights(particle) = std::exp(log_weights(particle) - maximum);
    }
  }
  const double total = arma::sum(weights);
  if (!std::isfinite(total) || total <= 0.0) {
    stop("Particle weights could not be normalized.");
  }

  WeightUpdate result;
  result.weights = weights / total;
  result.log_normalizer = maximum + std::log(total);
  return result;
}

arma::vec branch_particles(const arma::vec& weights) {
  const int num_particles = weights.n_elem;
  if (num_particles == 1) {
    return arma::ones<arma::vec>(1);
  }
  if (!weights.is_finite() || arma::any(weights < 0.0)) {
    stop("Branching weights must be finite and non-negative.");
  }

  const double total = arma::sum(weights);
  if (!std::isfinite(total) || total <= 0.0) {
    stop("The sum of branching weights must be positive and finite.");
  }
  const arma::vec normalized_weights = weights / total;

  NumericVector uniforms = Rcpp::runif(num_particles - 1);
  double remaining_mean = static_cast<double>(num_particles);
  int remaining_offspring = num_particles;
  arma::vec offspring(num_particles, arma::fill::zeros);

  // Minimal-variance, fixed-population branching from Section 9.2.1.
  for (int particle = 0; particle < num_particles - 1; ++particle) {
    const double mean_integer = std::floor(remaining_mean);
    const double mean_fraction = remaining_mean - mean_integer;
    const double expected_offspring =
        normalized_weights(particle) * num_particles;
    const double expected_integer = std::floor(expected_offspring);
    const double expected_fraction = expected_offspring - expected_integer;
    const double uniform = uniforms[particle];

    int count;
    if (expected_fraction <= mean_fraction) {
      if (mean_fraction > 0.0 &&
          uniform < expected_fraction / mean_fraction) {
        count = static_cast<int>(expected_integer +
                                 remaining_offspring - mean_integer);
      } else {
        count = static_cast<int>(expected_integer);
      }
    } else {
      const double probability_special =
          (1.0 - expected_fraction) / (1.0 - mean_fraction);
      if (uniform < probability_special) {
        count = static_cast<int>(expected_integer +
                                 remaining_offspring - mean_integer);
      } else {
        count = static_cast<int>(expected_integer + 1.0);
      }
    }

    if (count < 0 || count > remaining_offspring) {
      stop("The branching algorithm produced an invalid offspring count.");
    }
    offspring(particle) = count;
    remaining_mean -= expected_offspring;
    remaining_offspring -= count;
  }
  offspring(num_particles - 1) = remaining_offspring;
  return offspring;
}

double weighted_quantile(const arma::vec& values,
                         const arma::vec& weights,
                         double probability) {
  const arma::uvec order = arma::sort_index(values);
  double cumulative_weight = 0.0;
  for (arma::uword i = 0; i < order.n_elem; ++i) {
    cumulative_weight += weights(order(i));
    if (cumulative_weight >= probability) {
      return values(order(i));
    }
  }
  return values(order(order.n_elem - 1));
}

void summarize_particles(const arma::mat& particles,
                         const arma::vec& weights,
                         int time_index,
                         const arma::vec& interval_levels,
                         arma::mat& means,
                         arma::cube& covariances,
                         arma::vec& ess,
                         NumericVector& intervals) {
  const int state_dimension = particles.n_cols;
  const int num_times = means.n_cols;

  const arma::vec current_mean = particles.t() * weights;
  means.col(time_index) = current_mean;
  const arma::mat centered = particles.each_row() - current_mean.t();
  covariances.slice(time_index) =
      centered.t() * (centered.each_col() % weights);
  ess(time_index) = 1.0 / arma::dot(weights, weights);

  for (arma::uword level_index = 0;
       level_index < interval_levels.n_elem; ++level_index) {
    const double lower_probability = (1.0 - interval_levels(level_index)) / 2.0;
    const double upper_probability = 1.0 - lower_probability;
    for (int state = 0; state < state_dimension; ++state) {
      const arma::vec state_values = particles.col(state);
      const R_xlen_t lower_index =
          state + state_dimension *
                      (0 + 2 * (time_index + num_times * level_index));
      const R_xlen_t upper_index =
          state + state_dimension *
                      (1 + 2 * (time_index + num_times * level_index));
      intervals[lower_index] = weighted_quantile(
          state_values, weights, lower_probability);
      intervals[upper_index] = weighted_quantile(
          state_values, weights, upper_probability);
    }
  }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List euler_multi_particles_with_weights_and_branching(
    const arma::mat& xinits,
    const arma::vec& weight_init,
    double t0,
    int noise_dimension,
    int observed_dimension,
    double dt,
    int steps,
    const std::string& time_var,
    const CharacterVector& unobserved_vars,
    const CharacterVector& observed_vars,
    int simulations_per_weight_update,
    int weight_updates_per_branching,
    const ExpressionVector& observed_drift,
    const ExpressionVector& unobserved_drift,
    const ExpressionVector& observed_diffusion,
    const ExpressionVector& unobserved_diffusion,
    const arma::mat& delta_y,
    const arma::mat& observed_values,
    Environment eval_env,
    int particle_storage,
    bool keep_ancestors,
    const arma::vec& interval_levels) {
  if (steps <= 0 || simulations_per_weight_update <= 0 ||
      weight_updates_per_branching <= 0) {
    stop("Step and branching counts must be positive.");
  }
  if (noise_dimension <= 0 || observed_dimension <= 0) {
    stop("Noise and observation dimensions must be positive.");
  }
  if (particle_storage < 0 || particle_storage > 3) {
    stop("Unknown particle storage mode.");
  }
  if (xinits.n_rows == 0 || xinits.n_cols == 0 ||
      weight_init.n_elem != xinits.n_rows) {
    stop("Initial particles and weights have incompatible dimensions.");
  }
  if (unobserved_vars.size() != static_cast<int>(xinits.n_cols) ||
      observed_vars.size() != observed_dimension) {
    stop("State-variable names have incompatible dimensions.");
  }
  if (delta_y.n_rows != static_cast<arma::uword>(observed_dimension) ||
      delta_y.n_cols != static_cast<arma::uword>(steps) ||
      observed_values.n_rows != static_cast<arma::uword>(observed_dimension) ||
      observed_values.n_cols != static_cast<arma::uword>(steps)) {
    stop("Observation matrices have incompatible dimensions.");
  }
  if (!xinits.is_finite() || !weight_init.is_finite() ||
      !delta_y.is_finite() || !observed_values.is_finite() ||
      arma::any(weight_init < 0.0) || arma::sum(weight_init) <= 0.0) {
    stop("Particle-filter inputs must be finite and weights non-negative.");
  }
  if (!interval_levels.is_finite() || arma::any(interval_levels <= 0.0) ||
      arma::any(interval_levels >= 1.0)) {
    stop("Interval levels must lie strictly between zero and one.");
  }

  const int num_particles = xinits.n_rows;
  const int state_dimension = xinits.n_cols;
  const int num_times = steps + 1;
  const int total_simulation_steps = steps * simulations_per_weight_update;
  arma::vec current_weights = weight_init / arma::sum(weight_init);
  arma::mat current_particles = xinits;

  arma::mat means(state_dimension, num_times);
  arma::cube covariances(state_dimension, state_dimension, num_times);
  arma::vec ess(num_times);
  arma::vec log_likelihood_increment(num_times, arma::fill::zeros);
  LogicalVector branched(num_times, false);
  NumericVector intervals(
      static_cast<R_xlen_t>(state_dimension) * 2 * num_times *
      interval_levels.n_elem);
  intervals.attr("dim") = IntegerVector::create(
      state_dimension, 2, num_times, interval_levels.n_elem);

  summarize_particles(current_particles, current_weights, 0, interval_levels,
                      means, covariances, ess, intervals);

  const bool store_observations = particle_storage >= 2;
  const bool store_last = particle_storage == 1;
  const bool store_paths = particle_storage == 3;
  arma::cube stored_particles;
  arma::mat stored_weights;
  if (store_observations) {
    stored_particles.set_size(num_particles, state_dimension, num_times);
    stored_weights.set_size(num_particles, num_times);
    stored_particles.slice(0) = current_particles;
    stored_weights.col(0) = current_weights;
  } else if (store_last) {
    stored_particles.set_size(num_particles, state_dimension, 1);
    stored_weights.set_size(num_particles, 1);
    stored_particles.slice(0) = current_particles;
    stored_weights.col(0) = current_weights;
  }

  arma::cube paths;
  if (store_paths) {
    paths.set_size(num_particles, state_dimension,
                   total_simulation_steps + 1);
    paths.slice(0) = current_particles;
  }

  arma::umat ancestors;
  if (keep_ancestors) {
    ancestors.set_size(num_particles, num_times);
    for (int particle = 0; particle < num_particles; ++particle) {
      ancestors(particle, 0) = particle + 1;
    }
  }

  for (int i = 0; i < steps; ++i) {
    arma::cube noises(noise_dimension, num_particles,
                      simulations_per_weight_update);
    noises.randn();
    noises *= std::sqrt(dt / simulations_per_weight_update);

    arma::cube values = euler_multi_particles(
        current_particles, t0 + i * dt,
        dt / simulations_per_weight_update, simulations_per_weight_update,
        noise_dimension, time_var, unobserved_vars, observed_vars,
        observed_values.col(i), unobserved_drift, unobserved_diffusion,
        eval_env, noises);
    const arma::mat endpoint = values.slice(simulations_per_weight_update);

    if (store_paths) {
      for (int substep = 1; substep <= simulations_per_weight_update;
           ++substep) {
        paths.slice(i * simulations_per_weight_update + substep) =
            values.slice(substep);
      }
    }

    const WeightUpdate update = update_weights(
        values, noise_dimension, observed_dimension, current_weights,
        time_var, t0 + i * dt, dt, unobserved_vars, observed_vars,
        observed_values.col(i), observed_drift, observed_diffusion,
        delta_y.col(i), eval_env);
    const arma::vec new_weights = update.weights;
    log_likelihood_increment(i + 1) = update.log_normalizer;
    summarize_particles(endpoint, new_weights, i + 1, interval_levels,
                        means, covariances, ess, intervals);

    if (store_observations) {
      stored_particles.slice(i + 1) = endpoint;
      stored_weights.col(i + 1) = new_weights;
    } else if (store_last) {
      stored_particles.slice(0) = endpoint;
      stored_weights.col(0) = new_weights;
    }

    if ((i + 1) % weight_updates_per_branching == 0) {
      branched[i + 1] = true;
      const arma::vec offspring = branch_particles(new_weights);
      arma::mat branched_particles(num_particles, state_dimension);
      int next_particle = 0;
      for (int parent = 0; parent < num_particles; ++parent) {
        const int count = static_cast<int>(offspring(parent));
        for (int child = 0; child < count; ++child) {
          branched_particles.row(next_particle) = endpoint.row(parent);
          if (keep_ancestors) {
            ancestors(next_particle, i + 1) = parent + 1;
          }
          ++next_particle;
        }
      }
      if (next_particle != num_particles) {
        stop("Branching did not preserve the number of particles.");
      }
      current_particles = branched_particles;
      current_weights.fill(1.0 / num_particles);
    } else {
      current_particles = endpoint;
      current_weights = new_weights;
      if (keep_ancestors) {
        for (int particle = 0; particle < num_particles; ++particle) {
          ancestors(particle, i + 1) = particle + 1;
        }
      }
    }
  }

  RObject particle_result = R_NilValue;
  RObject weight_result = R_NilValue;
  RObject path_result = R_NilValue;
  RObject ancestor_result = R_NilValue;
  if (particle_storage != 0) {
    particle_result = wrap(stored_particles);
    weight_result = wrap(stored_weights);
  }
  if (store_paths) {
    path_result = wrap(paths);
  }
  if (keep_ancestors) {
    ancestor_result = wrap(ancestors);
  }

  return Rcpp::List::create(
      Rcpp::Named("particles") = particle_result,
      Rcpp::Named("weights") = weight_result,
      Rcpp::Named("paths") = path_result,
      Rcpp::Named("ancestors") = ancestor_result,
      Rcpp::Named("mean") = means,
      Rcpp::Named("vcov") = covariances,
      Rcpp::Named("ess") = ess,
      Rcpp::Named("branched") = branched,
      Rcpp::Named("intervals") = intervals,
      Rcpp::Named("logLik") = arma::sum(log_likelihood_increment),
      Rcpp::Named("logLik_increment") = log_likelihood_increment);
}
