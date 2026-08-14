particle_filter_cpp <- function(
    xinits,
    weight_init,
    t0 = 0,
    noise_dimension = 1L,
    observed_dimension = 1L,
    dt = 0.1,
    steps = 1L,
    time_var = "t",
    unobserved_vars = "X",
    observed_vars = "Y",
    simulations_per_weight_update = 1L,
    weight_updates_per_branching = 2L,
    observed_drift = parse(text = "c(0)"),
    unobserved_drift = parse(text = "c(0)"),
    observed_diffusion = parse(text = "c(1)"),
    unobserved_diffusion = parse(text = "c(0)"),
    delta_y = matrix(0, nrow = observed_dimension, ncol = steps),
    observed_values = matrix(0, nrow = observed_dimension, ncol = steps),
    particle_storage = 2L,
    keep_ancestors = FALSE,
    interval_levels = numeric()) {
  euler_multi_particles_with_weights_and_branching(
    xinits = xinits,
    weight_init = weight_init,
    t0 = t0,
    noise_dimension = noise_dimension,
    observed_dimension = observed_dimension,
    dt = dt,
    steps = steps,
    time_var = time_var,
    unobserved_vars = unobserved_vars,
    observed_vars = observed_vars,
    simulations_per_weight_update = simulations_per_weight_update,
    weight_updates_per_branching = weight_updates_per_branching,
    observed_drift = observed_drift,
    unobserved_drift = unobserved_drift,
    observed_diffusion = observed_diffusion,
    unobserved_diffusion = unobserved_diffusion,
    delta_y = delta_y,
    observed_values = observed_values,
    eval_env = new.env(parent = baseenv()),
    particle_storage = particle_storage,
    keep_ancestors = keep_ancestors,
    interval_levels = interval_levels
  )
}

make_multidimensional_particle_filter_object <- function() {
  model <- setModel(
    drift = c(
      "-0.7 * X1",
      "-0.4 * X2",
      "X1 + 0.5 * X2",
      "-0.25 * X1 + X2"
    ),
    diffusion = matrix(c(
      0.35, 0, 0, 0,
      0, 0.25, 0, 0,
      0, 0, 0.30, 0.10,
      0, 0, 0.05, 0.40
    ), nrow = 4, byrow = TRUE),
    solve.variable = c("X1", "X2", "Y1", "Y2"),
    state.variable = c("X1", "X2", "Y1", "Y2"),
    observed.variable = c("Y1", "Y2")
  )
  sampling <- suppressWarnings(setSampling(n = 20, Terminal = 0.2))
  set.seed(42)
  simulated <- suppressWarnings(
    simulate(model, sampling = sampling, xinit = rep(0, 4))
  )
  suppressWarnings(
    setYuima(
      model = model,
      data = simulated@data,
      sampling = sampling,
      variable_data_mapping = list(X1 = 1, X2 = 2, Y1 = 3, Y2 = 4)
    )
  )
}

test_that("multidimensional log weights follow Section 9.6", {
  result <- particle_filter_cpp(
    xinits = matrix(c(0, 1), ncol = 1),
    weight_init = c(0.5, 0.5),
    noise_dimension = 2L,
    observed_dimension = 2L,
    observed_vars = c("Y1", "Y2"),
    observed_drift = parse(text = "c(X, 2 * X)"),
    observed_diffusion = parse(text = "c(1, 0, 0, 1)"),
    unobserved_diffusion = parse(text = "c(0, 0)"),
    delta_y = matrix(c(1, 1), nrow = 2),
    observed_values = matrix(c(0, 0), nrow = 2)
  )

  expected <- exp(c(0, 2.75) - 2.75)
  expected <- expected / sum(expected)
  expect_equal(result$weights[, 2], expected, tolerance = 1e-12)
  expect_equal(sum(result$weights[, 2]), 1, tolerance = 1e-14)
})

test_that("the weight update uses the interval start time", {
  result <- particle_filter_cpp(
    xinits = matrix(c(0, 1), ncol = 1),
    weight_init = c(0.5, 0.5),
    t0 = 2,
    observed_drift = parse(text = "c(t * X)"),
    delta_y = matrix(1, nrow = 1)
  )

  expected <- exp(c(0, 1.8) - 1.8)
  expected <- expected / sum(expected)
  expect_equal(result$weights[, 2], expected, tolerance = 1e-12)
})

test_that("branching copies the propagated endpoint", {
  result <- particle_filter_cpp(
    xinits = matrix(0, nrow = 20, ncol = 1),
    weight_init = rep(1 / 20, 20),
    weight_updates_per_branching = 1L,
    unobserved_drift = parse(text = "c(1)"),
    particle_storage = 3L,
    keep_ancestors = TRUE
  )

  expect_equal(result$particles[, 1, 2], rep(0.1, 20), tolerance = 1e-14)
  expect_equal(result$paths[, 1, 2], rep(0.1, 20), tolerance = 1e-14)
  expect_equal(result$weights[, 2], rep(1 / 20, 20), tolerance = 1e-14)
  expect_true(result$branched[[2]])
  expect_equal(result$ancestors[, 2], seq_len(20))
})

test_that("weighted summaries and intervals are computed before branching", {
  result <- particle_filter_cpp(
    xinits = matrix(c(0, 1), ncol = 1),
    weight_init = c(0.25, 0.75),
    interval_levels = 0.5,
    particle_storage = 0L
  )

  expect_equal(result$mean[, 1], 0.75)
  expect_equal(as.numeric(result$vcov[, , 1]), 0.1875, tolerance = 1e-14)
  expect_equal(result$ess[[1]], 1 / (0.25^2 + 0.75^2), tolerance = 1e-14)
  expect_equal(result$intervals[1, , 1, 1], c(0, 1))
  expect_null(result$particles)
})

test_that("the public filter supports multidimensional observations", {
  object <- make_multidimensional_particle_filter_object()
  result <- particleFilter(
    object,
    xinits = matrix(0, nrow = 100, ncol = 2),
    simulations_per_weight_update = 2,
    weight_updates_per_branching = 5,
    store_particles = "observation",
    keep_ancestors = TRUE,
    interval_levels = c(0.5, 0.95),
    seed = 123
  )

  expect_s4_class(result, "yuima.particleFilter")
  expect_equal(dim(result@particles), c(100, 2, 21))
  expect_equal(dim(result@weights), c(100, 21))
  expect_equal(dim(result@mean), c(21, 2))
  expect_equal(dim(result@vcov), c(2, 2, 21))
  expect_equal(dim(result@intervals), c(2, 2, 21, 2))
  expect_equal(dim(result@ancestors), c(100, 21))
  expect_equal(unname(colSums(result@weights)), rep(1, 21), tolerance = 1e-12)
  expect_equal(sum(result@branched), 4)
  expect_true(all(is.finite(result@ess)))
  expect_true(is.finite(result@logLik))
  expect_identical(mean(result), result@mean)
  expect_identical(vcov(result), result@vcov)
})

test_that("particle storage modes have stable result shapes", {
  object <- make_multidimensional_particle_filter_object()
  common <- list(
    yuima = object,
    xinits = matrix(0, nrow = 25, ncol = 2),
    steps = 3L,
    simulations_per_weight_update = 2L,
    seed = 10
  )

  last <- do.call(particleFilter, common)
  none <- do.call(particleFilter, c(common, list(store_particles = "none")))
  all <- do.call(particleFilter, c(common, list(store_particles = "all")))

  expect_equal(dim(last@particles), c(25, 2, 1))
  expect_equal(dim(last@weights), c(25, 1))
  expect_equal(dim(last@paths), c(0, 0, 0))
  expect_equal(dim(none@particles), c(0, 0, 0))
  expect_equal(dim(none@weights), c(0, 0))
  expect_equal(dim(all@particles), c(25, 2, 4))
  expect_equal(dim(all@paths), c(25, 2, 7))
  expect_equal(dim(none@mean), c(4, 2))
  expect_equal(dim(none@intervals), c(2, 2, 4, 1))

  no_intervals <- do.call(
    particleFilter,
    c(common, list(store_particles = "none", interval_levels = NULL))
  )
  expect_equal(dim(no_intervals@intervals), c(2, 2, 4, 0))
  expect_length(no_intervals@interval.levels, 0L)
})

test_that("show and summary follow the Kalman-Bucy result style", {
  object <- make_multidimensional_particle_filter_object()
  result <- particleFilter(
    object,
    xinits = matrix(0, nrow = 20, ncol = 2),
    steps = 2L,
    seed = 1
  )

  shown <- capture.output(show(result))
  summarized <- capture.output(summary(result))
  expect_match(shown[[1]], "Particle Filter")
  expect_true(any(grepl("Mean values", shown, fixed = TRUE)))
  expect_true(any(grepl("Variance-covariance matrices", shown, fixed = TRUE)))
  expect_match(summarized[[1]], "Summary of estimation by Particle Filter")
  expect_true(any(grepl("Storage mode: last", summarized, fixed = TRUE)))
  expect_true(any(grepl("Minimum ESS", summarized, fixed = TRUE)))
})

test_that("the compatibility wrapper returns the new result object", {
  object <- make_multidimensional_particle_filter_object()
  result <- simulate_multi_particles_with_weights(
    object,
    xinits = matrix(0, nrow = 10, ncol = 2),
    steps = 1L,
    seed = 2
  )
  expect_s4_class(result, "yuima.particleFilter")
})

test_that("the public filter validates unsupported inputs", {
  correlated_model <- setModel(
    drift = c("-X", "X"),
    diffusion = matrix(c(1, 0.2), nrow = 2),
    solve.variable = c("X", "Y"),
    state.variable = c("X", "Y"),
    observed.variable = "Y"
  )
  sampling <- suppressWarnings(setSampling(n = 2, Terminal = 0.2))
  simulated <- suppressWarnings(
    simulate(correlated_model, sampling = sampling, xinit = c(0, 0))
  )
  object <- suppressWarnings(
    setYuima(
      model = correlated_model,
      data = simulated@data,
      sampling = sampling,
      variable_data_mapping = list(X = 1, Y = 2)
    )
  )

  expect_error(
    particleFilter(
      object,
      xinits = matrix(0, nrow = 10),
      weight_updates_per_branching = 0
    ),
    "positive integer"
  )
  expect_error(
    particleFilter(object, xinits = matrix(0, nrow = 10)),
    "independent"
  )
  expect_error(
    particleFilter(
      make_multidimensional_particle_filter_object(),
      xinits = matrix(0, nrow = 10, ncol = 2),
      store_particles = "last",
      keep_ancestors = TRUE
    ),
    "requires store_particles"
  )
})
