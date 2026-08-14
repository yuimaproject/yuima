make_particle_filter_yuima <- function(model, sampling, xinit, parameters) {
  simulated <- suppressWarnings(
    simulate(
      model,
      sampling = sampling,
      xinit = xinit,
      true.parameter = parameters
    )
  )
  variables <- model@state.variable
  mapping <- as.list(seq_along(variables))
  names(mapping) <- variables
  suppressWarnings(
    setYuima(
      model = model,
      data = simulated@data,
      sampling = sampling,
      variable_data_mapping = mapping
    )
  )
}

make_gaussian_particle_filter_case <- function(
    n = 30L,
    particles = 2000L,
    interval_levels = c(0.9, 0.95)) {
  model <- setModel(
    drift = c("-a * X", "c * X"),
    diffusion = matrix(
      c("b", "0", "0", "obs_sd"),
      nrow = 2L,
      byrow = TRUE
    ),
    solve.variable = c("X", "Y"),
    state.variable = c("X", "Y"),
    observed.variable = "Y"
  )
  sampling <- suppressWarnings(setSampling(n = n, Terminal = n * 0.01))
  parameters <- list(a = 0.7, b = 0.3, c = 1, obs_sd = 0.4)
  set.seed(1001)
  object <- make_particle_filter_yuima(
    model, sampling, c(0, 0), parameters
  )
  initial_particles <- matrix(
    rnorm(particles, mean = 0, sd = 0.5),
    ncol = 1L
  )
  result <- particleFilter(
    object,
    xinits = initial_particles,
    params = parameters,
    simulations_per_weight_update = 2L,
    weight_updates_per_branching = 10L,
    store_particles = "observation",
    interval_levels = interval_levels,
    seed = 1002
  )
  list(
    model = model,
    object = object,
    parameters = parameters,
    result = result
  )
}

make_lorenz_particle_filter_case <- function(
    n = 100L,
    particles = 600L,
    store_particles = "observation") {
  model <- setModel(
    drift = c(
      "sigma * (X2 - X1)",
      "X1 * (rho - X3) - X2",
      "X1 * X2 - beta * X3",
      "X1",
      "X3"
    ),
    diffusion = diag(c(0.10, 0.10, 0.10, 0.60, 0.60)),
    solve.variable = c("X1", "X2", "X3", "Y1", "Y2"),
    state.variable = c("X1", "X2", "X3", "Y1", "Y2"),
    observed.variable = c("Y1", "Y2")
  )
  sampling <- suppressWarnings(setSampling(delta = 0.01, n = n))
  parameters <- list(sigma = 10, rho = 28, beta = 8 / 3)
  set.seed(3001)
  object <- make_particle_filter_yuima(
    model, sampling, c(1, 1, 1, 0, 0), parameters
  )
  initial_particles <- cbind(
    rnorm(particles, 1, 0.2),
    rnorm(particles, 1, 0.2),
    rnorm(particles, 1, 0.2)
  )
  result <- particleFilter(
    object,
    xinits = initial_particles,
    params = parameters,
    simulations_per_weight_update = 1L,
    weight_updates_per_branching = 8L,
    store_particles = store_particles,
    interval_levels = 0.9,
    seed = 3002
  )
  truth <- vapply(
    seq_len(3L),
    function(i) as.numeric(object@data@zoo.data[[i]]),
    numeric(n + 1L)
  )
  lower <- t(result@intervals[, 1L, , 1L])
  upper <- t(result@intervals[, 2L, , 1L])

  list(
    model = model,
    object = object,
    parameters = parameters,
    result = result,
    truth = truth,
    rmse = sqrt(colMeans((as.matrix(result@mean) - truth)^2)),
    coverage = colMeans(truth >= lower & truth <= upper)
  )
}

test_that("Gaussian particle filtering agrees with Kalman-Bucy filtering", {
  case <- make_gaussian_particle_filter_case()
  kalman <- kalmanBucyFilter(
    case$object,
    params = case$parameters,
    mean_init = 0,
    vcov_init = matrix(0.25)
  )

  particle_mean <- as.numeric(mean(case$result))
  kalman_mean <- as.numeric(mean(kalman))
  expect_equal(length(particle_mean), length(kalman_mean))
  expect_lt(sqrt(mean((particle_mean - kalman_mean)^2)), 0.05)
  expect_lt(abs(tail(particle_mean, 1L) - tail(kalman_mean, 1L)), 0.05)
  expect_true(all(is.finite(case$result@vcov)))
})

test_that("a nonlinear multidimensional model has stable filter summaries", {
  model <- setModel(
    drift = c(
      "-0.8 * X1 + 0.3 * X2",
      "-0.4 * X2 - 0.15 * X1^3",
      "sin(X1) + 0.2 * X2",
      "X1 * X2 / (1 + X1^2)"
    ),
    diffusion = matrix(c(
      0.25, 0, 0, 0,
      0, 0.20, 0, 0,
      0, 0, 0.35, 0.05,
      0, 0, 0.02, 0.30
    ), nrow = 4L, byrow = TRUE),
    solve.variable = c("X1", "X2", "Y1", "Y2"),
    state.variable = c("X1", "X2", "Y1", "Y2"),
    observed.variable = c("Y1", "Y2")
  )
  sampling <- suppressWarnings(setSampling(n = 25L, Terminal = 0.25))
  set.seed(2001)
  object <- make_particle_filter_yuima(
    model, sampling, c(0.2, -0.1, 0, 0), list()
  )
  initial_particles <- cbind(
    rnorm(500L, 0.2, 0.25),
    rnorm(500L, -0.1, 0.25)
  )
  result <- particleFilter(
    object,
    xinits = initial_particles,
    simulations_per_weight_update = 2L,
    weight_updates_per_branching = 5L,
    store_particles = "observation",
    interval_levels = c(0.5, 0.9),
    seed = 2002
  )

  expect_equal(dim(result@mean), c(26L, 2L))
  expect_equal(dim(result@vcov), c(2L, 2L, 26L))
  expect_equal(dim(result@intervals), c(2L, 2L, 26L, 2L))
  expect_true(all(is.finite(result@mean)))
  expect_true(all(is.finite(result@vcov)))
  expect_true(all(result@intervals[, 1L, , ] <= result@intervals[, 2L, , ]))
  expect_equal(unname(colSums(result@weights)), rep(1, 26L), tolerance = 1e-12)
  expect_true(all(result@ess > 0 & result@ess <= 500L + 1e-10))
})

test_that("Lorenz-63 latent dynamics can be filtered from two observations", {
  case <- make_lorenz_particle_filter_case()
  result <- case$result

  expect_equal(dim(result@mean), c(101L, 3L))
  expect_equal(dim(result@particles), c(600L, 3L, 101L))
  expect_true(all(is.finite(result@mean)))
  expect_true(all(is.finite(result@vcov)))
  expect_true(is.finite(result@logLik))
  expect_equal(sum(result@branched), 12L)
  expect_true(all(result@ess > 0 & result@ess <= 600L + 1e-10))
  expect_lt(max(case$rmse), 2)
  expect_gt(min(case$coverage), 0.6)

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  expect_invisible(plot(result, plot_truth = TRUE, level = 0.9))
  grDevices::dev.off()
  expect_gt(file.info(plot_file)$size, 0)
  unlink(plot_file)
})

test_that("Lorenz-63 long-horizon particle filtering remains calibrated", {
  skip_on_cran()
  if (!identical(Sys.getenv("YUIMA_RUN_LONG_TESTS"), "true")) {
    skip("Set YUIMA_RUN_LONG_TESTS=true to run long particle-filter tests.")
  }

  case <- make_lorenz_particle_filter_case(
    n = 1000L,
    particles = 1000L,
    store_particles = "none"
  )

  expect_lt(max(case$rmse), 2)
  expect_gt(min(case$coverage), 0.75)
  expect_true(all(case$result@ess > 0))
  expect_equal(sum(case$result@branched), 125L)
})

test_that("particle-filter plotting restores par and validates interval levels", {
  case <- make_gaussian_particle_filter_case(
    n = 10L,
    particles = 400L,
    interval_levels = c(0.9, 0.95)
  )
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  graphics::par(
    mar = c(3, 3, 2, 1),
    oma = c(0, 1, 0, 1),
    mgp = c(2, 0.7, 0),
    mfrow = c(1, 1)
  )
  before <- graphics::par(no.readonly = TRUE)

  expect_invisible(
    plot(case$result, plot_truth = TRUE, level = 0.95)
  )
  after <- graphics::par(no.readonly = TRUE)
  expect_equal(after, before, tolerance = 1e-14)

  before_error <- graphics::par(no.readonly = TRUE)
  expect_error(
    plot(case$result, level = 0.8),
    "requested interval level 0.8 was not stored"
  )
  expect_equal(
    graphics::par(no.readonly = TRUE),
    before_error,
    tolerance = 1e-14
  )

  missing_truth <- case$result
  missing_truth@data@zoo.data[[1L]][] <- NA_real_
  before_drawing_error <- graphics::par(no.readonly = TRUE)
  expect_error(
    plot(missing_truth, plot_truth = TRUE, level = 0.95),
    "No data for X is found"
  )
  expect_equal(
    graphics::par(no.readonly = TRUE),
    before_drawing_error,
    tolerance = 1e-14
  )

  grDevices::dev.off()
  expect_gt(file.info(plot_file)$size, 0)
  unlink(plot_file)
})
