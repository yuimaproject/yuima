make_time_qmle_data <- function(model) {
  values <- c(0.2, 0.1, 0.4, -0.1, 0.3, 0.5)
  suppressWarnings(setYuima(
    model = model,
    data = setData(ts(values, start = 2, deltat = 0.1)),
    variable_data_mapping = list(x = 1)
  ))
}

test_that("pointwise likelihood uses observation times in drift and diffusion", {
  for (clock in c("t", "s")) {
    model <- setModel(
      drift = paste0("-alpha*x + ", clock),
      diffusion = paste0("exp(gamma * ", clock, " / 2)"),
      time.variable = clock, state.variable = "x", solve.variable = "x"
    )
    object <- make_time_qmle_data(model)
    theta <- list(alpha = 0.7, gamma = 0.3)
    observed <- object@data@zoo.data[[1]]
    x <- as.numeric(observed)
    times <- as.numeric(index(observed))
    h <- deltat(observed)
    drift <- -theta$alpha * head(x, -1) + head(times, -1)
    variance <- h * exp(theta$gamma * head(times, -1))
    expected <- sum(dnorm(diff(x), mean = h * drift,
                          sd = sqrt(variance), log = TRUE))
    expect_equal(quasilogl(object, theta, rcpp = FALSE), expected)
    expect_equal(quasilogl(object, theta, rcpp = TRUE), expected)
  }
})

test_that("pointwise coefficients support scalar functions and custom time names", {
  model <- setModel(
    drift = "-alpha*x + scalar_clock(s)",
    diffusion = "exp(gamma * scalar_clock(s) / 2)",
    time.variable = "s", state.variable = "x", solve.variable = "x"
  )
  object <- make_time_qmle_data(model)
  observed <- object@data@zoo.data[[1]]
  times <- as.numeric(index(observed))
  env <- new.env()
  env$X <- matrix(as.numeric(observed), ncol = 1)
  env$scalar_clock <- function(s) {
    if (length(s) != 1L) stop("scalar input required")
    if (s > 2.2) s^2 else s
  }
  env$s <- -999
  theta <- list(alpha = 0.7, gamma = 0.3)
  clock_values <- vapply(times, env$scalar_clock, numeric(1))
  expect_equal(as.numeric(yuima:::drift.term(object, theta, env)),
               -theta$alpha * as.numeric(observed) + clock_values)
  expect_equal(as.numeric(yuima:::diffusion.term(object, theta, env)),
               exp(theta$gamma * clock_values / 2))
  expect_identical(env$s, -999)
})

test_that("qmle defaults to the C++ likelihood and both paths agree", {
  set.seed(123)
  x1.func <- function(t, x = 0) {
    if (length(t) != 1L || length(x) != 1L) stop("scalar input required")
    cos(2*pi*t)
  }
  x2.func <- function(t, x = 0) {
    if (length(t) != 1L || length(x) != 1L) stop("scalar input required")
    sin(2*pi*t)
  }
  model <- setModel(
    drift = "-alpha*x",
    diffusion = "exp((gamma1*x1.func(t,x)+gamma2*x2.func(t,x))/2)",
    time.variable = "t", state.variable = "x", solve.variable = "x"
  )
  simulation_model <- setModel(drift = "-alpha*x",
    diffusion = "exp((gamma1*cos(2*pi*t)+gamma2*sin(2*pi*t))/2)",
    time.variable = "t", state.variable = "x", solve.variable = "x")
  object <- suppressWarnings(simulate(
    simulation_model, sampling = setSampling(Terminal = 1, n = 400), xinit = 0,
    true.parameter = list(alpha = 3, gamma1 = -2, gamma2 = 3)
  ))
  object@model <- model
  args <- list(
    yuima = object, start = list(alpha = 1, gamma1 = 0, gamma2 = 0),
    lower = list(alpha = 0.0001, gamma1 = -10, gamma2 = -10),
    upper = list(alpha = 10, gamma1 = 10, gamma2 = 10), envir = environment()
  )
  default <- do.call(qmle, args)
  cpp <- do.call(qmle, c(args, list(rcpp = TRUE)))
  pointwise <- do.call(qmle, c(args, list(rcpp = FALSE)))
  expect_identical(formals(qmle)$rcpp, TRUE)
  expect_equal(coef(default), coef(cpp))
  expect_equal(coef(pointwise), coef(cpp), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(pointwise)), as.numeric(logLik(cpp)),
               tolerance = 1e-7)
})

test_that("C++ likelihood preserves state-by-noise diffusion ordering", {
  model <- setModel(
    drift = c("-alpha*x + s", "-alpha*y - s"),
    diffusion = matrix(c("1+s", "0.3", "0.2", "2+s", "0.4", "0.1"), 2, 3),
    state.variable = c("x", "y"), solve.variable = c("x", "y"),
    time.variable = "s"
  )
  values <- cbind(c(0.2, 0.1, 0.4, -0.1, 0.3, 0.5),
                  c(-0.1, 0.3, 0.2, 0.5, 0.1, -0.2))
  object <- suppressWarnings(setYuima(
    model = model, data = setData(ts(values, start = 2, deltat = 0.1)),
    variable_data_mapping = list(x = 1, y = 2)
  ))
  times <- as.numeric(index(object@data@zoo.data[[1]]))
  h <- deltat(object@data@zoo.data[[1]])
  expected <- sum(vapply(seq_len(nrow(values)-1), function(i) {
    a <- matrix(c(1+times[i], 0.3, 0.2, 2+times[i], 0.4, 0.1), 2, 3)
    b <- -0.7 * values[i, ] + c(times[i], -times[i])
    mvtnorm::dmvnorm(values[i+1, ] - values[i, ], mean = h*b,
                     sigma = h*tcrossprod(a), log = TRUE)
  }, numeric(1)))
  expect_equal(quasilogl(object, list(alpha = 0.7), rcpp = TRUE), expected)
  expect_equal(quasilogl(object, list(alpha = 0.7), rcpp = FALSE), expected)
})
