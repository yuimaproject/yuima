test_that("adaBayes supports scalar time-dependent coefficients at reduced rates", {
  scalar_wave <- function(s, x) {
    if (length(s) != 1L || length(x) != 1L) stop("scalar input required")
    cos(2*pi*s)
  }
  model <- setModel(drift = "-alpha*x + s",
    diffusion = "exp(gamma*scalar_wave(s,x)/2)",
    time.variable = "s", state.variable = "x", solve.variable = "x")
  set.seed(42)
  values <- cumsum(c(0, rnorm(80, sd = 0.1)))
  object <- suppressWarnings(setYuima(model = model,
    data = setData(ts(values, start = 2, deltat = 0.025)),
    variable_data_mapping = list(x = 1)))
  for (rate in c(1, 0.7)) {
    args <- list(yuima = object, start = list(alpha = 1, gamma = 0.1),
      lower = list(alpha = 0.01, gamma = -3),
      upper = list(alpha = 10, gamma = 3),
      envir = environment(), rate = rate, mcmc = 20,
      sd = c(0.02, 0.02), path = TRUE)
    set.seed(123)
    cpp <- do.call(adaBayes, c(args, list(rcpp = TRUE)))
    set.seed(123)
    r <- do.call(adaBayes, c(args, list(rcpp = FALSE)))
    expect_true(all(is.finite(coef(cpp))))
    expect_equal(coef(cpp), coef(r), tolerance = 1e-5)
    expect_equal(cpp@mcmc, r@mcmc, tolerance = 1e-5)
  }
})

test_that("reduced likelihood evaluates only the matching observation times", {
  model <- setModel(drift = "-alpha*x + s", diffusion = "exp(gamma*s/2)",
    time.variable = "s", state.variable = "x", solve.variable = "x")
  values <- c(0.2, 0.1, 0.4, -0.1, 0.3, 0.5)
  object <- suppressWarnings(setYuima(model = model,
    data = setData(ts(values, start = 2, deltat = 0.1)),
    variable_data_mapping = list(x = 1)))
  env <- new.env()
  env$X <- matrix(values[1:4], ncol = 1)
  env$deltaX <- matrix(diff(values[1:4]), ncol = 1)
  env$time <- as.numeric(index(object@data@zoo.data[[1]]))[1:4]
  env$h <- deltat(object@data@zoo.data[[1]])
  env$Cn.r <- rep(1, 3)
  theta <- list(alpha = 0.7, gamma = 0.3)
  b <- -theta$alpha * values[1:3] + env$time[1:3]
  variance <- env$h * exp(theta$gamma * env$time[1:3])
  expected <- -sum(dnorm(diff(values[1:4]), env$h*b, sqrt(variance), log = TRUE))
  for (rcpp in c(TRUE, FALSE)) {
    expect_equal(yuima:::minusquasilogl(object, theta, env = env, rcpp = rcpp),
                 expected)
  }
})
