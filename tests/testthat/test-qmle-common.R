test_that("log Black-Scholes estimates shared diffusion parameters first", {
  increments <- c(0.05, -0.01, 0.04, 0.02, -0.03, 0.06, 0.01, 0.03)
  h <- 0.1
  sigma <- sqrt(mean(increments^2) / h)
  mu <- mean(increments) / h + sigma^2 / 2
  env <- new.env()
  env$sigma <- 99
  for (drift in c("mu - sigma^2/2", "-sigma^2/2")) {
    model <- setModel(drift = drift, diffusion = "sigma",
      state.variable = "x", solve.variable = "x")
    object <- suppressWarnings(setYuima(model = model,
      data = setData(ts(c(0, cumsum(increments)), deltat = h)),
      variable_data_mapping = list(x = 1)))
    start <- list(sigma = 0.3)
    lower <- list(sigma = 0.001)
    upper <- list(sigma = 2)
    expected <- c(sigma = sigma)
    if (grepl("mu", drift)) {
      start$mu <- 0.1
      lower$mu <- -2
      upper$mu <- 2
      expected <- c(expected, mu = mu)
    }
    for (rcpp in c(TRUE, FALSE)) {
      fit <- qmle(object, start = start, lower = lower, upper = upper,
                  envir = env, joint = FALSE, rcpp = rcpp)
      expect_equal(coef(fit)[names(expected)], expected, tolerance = 1e-4)
      expect_true(all(is.finite(diag(vcov(fit)))))
      expect_true(all(diag(vcov(fit)) > 0))
    }
  }
  expect_identical(env$sigma, 99)
})

test_that("log Black-Scholes still permits explicit joint estimation", {
  increments <- c(0.05, -0.01, 0.04, 0.02, -0.03, 0.06, 0.01, 0.03)
  h <- 0.1
  sigma <- sqrt(mean((increments - mean(increments))^2) / h)
  mu <- mean(increments) / h + sigma^2 / 2
  model <- setModel(drift = "mu - sigma^2/2", diffusion = "sigma",
    state.variable = "x", solve.variable = "x")
  object <- suppressWarnings(setYuima(model = model,
    data = setData(ts(c(0, cumsum(increments)), deltat = h)),
    variable_data_mapping = list(x = 1)))
  for (rcpp in c(TRUE, FALSE)) {
    fit <- qmle(object, start = list(mu = 0.1, sigma = 0.3),
      lower = list(mu = -2, sigma = 0.001),
      upper = list(mu = 2, sigma = 2), joint = TRUE, rcpp = rcpp)
    expect_equal(coef(fit)[c("sigma", "mu")], c(sigma = sigma, mu = mu),
                 tolerance = 1e-4)
  }
})
