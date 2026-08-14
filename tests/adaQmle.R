library(yuima)

initialChoice <- yuima:::.adaQmle_normalize_initial(c("W", "V"))
stopifnot(identical(initialChoice[["drift"]], "least_squares"),
          identical(initialChoice[["diffusion"]], "gaussian"),
          identical(yuima:::.adaQmle_normalize_schedule(c("b", "m"), 4L),
                    c("bayes", "optim", "optim", "optim")))

# The C++ Gaussian and least-squares kernels agree with their direct formulas.
dx <- matrix(c(0.1, -0.2, 0.3), ncol = 1)
drift <- matrix(0, nrow(dx), 1)
diffusion <- matrix(2, nrow(dx), 1)
h <- 0.1
v <- yuima:::adaContrastCpp(dx, drift, diffusion, list(), list(), h,
                            "initial_gaussian_diffusion", 1L)
stopifnot(isTRUE(all.equal(v, sum(log(4) + dx[, 1]^2 / (4 * h)),
                           tolerance = 1e-12)))
w <- yuima:::adaContrastCpp(dx, drift, diffusion, list(), list(), h,
                            "initial_least_squares_diffusion", 1L)
stopifnot(isTRUE(all.equal(w, sum((dx[, 1]^2 - 4 * h)^2 / h^2),
                           tolerance = 1e-12)))

# p determines k0 and l0; there is no public order argument.
model <- setModel(drift = "-b*x", diffusion = "a",
                  state.variable = "x", solve.variable = "x")
symbolic <- yuima:::.adaQmle_symbolic_terms(model, 4L)
stopifnot(symbolic$k0 == 2L, symbolic$l0 == 1L)

# U4 agrees with the closed formal expansion for a one-dimensional OU model.
x <- matrix(c(0.2, -0.1, 0.4), ncol = 1)
thetaEnvironment <- new.env(parent = globalenv())
thetaEnvironment$a <- 1.3
thetaEnvironment$b <- 0.7
bValues <- yuima:::adaEvalTermsCpp(symbolic$drift, "x", x, "t", 0:2,
                                   thetaEnvironment)
sigmaValues <- yuima:::adaEvalTermsCpp(symbolic$diffusion, "x", x, "t", 0:2,
                                       thetaEnvironment)
meanTerms <- lapply(symbolic$mean, yuima:::adaEvalTermsCpp, state = "x", data = x,
                    timeVariable = "t", time = 0:2, env = thetaEnvironment)
momentTerms <- lapply(symbolic$moment, yuima:::adaEvalTermsCpp, state = "x", data = x,
                      timeVariable = "t", time = 0:2, env = thetaEnvironment)
u <- yuima:::adaContrastCpp(dx, bValues, sigmaValues, meanTerms, momentTerms, h,
                            "coupled", 4L)
beta <- thetaEnvironment$b
sigma <- thetaEnvironment$a
uDirect <- sum(log(sigma^2) - beta * h + beta^2 * h^2 / 6 +
  (dx[, 1] - h * (-beta * x[, 1]) - h^2 * beta^2 * x[, 1] / 2)^2 / h *
  (1 + beta * h + beta^2 * h^2 / 3) / sigma^2)
stopifnot(isTRUE(all.equal(u, uDirect, tolerance = 1e-11)))

# A common parameter is assigned only to the diffusion block and is plugged
# into the drift block, preventing the shared-parameter bug in qmle.
commonModel <- setModel(drift = "-c*x", diffusion = "sqrt(c)",
                        state.variable = "x", solve.variable = "x")
sampling <- setSampling(Terminal = 1, n = 100)
commonYuima <- setYuima(model = commonModel, sampling = sampling)
set.seed(2718)
commonYuima <- simulate(commonYuima, xinit = 1,
                        true.parameter = list(c = 1), sampling = sampling)
fit <- adaQmle(commonYuima, start = list(c = 0.8), p = 2,
               lower = list(c = 0.05), upper = list(c = 3),
               control = list(maxit = 50))
stopifnot(is(fit, "yuima.qmle"), identical(fit@details$p, 2L),
          identical(unname(fit@details$order), c(0, 1)),
          identical(fit@details$common.parameters, "c"),
          isTRUE(fit@details$stages[[2L]]$skipped),
          is.finite(fit@fullcoef[["c"]]))

ouModel <- setModel(drift = "-b*x", diffusion = "a",
                    state.variable = "x", solve.variable = "x")
ouYuima <- setYuima(model = ouModel, sampling = sampling)
set.seed(31415)
ouYuima <- simulate(ouYuima, xinit = 0.5,
                    true.parameter = list(a = 1, b = 0.8), sampling = sampling)
fitPlugin <- adaQmle(ouYuima, start = list(a = 0.9, b = 0.7), p = 3,
                     lower = list(a = 0.05, b = 0.05),
                     upper = list(a = 3, b = 3),
                     refinement = "plugin", expansion = "progressive",
                     estimator = "m", control = list(maxit = 50))
stopifnot(identical(fitPlugin@details$contrast.order, 1:3),
          identical(unname(fitPlugin@details$order), c(1, 1)),
          all(is.finite(fitPlugin@coef)),
          is(summary(fitPlugin), "summary.yuima.qmle"))

fitCoupled <- adaQmle(ouYuima, start = list(a = 0.9, b = 0.7), p = 4,
                      lower = list(a = 0.05, b = 0.05),
                      upper = list(a = 3, b = 3),
                      refinement = "U", expansion = "full",
                      estimator = "m", control = list(maxit = 50))
stopifnot(identical(fitCoupled@details$contrast.order, c(1L, 2L, 4L, 4L)),
          identical(unname(fitCoupled@details$order), c(1, 2)),
          all(is.finite(fitCoupled@coef)))

prior <- list(a = list(measure.type = "code", df = "dunif(z,0.05,3)"),
              b = list(measure.type = "code", df = "dunif(z,0.05,3)"))
set.seed(1618)
fitBayes <- adaQmle(ouYuima, start = list(a = 0.9, b = 0.7), p = 2,
                    lower = list(a = 0.05, b = 0.05),
                    upper = list(a = 3, b = 3), estimator = c("b", "m"),
                    prior = prior, mcmc = 20, path = TRUE,
                    sd = list(a = 0.03, b = 0.03),
                    control = list(maxit = 50))
stopifnot(identical(fitBayes@details$estimator, c("bayes", "optim")),
          is.matrix(fitBayes@details$stages[[1L]]$path),
          all(is.finite(fitBayes@coef)))

fitFixed <- adaQmle(ouYuima, start = list(a = 0.9), p = 2,
                    fixed = list(b = 0.8),
                    lower = list(a = 0.05, b = 0.05),
                    upper = list(a = 3, b = 3), control = list(maxit = 50))
stopifnot(identical(names(fitFixed@fixed), "b"), fitFixed@fixed[["b"]] == 0.8,
          identical(names(fitFixed@coef), "a"),
          isTRUE(fitFixed@details$stages[[2L]]$skipped))
