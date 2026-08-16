library(yuima)

adaQmleArguments <- names(formals(adaQmle))
stopifnot("..." %in% adaQmleArguments,
          !"control" %in% adaQmleArguments)

initialChoice <- yuima:::.adaQmle_normalize_initial(c("W", "V"))
stopifnot(identical(initialChoice[["drift"]], "least_squares"),
          identical(initialChoice[["diffusion"]], "gaussian"),
          identical(yuima:::.adaQmle_normalize_schedule(c("b", "m"), 4L),
                    c("bayes", "optim", "optim", "optim")))

# Initial Bayes temperatures use q=max(p,2/G), and temper only the
# quasi-likelihood contribution rather than the prior.
temperature <- yuima:::.adaQmle_bayes_temperature(
  n0 = 20000, h = 1 / 390, p = 4, rate = 0.86
)
stopifnot(temperature$q == 4,
          isTRUE(all.equal(temperature$exponent, 2 / (4 * 0.86) - 1)),
          isTRUE(all.equal(temperature$diffusion,
                           20000^(2 / (4 * 0.86) - 1))),
          isTRUE(all.equal(temperature$drift,
                           (20000 / 390)^(2 / (4 * 0.86) - 1))),
          yuima:::.adaQmle_tempered_log_target(
            contrast = 10, logPrior = 3, temperature = 0.2
          ) == 2)
temperatureBoundary <- yuima:::.adaQmle_bayes_temperature(
  n0 = 100, h = 0.01, p = 2, rate = 0.5
)
stopifnot(temperatureBoundary$q == 4,
          temperatureBoundary$exponent == 0,
          temperatureBoundary$diffusion == 1,
          temperatureBoundary$drift == 1)

# MpCN uses the Gamma scale mixture and the radial Hastings correction used
# by adaBayes.  The deterministic comparison checks the exact proposal law.
mpcnCurrent <- c(a = 1.1, b = -0.4)
mpcnCenter <- c(a = 0.2, b = 0.3)
mpcnPreconditioner <- c(a = 1.5, b = 0.7)
mpcnRho <- 0.8
set.seed(1414)
currentRadius <- max(sum(((mpcnCurrent - mpcnCenter) * mpcnPreconditioner)^2),
                     1e-7)
mixingPrecision <- rgamma(1L, shape = length(mpcnCurrent) / 2,
                          scale = 2 / currentRadius)
expectedProposal <- mpcnCenter + sqrt(mpcnRho) *
  (mpcnCurrent - mpcnCenter) + rnorm(length(mpcnCurrent)) *
  sqrt((1 - mpcnRho) / mixingPrecision)
expectedRadius <- max(sum(((expectedProposal - mpcnCenter) *
                           mpcnPreconditioner)^2), 1e-7)
expectedLogHastings <- length(mpcnCurrent) / 2 *
  (log(expectedRadius) - log(currentRadius))
set.seed(1414)
mpcnProposal <- yuima:::.adaQmle_mpcn_proposal(
  mpcnCurrent, mpcnCenter, mpcnPreconditioner, mpcnRho
)
stopifnot(isTRUE(all.equal(mpcnProposal$par, expectedProposal)),
          isTRUE(all.equal(mpcnProposal$radius, expectedRadius)),
          isTRUE(all.equal(mpcnProposal$log.hastings, expectedLogHastings)))

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

# Model expressions are evaluated by vector chunks, including constant terms.
evaluationEnvironment <- new.env(parent = globalenv())
evaluationEnvironment$a <- 2
evaluationData <- matrix(c(-2, -1, 0, 1, 2), ncol = 1)
vectorizedTerms <- yuima:::adaEvalTermsCpp(
  expression(a * x^2, 3), "x", evaluationData,
  character(0), numeric(), evaluationEnvironment
)
stopifnot(isTRUE(all.equal(
  vectorizedTerms,
  cbind(2 * evaluationData[, 1]^2, rep(3, nrow(evaluationData)))
)))

# A scalar-only user function falls back to observation-by-observation eval.
evaluationEnvironment$scalarOnly <- function(x) {
  if (length(x) != 1L) stop("scalar input required")
  x + 1
}
fallbackTerms <- yuima:::adaEvalTermsCpp(
  expression(scalarOnly(x)), "x", evaluationData,
  character(0), numeric(), evaluationEnvironment
)
stopifnot(isTRUE(all.equal(fallbackTerms[, 1], evaluationData[, 1] + 1)))

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

set.seed(1619)
fitMpCN <- adaQmle(ouYuima, start = list(a = 0.9, b = 0.7), p = 2,
                   lower = list(a = 0.05, b = 0.05),
                   upper = list(a = 3, b = 3), estimator = "b",
                   prior = prior, mcmc = 20, rate = 0.5,
                   algorithm = "mpcn", center = list(a = 1.2, b = 0.7),
                   sd = list(a = 1, b = 1), path = TRUE,
                   control = list(maxit = 50))
mpcnStage <- fitMpCN@details$stages[[1L]]
mpcnDriftStage <- fitMpCN@details$stages[[2L]]
stopifnot(identical(fitMpCN@details$estimator, c("bayes", "bayes")),
          fitMpCN@details$q == 4,
          fitMpCN@details$rate == 0.5,
          fitMpCN@details$bayes.n0 ==
            floor((nrow(yuima:::onezoo(ouYuima)))^0.5),
          mpcnStage$q == 4,
          mpcnStage$temperature.exponent == 0,
          mpcnStage$temperature == 1,
          mpcnStage$n0 == floor((nrow(yuima:::onezoo(ouYuima)))^0.5),
          mpcnStage$accept.rate >= 0, mpcnStage$accept.rate <= 1,
          is.matrix(mpcnStage$path),
          mpcnDriftStage$q == 4,
          mpcnDriftStage$temperature == 1,
          mpcnDriftStage$accept.rate >= 0,
          mpcnDriftStage$accept.rate <= 1,
          is.matrix(mpcnDriftStage$path),
          all(is.finite(fitMpCN@coef)))

fitFixed <- adaQmle(ouYuima, start = list(a = 0.9), p = 2,
                    fixed = list(b = 0.8),
                    lower = list(a = 0.05, b = 0.05),
                    upper = list(a = 3, b = 3), control = list(maxit = 50))
stopifnot(identical(names(fitFixed@fixed), "b"), fitFixed@fixed[["b"]] == 0.8,
          identical(names(fitFixed@coef), "a"),
          isTRUE(fitFixed@details$stages[[2L]]$skipped))
