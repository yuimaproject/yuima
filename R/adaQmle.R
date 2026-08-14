# Adaptive quasi-likelihood estimation for diffusion processes

.adaQmle_scalar_expression <- function(x) {
  if (is.expression(x)) x[[1L]] else x
}

.adaQmle_is_zero <- function(x) {
  is.numeric(x) && length(x) == 1L && isTRUE(x == 0)
}

.adaQmle_add <- function(x, y) {
  if (.adaQmle_is_zero(x)) return(y)
  if (.adaQmle_is_zero(y)) return(x)
  call("+", x, y)
}

.adaQmle_multiply <- function(x, y) {
  if (.adaQmle_is_zero(x) || .adaQmle_is_zero(y)) return(0)
  if (is.numeric(x) && length(x) == 1L && x == 1) return(y)
  if (is.numeric(y) && length(y) == 1L && y == 1) return(x)
  call("*", x, y)
}

.adaQmle_derivative <- function(x, variable) {
  calculus::derivative(as.expression(x), var = variable, deparse = FALSE)[[1L]]
}

.adaQmle_generator <- function(x, drift, covariance, state) {
  dimension <- length(state)
  value <- 0
  for (i in seq_len(dimension)) {
    value <- .adaQmle_add(
      value,
      .adaQmle_multiply(drift[[i]], .adaQmle_derivative(x, state[[i]]))
    )
  }
  second <- 0
  for (column in seq_len(dimension)) {
    firstDerivative <- .adaQmle_derivative(x, state[[column]])
    for (row in seq_len(dimension)) {
      second <- .adaQmle_add(
        second,
        .adaQmle_multiply(
          covariance[[row, column]],
          .adaQmle_derivative(firstDerivative, state[[row]])
        )
      )
    }
  }
  .adaQmle_add(value, .adaQmle_multiply(0.5, second))
}

.adaQmle_symbolic_terms <- function(model, p) {
  state <- model@state.variable
  dimension <- length(state)
  drift <- lapply(as.list(model@drift), .adaQmle_scalar_expression)
  diffusion <- lapply(model@diffusion, function(row) {
    lapply(as.list(row), .adaQmle_scalar_expression)
  })
  wienerDimension <- length(diffusion[[1L]])

  covariance <- matrix(vector("list", dimension * dimension), dimension, dimension)
  for (column in seq_len(dimension)) {
    for (row in seq_len(dimension)) {
      entry <- 0
      for (noise in seq_len(wienerDimension)) {
        entry <- .adaQmle_add(
          entry,
          .adaQmle_multiply(diffusion[[row]][[noise]], diffusion[[column]][[noise]])
        )
      }
      covariance[[row, column]] <- entry
    }
  }

  # k0 and l0 are the orders in Uchida--Yoshida (2012).  Computing
  # moments through k0 + 1 supplies both the U and V constructions.
  k0 <- floor(p / 2)
  l0 <- floor((p - 1) / 2)
  highest <- k0 + 1L
  meanTerms <- vector("list", max(0L, highest - 1L))
  momentTerms <- vector("list", max(0L, highest - 1L))

  meanCurrent <- lapply(state, as.name)
  for (j in seq_len(highest)) {
    meanCurrent <- lapply(
      meanCurrent,
      .adaQmle_generator,
      drift = drift,
      covariance = covariance,
      state = state
    )
    if (j >= 2L) {
      meanTerms[[j - 1L]] <- as.expression(lapply(meanCurrent, function(x) call("/", x, factorial(j))))
    }
  }

  anchor <- paste0(".adaQmle_anchor_", seq_len(dimension))
  replacement <- setNames(lapply(state, as.name), anchor)
  momentCurrent <- matrix(vector("list", dimension * dimension), dimension, dimension)
  for (column in seq_len(dimension)) {
    for (row in seq_len(dimension)) {
      momentCurrent[[row, column]] <- .adaQmle_multiply(
        call("-", as.name(state[[row]]), as.name(anchor[[row]])),
        call("-", as.name(state[[column]]), as.name(anchor[[column]]))
      )
    }
  }
  for (j in seq_len(highest)) {
    for (column in seq_len(dimension)) {
      for (row in seq_len(dimension)) {
        momentCurrent[[row, column]] <- .adaQmle_generator(
          momentCurrent[[row, column]], drift, covariance, state
        )
      }
    }
    if (j >= 2L) {
      flattened <- vector("list", dimension * dimension)
      for (column in seq_len(dimension)) {
        for (row in seq_len(dimension)) {
          term <- do.call(substitute, list(momentCurrent[[row, column]], replacement))
          flattened[[row + dimension * (column - 1L)]] <- call("/", term, factorial(j))
        }
      }
      momentTerms[[j - 1L]] <- as.expression(flattened)
    }
  }

  diffusionTerms <- vector("list", dimension * wienerDimension)
  for (noise in seq_len(wienerDimension)) {
    for (row in seq_len(dimension)) {
      diffusionTerms[[row + dimension * (noise - 1L)]] <- diffusion[[row]][[noise]]
    }
  }

  list(
    drift = as.expression(drift),
    diffusion = as.expression(diffusionTerms),
    mean = meanTerms,
    moment = momentTerms,
    k0 = k0,
    l0 = l0
  )
}

.adaQmle_normalize_initial <- function(x) {
  if (is.list(x)) x <- unlist(x, use.names = TRUE)
  if (!length(x)) x <- "gaussian"
  if (length(x) == 1L) x <- rep(x, 2L)
  if (length(x) != 2L) yuima.stop("'initial.contrast' must have length one or two.")
  if (!is.null(names(x)) && all(c("drift", "diffusion") %in% names(x))) {
    x <- x[c("drift", "diffusion")]
  }
  aliases <- c(v = "gaussian", gaussian = "gaussian",
               w = "least_squares", least_squares = "least_squares",
               `least-squares` = "least_squares", lse = "least_squares")
  normalized <- unname(aliases[tolower(x)])
  if (anyNA(normalized)) {
    yuima.stop("'initial.contrast' must be 'gaussian' (V) or 'least_squares' (W).")
  }
  # The public two-entry convention is c(drift, diffusion), matching the
  # mathematical coefficient pair. Adaptive execution still starts with the
  # diffusion block.
  setNames(normalized, c("drift", "diffusion"))
}

.adaQmle_normalize_schedule <- function(estimator, stages) {
  if (is.list(estimator)) estimator <- unlist(estimator, use.names = FALSE)
  if (!length(estimator)) estimator <- "optim"
  aliases <- c(m = "optim", mle = "optim", optim = "optim",
               b = "bayes", bayes = "bayes")
  estimator <- unname(aliases[tolower(estimator)])
  if (anyNA(estimator)) yuima.stop("'estimator' entries must be 'm'/'optim' or 'b'/'bayes'.")
  if (length(estimator) > stages) yuima.stop("'estimator' has more entries than adaptive steps.")
  c(estimator, rep(tail(estimator, 1L), stages - length(estimator)))
}

.adaQmle_named_numeric <- function(x, allowed, argument, default) {
  result <- setNames(rep(default, length(allowed)), allowed)
  if (missing(x) || is.null(x) || !length(x)) return(result)
  if (is.list(x)) x <- unlist(x)
  if (is.null(names(x)) || any(names(x) == "")) {
    yuima.stop(paste0("'", argument, "' must be a named list or named vector."))
  }
  if (any(!names(x) %in% allowed)) {
    yuima.stop(paste0("Some names in '", argument, "' are not model parameters."))
  }
  result[names(x)] <- as.numeric(x)
  result
}

.adaQmle_log_prior <- function(prior, parameterNames, envir) {
  if (missing(prior) || is.null(prior)) return(function(theta) 0)
  if (!is.list(prior) || any(!parameterNames %in% names(prior))) {
    yuima.stop("'prior' must contain a named density specification for every estimated parameter.")
  }
  expressions <- lapply(parameterNames, function(parameter) {
    specification <- prior[[parameter]]
    if (!is.list(specification) || !identical(specification$measure.type, "code") ||
        !is.character(specification$df) || length(specification$df) != 1L) {
      yuima.stop("Each prior must have measure.type='code' and a character 'df'.")
    }
    parse(text = specification$df)[[1L]]
  })
  names(expressions) <- parameterNames
  function(theta) {
    value <- 0
    for (parameter in parameterNames) {
      evaluationEnvironment <- new.env(parent = envir)
      evaluationEnvironment$z <- unname(theta[[parameter]])
      assign(parameter, unname(theta[[parameter]]), envir = evaluationEnvironment)
      density <- eval(expressions[[parameter]], envir = evaluationEnvironment)
      if (length(density) != 1L || !is.finite(density) || density <= 0) return(-Inf)
      value <- value + log(density)
    }
    value
  }
}

.adaQmle_mcmc <- function(start, objective, logPrior, lower, upper, mcmc,
                          algorithm, center, proposalSd, rho, temperature, path) {
  dimension <- length(start)
  chain <- matrix(NA_real_, mcmc, dimension, dimnames = list(NULL, names(start)))
  chain[1L, ] <- start
  current <- start
  logTarget <- function(x) {
    names(x) <- names(start)
    priorValue <- logPrior(x)
    if (!is.finite(priorValue)) return(-Inf)
    contrast <- objective(x)
    if (!is.finite(contrast)) return(-Inf)
    -0.5 * temperature * contrast + priorValue
  }
  currentTarget <- logTarget(current)
  if (!is.finite(currentTarget)) yuima.stop("The initial value has zero posterior density.")
  accepted <- 0L
  covariance <- diag(proposalSd^2, dimension)
  sqrtRho <- sqrt(rho)
  proposalCovariance <- (1 - rho) * covariance

  for (iteration in 2:mcmc) {
    if (algorithm == "randomwalk") {
      proposal <- current + as.numeric(mvtnorm::rmvnorm(1L, sigma = covariance))
      logHastings <- 0
    } else {
      proposalMean <- center + sqrtRho * (current - center)
      proposal <- as.numeric(mvtnorm::rmvnorm(1L, mean = proposalMean, sigma = proposalCovariance))
      reverseMean <- center + sqrtRho * (proposal - center)
      logHastings <- mvtnorm::dmvnorm(current, mean = reverseMean,
                                      sigma = proposalCovariance, log = TRUE) -
        mvtnorm::dmvnorm(proposal, mean = proposalMean,
                         sigma = proposalCovariance, log = TRUE)
    }
    names(proposal) <- names(start)
    if (all(proposal >= lower) && all(proposal <= upper)) {
      proposalTarget <- logTarget(proposal)
      if (is.finite(proposalTarget) &&
          log(stats::runif(1L)) < proposalTarget - currentTarget + logHastings) {
        current <- proposal
        currentTarget <- proposalTarget
        accepted <- accepted + 1L
      }
    }
    chain[iteration, ] <- current
  }
  retained <- seq.int(floor(mcmc / 2) + 1L, mcmc)
  estimate <- colMeans(chain[retained, , drop = FALSE])
  covarianceEstimate <- if (length(retained) > 1L) {
    stats::cov(chain[retained, , drop = FALSE])
  } else {
    matrix(NA_real_, dimension, dimension)
  }
  list(
    par = estimate,
    vcov = covarianceEstimate,
    accept.rate = accepted / (mcmc - 1L),
    path = if (path) chain else NULL
  )
}

#' Adaptive quasi-maximum likelihood estimation
#'
#' Estimates diffusion parameters first and then plugs them into each drift
#' step.  Higher-order U or V contrasts are selected independently from the
#' initial Gaussian (V) or least-squares (W) contrasts.
#'
#' @param yuima A `yuima` object for a diffusion model.
#' @param start Named starting values.
#' @param p Integer exponent in the sampling condition `n h^p -> 0`.
#' @param lower,upper Named parameter bounds.
#' @param fixed Named values held fixed.
#' @param initial.contrast Initial contrast for diffusion and drift.  Use
#'   `"gaussian"`/`"V"` or `"least_squares"`/`"W"`.
#' @param refinement `"plugin"`/`"V"` or `"coupled"`/`"U"`.
#' @param expansion `"progressive"` or `"terminal"`/`"full"`.
#' @param estimator Per-step `"m"`/`"optim"` or `"b"`/`"bayes"` schedule.
#' @param prior Prior specifications as accepted by `adaBayes`.
#' @param method Optimization method passed to [stats::optim()].
#' @param control Control list passed to [stats::optim()].
#' @param envir Parent environment for model evaluation.
#' @param mcmc Number of MCMC draws for a Bayes step.
#' @param iteration Optional alias for `mcmc`, retained for `adaBayes`
#'   compatibility.
#' @param rate Initial Bayes sample exponent, as in `adaBayes`.
#' @param algorithm `"randomwalk"` or `"mpcn"`.
#' @param center Optional named centers for MpCN.
#' @param sd Optional named proposal standard deviations.
#' @param rho MpCN persistence parameter.
#' @param path Whether to retain MCMC paths in `details`.
#' @details
#' Write `Delta_i X = X_{t_i} - X_{t_{i-1}}`, let `h` be the sampling
#' interval, and put `a(x, alpha) = sigma(x, alpha) sigma(x, alpha)'`.
#' Up to parameter-independent constants, the initial Gaussian contrasts are
#' \deqn{V_{alpha,0} = sum_i [log det(a_i) +
#' h^{-1} (Delta_i X)' a_i^{-1} (Delta_i X)]}
#' and
#' \deqn{V_{beta,0} = sum_i h^{-1} (Delta_i X-h b_i)'
#' a_i^{-1} (Delta_i X-h b_i),}
#' where the diffusion estimate is plugged into `a_i` in the drift contrast.
#' The corresponding least-squares contrasts are
#' \deqn{W_{alpha,0} = sum_i ||h^{-1} Delta_i X (Delta_i X)' - a_i||_F^2}
#' and
#' \deqn{W_{beta,0} = sum_i h^{-1} ||Delta_i X-h b_i||^2.}
#'
#' For higher-order steps, let `mu_i^{(r)}` and `C_i^{(r)}` denote the
#' order-`r` expansions of the conditional mean of `Delta_i X` and its
#' conditional covariance, and write `Cbar_i^{(r)} = C_i^{(r)}/h`.
#' The coupled U contrast has the Gaussian form
#' \deqn{U^{(r)}(theta) = sum_i [log det(Cbar_i^{(r)}(theta)) +
#' h^{-1} (Delta_i X-mu_i^{(r)}(theta))'
#' Cbar_i^{(r)}(theta)^{-1} (Delta_i X-mu_i^{(r)}(theta))],}
#' interpreted as the required truncated power series in `h`.  With
#' `refinement = "coupled"`, all correction terms are recomputed at the
#' current candidate.  With `refinement = "plugin"`, the leading coefficient
#' for the parameter block currently being optimized remains variable, while
#' higher-order mean or second-moment corrections are fixed at the preceding
#' adaptive estimates.  The implemented contrast is twice the negative
#' quasi-log likelihood, up to parameter-independent constants, and is
#' therefore minimized.
#' @references
#' Uchida, M. and Yoshida, N. (2012). Adaptive estimation of an ergodic
#' diffusion process based on sampled data. *Stochastic Processes and their
#' Applications*, 122, 2885--2924.
#'
#' Kamatani, K. and Uchida, M. (2015). Hybrid multi-step estimators for
#' stochastic differential equations based on sampled data. *Statistical
#' Inference for Stochastic Processes*, 18, 177--204.
#'
#' Kaino, Y., Uchida, M. and Yoshida, Y. (2017). Hybrid estimation for an
#' ergodic diffusion process based on reduced data. *Bulletin of Informatics
#' and Cybernetics*, 49, 89--118.
#' @return A `yuima.qmle` object.  Its `details` slot records `p`, the derived
#'   orders, all adaptive stages and any retained MCMC paths.
#' @export
adaQmle <- function(yuima, start, p, lower, upper, fixed = list(),
                    initial.contrast = c(drift = "gaussian", diffusion = "gaussian"),
                    refinement = "plugin",
                    expansion = c("progressive", "terminal"),
                    estimator = "optim", prior = NULL,
                    method = "L-BFGS-B", control = list(),
                    envir = globalenv(), mcmc = 1000L, iteration = NULL, rate = 1,
                    algorithm = c("randomwalk", "mpcn"), center = NULL,
                    sd = NULL, rho = 0.8, path = FALSE) {
  call <- match.call()
  if (missing(yuima) || !is(yuima, "yuima")) yuima.stop("'yuima' must be a yuima object.")
  if (missing(start) || !is.list(start) || is.null(names(start))) {
    yuima.stop("'start' must be a named list.")
  }
  if (missing(p) || length(p) != 1L || !is.finite(p) || p < 2 || p != as.integer(p)) {
    yuima.stop("'p' must be an integer greater than or equal to two.")
  }
  p <- as.integer(p)
  if (!length(refinement)) yuima.stop("'refinement' cannot be empty.")
  refinementAliases <- c(v = "plugin", plugin = "plugin",
                         u = "coupled", coupled = "coupled")
  refinement <- unname(refinementAliases[tolower(refinement[[1L]])])
  if (is.na(refinement)) yuima.stop("'refinement' must be 'plugin' (V) or 'coupled' (U).")
  if (!length(expansion)) yuima.stop("'expansion' cannot be empty.")
  expansionAliases <- c(progressive = "progressive", increasing = "progressive",
                        stepwise = "progressive", terminal = "terminal",
                        full = "terminal")
  expansion <- unname(expansionAliases[tolower(expansion[[1L]])])
  if (is.na(expansion)) yuima.stop("'expansion' must be 'progressive' or 'terminal'/'full'.")
  algorithm <- match.arg(tolower(algorithm), c("randomwalk", "mpcn"))
  initial.contrast <- .adaQmle_normalize_initial(initial.contrast)
  schedule <- .adaQmle_normalize_schedule(estimator, p)

  model <- yuima@model
  if (is(model, "yuima.carma") || is.COGARCH(yuima) || is.PPR(yuima) || is.Poisson(yuima) ||
      length(model@parameter@jump) || length(model@parameter@measure) || length(model@measure.type)) {
    yuima.stop("'adaQmle' currently supports diffusion models without jumps or point processes.")
  }
  diffusionParameters <- model@parameter@diffusion
  driftParameters <- model@parameter@drift
  if (!length(diffusionParameters) || !length(driftParameters)) {
    yuima.stop("Both drift and diffusion parameters are required.")
  }
  commonParameters <- unique(c(model@parameter@common,
                               intersect(diffusionParameters, driftParameters)))
  allParameters <- unique(c(diffusionParameters, driftParameters))

  if (!is.list(fixed) && is.atomic(fixed)) fixed <- as.list(fixed)
  if (!is.list(fixed) || (length(fixed) && is.null(names(fixed)))) {
    yuima.stop("'fixed' must be a named list.")
  }
  if (any(!names(fixed) %in% allParameters)) {
    yuima.stop("Some named arguments in 'fixed' are not model parameters.")
  }
  freeParameters <- setdiff(allParameters, names(fixed))
  if (any(!names(start) %in% allParameters) || any(!freeParameters %in% names(start))) {
    yuima.stop("'start' must provide every non-fixed drift and diffusion parameter.")
  }
  current <- setNames(numeric(length(allParameters)), allParameters)
  current[names(start)] <- as.numeric(unlist(start))
  if (length(fixed)) current[names(fixed)] <- as.numeric(unlist(fixed))
  if (any(!is.finite(current))) yuima.stop("All starting and fixed values must be finite.")

  if (missing(lower)) {
    lowerValues <- setNames(rep(-Inf, length(freeParameters)), freeParameters)
  } else {
    lowerValues <- .adaQmle_named_numeric(lower, allParameters, "lower", -Inf)[freeParameters]
  }
  if (missing(upper)) {
    upperValues <- setNames(rep(Inf, length(freeParameters)), freeParameters)
  } else {
    upperValues <- .adaQmle_named_numeric(upper, allParameters, "upper", Inf)[freeParameters]
  }
  if (any(lowerValues >= upperValues)) yuima.stop("Every lower bound must be smaller than its upper bound.")
  if (any(current[freeParameters] < lowerValues | current[freeParameters] > upperValues)) {
    yuima.stop("Starting values must lie inside the parameter bounds.")
  }

  data <- as.matrix(onezoo(yuima))
  if (!is.numeric(data) || nrow(data) < 3L || anyNA(data) || any(!is.finite(data))) {
    yuima.stop("The observed diffusion data must be a finite numeric matrix with at least three rows.")
  }
  state <- model@state.variable
  if (ncol(data) != length(state)) yuima.stop("The data dimension does not match the model state dimension.")
  observationTime <- as.numeric(index(yuima@data@zoo.data[[1L]]))
  if (length(observationTime) != nrow(data)) yuima.stop("Observation times do not match the data.")
  stepSizes <- diff(observationTime)
  h <- stepSizes[[1L]]
  if (h <= 0 || max(abs(stepSizes - h)) > sqrt(.Machine$double.eps) * max(1, abs(h))) {
    yuima.stop("'adaQmle' currently requires an equidistant observation grid.")
  }
  increments <- data[-1L, , drop = FALSE] - data[-nrow(data), , drop = FALSE]
  stateData <- data[-nrow(data), , drop = FALSE]
  stateTime <- observationTime[-length(observationTime)]
  numberIncrements <- nrow(increments)

  if (!is.numeric(rate) || length(rate) != 1L || rate <= 0 || rate > 1) {
    yuima.stop("'rate' must satisfy 0 < rate <= 1.")
  }
  if (any(schedule == "bayes")) {
    if (!is.null(iteration)) mcmc <- iteration
    mcmc <- as.integer(mcmc)
    if (length(mcmc) != 1L || is.na(mcmc) || mcmc < 4L) yuima.stop("'mcmc' must be at least four.")
    if (!is.numeric(rho) || length(rho) != 1L || rho <= 0 || rho >= 1) {
      yuima.stop("'rho' must satisfy 0 < rho < 1.")
    }
  }

  symbolic <- .adaQmle_symbolic_terms(model, p)
  timeVariable <- model@time.variable
  if (!length(timeVariable)) timeVariable <- character(0)

  evaluateTerms <- function(theta, observations = numberIncrements, corrections = TRUE) {
    evaluationEnvironment <- new.env(parent = envir)
    list2env(as.list(theta), envir = evaluationEnvironment)
    rows <- seq_len(observations)
    drift <- adaEvalTermsCpp(symbolic$drift, state, stateData[rows, , drop = FALSE],
                             timeVariable, stateTime[rows], evaluationEnvironment)
    diffusion <- adaEvalTermsCpp(symbolic$diffusion, state, stateData[rows, , drop = FALSE],
                                 timeVariable, stateTime[rows], evaluationEnvironment)
    if (!corrections) return(list(drift = drift, diffusion = diffusion,
                                  mean = list(), moment = list()))
    meanTerms <- lapply(symbolic$mean, adaEvalTermsCpp, state = state,
                        data = stateData[rows, , drop = FALSE], timeVariable = timeVariable,
                        time = stateTime[rows], env = evaluationEnvironment)
    momentTerms <- lapply(symbolic$moment, adaEvalTermsCpp, state = state,
                          data = stateData[rows, , drop = FALSE], timeVariable = timeVariable,
                          time = stateTime[rows], env = evaluationEnvironment)
    list(drift = drift, diffusion = diffusion, mean = meanTerms, moment = momentTerms)
  }

  diffusionFree <- setdiff(diffusionParameters, names(fixed))
  # A common parameter is estimated with diffusion and thereafter plugged in.
  driftFree <- setdiff(setdiff(driftParameters, commonParameters), names(fixed))
  stageTargets <- ifelse(seq_len(p) %% 2L == 1L, "diffusion", "drift")
  stages <- vector("list", p)
  covariance <- matrix(0, length(freeParameters), length(freeParameters),
                       dimnames = list(freeParameters, freeParameters))
  bayesInitialObservations <- max(1L, floor((numberIncrements + 1L)^rate) - 1L)
  bayesInitialObservations <- min(numberIncrements, bayesInitialObservations)

  logPriorAll <- .adaQmle_log_prior(prior, freeParameters, envir)
  stageOrder <- integer(p)
  for (stage in seq_len(p)) {
    target <- stageTargets[[stage]]
    parameters <- if (target == "diffusion") diffusionFree else driftFree
    stageOrder[[stage]] <- if (stage <= 2L) {
      stage
    } else if (expansion == "progressive") {
      stage
    } else if (refinement == "coupled") {
      p
    } else if (target == "diffusion") {
      2L * symbolic$l0 + 1L
    } else {
      2L * symbolic$k0
    }
    if (!length(parameters)) {
      stages[[stage]] <- list(stage = stage, target = target, order = stageOrder[[stage]],
                              estimator = schedule[[stage]], skipped = TRUE,
                              reason = "all parameters in this block are fixed or plugged in")
      next
    }

    useBayesPilot <- schedule[[stage]] == "bayes" && stage <= 2L
    observations <- if (useBayesPilot) bayesInitialObservations else numberIncrements
    rows <- seq_len(observations)
    snapshot <- current
    fixedCorrection <- NULL
    if (stage > 2L && refinement == "plugin") {
      fixedCorrection <- evaluateTerms(snapshot, observations, corrections = TRUE)
    }

    contrastName <- if (stage == 1L) {
      paste0("initial_", initial.contrast[["diffusion"]], "_diffusion")
    } else if (stage == 2L) {
      paste0("initial_", initial.contrast[["drift"]], "_drift")
    } else if (refinement == "plugin") {
      paste0("plugin_", target)
    } else {
      "coupled"
    }

    objective <- function(candidate) {
      theta <- snapshot
      theta[parameters] <- as.numeric(candidate)
      core <- evaluateTerms(theta, observations, corrections = refinement == "coupled" && stage > 2L)
      corrections <- if (!is.null(fixedCorrection)) fixedCorrection else core
      adaContrastCpp(increments[rows, , drop = FALSE], core$drift, core$diffusion,
                     corrections$mean, corrections$moment, h, contrastName,
                     stageOrder[[stage]])
    }
    initial <- current[parameters]
    stageLower <- lowerValues[parameters]
    stageUpper <- upperValues[parameters]

    if (schedule[[stage]] == "optim") {
      fit <- stats::optim(initial, objective, method = method,
                          lower = stageLower, upper = stageUpper,
                          hessian = TRUE, control = control)
      estimate <- setNames(fit$par, parameters)
      # The C++ kernels return twice the negative quasi-log likelihood, as in
      # the pilot implementation. Convert its Hessian to the nll scale.
      stageCovariance <- tryCatch(2 * solve(fit$hessian), error = function(e) {
        matrix(NA_real_, length(parameters), length(parameters))
      })
      stages[[stage]] <- c(list(stage = stage, target = target,
                                order = stageOrder[[stage]], contrast = contrastName,
                                estimator = "optim", observations = observations), fit)
    } else {
      stageLogPrior <- function(candidate) {
        theta <- current[freeParameters]
        theta[parameters] <- candidate
        logPriorAll(theta)
      }
      proposalSd <- if (is.null(sd)) {
        pmax(abs(initial) * 0.05, 0.02)
      } else {
        suppliedSd <- unlist(sd)
        if (is.null(names(suppliedSd)) || any(!parameters %in% names(suppliedSd)) ||
            any(suppliedSd[parameters] <= 0)) {
          yuima.stop("'sd' must give a positive named proposal standard deviation for each Bayes parameter.")
        }
        suppliedSd[parameters]
      }
      stageCenter <- if (is.null(center)) initial else {
        suppliedCenter <- unlist(center)
        if (is.null(names(suppliedCenter)) || any(!parameters %in% names(suppliedCenter))) {
          yuima.stop("'center' must contain every parameter used by an MpCN step.")
        }
        suppliedCenter[parameters]
      }
      temperature <- if (stage == 1L) {
        (observations + 1L)^(2 / (p * rate) - 1)
      } else if (stage == 2L) {
        ((observations + 1L) * h)^(2 / (p * rate) - 1)
      } else 1
      fit <- .adaQmle_mcmc(initial, objective, stageLogPrior, stageLower, stageUpper,
                           mcmc, algorithm, stageCenter, proposalSd, rho, temperature, path)
      estimate <- setNames(fit$par, parameters)
      stageCovariance <- fit$vcov
      stages[[stage]] <- list(stage = stage, target = target,
                              order = stageOrder[[stage]], contrast = contrastName,
                              estimator = "bayes", observations = observations,
                              temperature = temperature, accept.rate = fit$accept.rate,
                              path = fit$path, value = objective(estimate))
    }
    current[parameters] <- estimate
    dimnames(stageCovariance) <- list(parameters, parameters)
    covariance[parameters, parameters] <- stageCovariance
  }

  coefficient <- current[freeParameters]
  finalValue <- tail(vapply(stages, function(x) {
    if (isTRUE(x$skipped)) return(NA_real_)
    if (!is.null(x$value)) return(as.numeric(x$value))
    NA_real_
  }, numeric(1L)), 1L)
  if (is.na(finalValue)) {
    usable <- which(!vapply(stages, function(x) isTRUE(x$skipped), logical(1L)))
    finalValue <- if (length(usable)) stages[[tail(usable, 1L)]]$value else NA_real_
  }
  details <- list(
    p = p,
    order = c(diffusion = symbolic$l0, drift = symbolic$k0),
    contrast.order = stageOrder,
    nhp = numberIncrements * h^p,
    initial.contrast = initial.contrast,
    refinement = refinement,
    expansion = expansion,
    estimator = schedule,
    common.parameters = commonParameters,
    stages = stages
  )
  minuslogl <- function(...) NA_real_
  new("yuima.qmle", call = call, coef = coefficient, fullcoef = current,
      fixed = if (length(fixed)) {
        setNames(as.numeric(unlist(fixed)), names(fixed))
      } else numeric(0),
      vcov = covariance, min = 0.5 * as.numeric(finalValue), details = details,
      minuslogl = minuslogl, method = "adaQmle",
      nobs = as.integer(numberIncrements), model = model)
}
