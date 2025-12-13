# QMLE for linear state space model

# main function
qmle.linear_state_space_model <- function(yuima, start, lower = NULL, upper = NULL, method = "L-BFGS-B", fixed = list(), envir = globalenv(), filter_mean_init, explicit = FALSE, drop_terms = 0, ...) {
  # validation
  if (drop_terms >= yuima@sampling@n[[1]] + 1) {
    yuima.stop("`drop_terms` must be smaller than the number of observations (=`yuima@sampleing@n[1] + 1`)")
  }
  
  time.vairable <- yuima@model@time.variable
  if (any(params_in_expr(params = time.vairable, expr = yuima@model@drift))||
      any(params_in_expr(params = time.vairable, expr = yuima@model@diffusion))){
    yuima.stop("Currently, the model must be time-homogeneous.")
  }
  
  observed.variables <- yuima@model@state.variable[yuima@model@is.observed]
  if (any(params_in_expr(params = observed.variables, expr = yuima@model@drift))||
      any(params_in_expr(params = observed.variables, expr = yuima@model@diffusion))){
    yuima.stop("Currently, the model must not include observed variables.")
  }

  #### fixed
  if (length(fixed) > 0 && !is.Poisson(yuima) && !is.CARMA(yuima) && !is.COGARCH(yuima)) {
    fixed_env <- new.env()
    for (params in names(fixed)) {
      assign(params, fixed[[params]], envir = fixed_env)
    }

    drift <- yuima@model@drift
    diffusion <- yuima@model@diffusion
    drift_slope <- yuima@model@drift_slope
    drift_intercept <- yuima@model@drift_intercept

    fixed_yuima <- yuima

    ## fix drift
    fixed_yuima@model@drift <- partial.eval(drift, fixed_env)

    ## fix diffusion
    for (i in 1:length(diffusion)) {
      fixed_yuima@model@diffusion[[i]] <- partial.eval(diffusion[[i]], fixed_env)
    }

    ## fix drift_slope
    for (i in 1:length(drift_slope)) {
      fixed_yuima@model@drift_slope[[i]] <- partial.eval(drift_slope[[i]], fixed_env)
    }

    ## fix drift_intercept
    for (i in 1:length(drift_intercept)) {
      fixed_yuima@model@drift_intercept[[i]] <- partial.eval(drift_intercept[[i]], fixed_env)
    }

    ## remove fixed parameter
    fixed_yuima@model@parameter@all <- yuima@model@parameter@all[!is.element(yuima@model@parameter@all, names(fixed))]
    fixed_yuima@model@parameter@common <- yuima@model@parameter@common[!is.element(yuima@model@parameter@common, names(fixed))]

    fixed_yuima@model@parameter@drift <- yuima@model@parameter@drift[!is.element(yuima@model@parameter@drift, names(fixed))]
    attr(fixed_yuima@model@parameter@drift, "observed") <- attr(yuima@model@parameter@drift, "observed")[!is.element(yuima@model@parameter@drift, names(fixed))]
    attr(fixed_yuima@model@parameter@drift, "unobserved") <- attr(yuima@model@parameter@drift, "unobserved")[!is.element(yuima@model@parameter@drift, names(fixed))]

    fixed_yuima@model@parameter@diffusion <- yuima@model@parameter@diffusion[!is.element(yuima@model@parameter@diffusion, names(fixed))]
    attr(fixed_yuima@model@parameter@diffusion, "observed") <- attr(yuima@model@parameter@diffusion, "observed")[!is.element(yuima@model@parameter@diffusion, names(fixed))]
    attr(fixed_yuima@model@parameter@diffusion, "unobserved") <- attr(yuima@model@parameter@diffusion, "unobserved")[!is.element(yuima@model@parameter@diffusion, names(fixed))]

    fixed_yuima@model@parameter@jump <- yuima@model@parameter@jump[!is.element(yuima@model@parameter@jump, names(fixed))]
    fixed_yuima@model@parameter@measure <- yuima@model@parameter@measure[!is.element(yuima@model@parameter@measure, names(fixed))]

    new.start <- start[!is.element(names(start), names(fixed))]
    new.lower <- lower[!is.element(names(lower), names(fixed))]
    new.upper <- upper[!is.element(names(upper), names(fixed))]

    # execute the function with new arguments
    res <- qmle.linear_state_space_model(fixed_yuima,
      start = new.start,
      method = method, fixed = list(),
      envir = envir, lower = new.lower, upper = new.upper,
      filter_mean_init = filter_mean_init, explicit = explicit,
      drop_terms = drop_terms, ...
    )

    # align object to return
    res@call <- match.call()
    res@model <- yuima@model
    fixed.res <- fixed
    mode(fixed.res) <- "numeric"
    res@fullcoef <- c(res@fullcoef, fixed.res)
    res@fixed <- fixed.res
    return(res)
  }

  # compute delta of observed data
  delta.observed.variable <- array(dim = c(length(observed.variables), yuima@sampling@n[[1]]), dimnames = list(observed.variables))
  for (variable in observed.variables) {
    delta.observed.variable[variable, ] <- diff(matrix(yuima@data@zoo.data[[which(yuima@model@state.variable == variable)]]))
  }

  # split parameters
  theta1_theta2_list <- split_parameters_into_theta1_and_theta2(yuima)

  # estimate theta1 (parameters in observed diffusion)
  theta1 <- theta1_theta2_list$theta1
  if (length(theta1) == 0) {
    estimate_theta1_res <- NULL
    theta1_coef <- NULL
  } else {
    estimate_theta1_res <- estimate.state_space.theta1(yuima, start = start, lower = lower, upper = upper, theta1 = theta1, delta.observed.variable = delta.observed.variable, method = method, envir = envir, drop_terms = drop_terms, ...)
    theta1_coef <- estimate_theta1_res@coef
  }

  # estimate theta2 (parameters in drift and unobserved diffusion)
  theta2 <- theta1_theta2_list$theta2
  if (length(theta2) == 0) {
    estimate_theta2_res <- NULL
  } else {
    estimate_theta2_res <- estimate.state_space.theta2(yuima, start = start, lower = lower, upper = upper, theta2 = theta2, theta1.est = theta1_coef, filter_mean_init = filter_mean_init, delta.observed.variable = delta.observed.variable, explicit = explicit, method = method, envir = envir, drop_terms = drop_terms, ...)
  }

  # align object to return
  call <- match.call()
  if (is.null(estimate_theta1_res)) {
    # no theta1 case
    coef <- estimate_theta2_res@coef
    vcov_names <- rownames(estimate_theta2_res@vcov)
    vcov <- estimate_theta2_res@vcov
    rownames(vcov) <- vcov_names
    colnames(vcov) <- vcov_names
    optim_min <- c(theta1 = NA, theta2 = estimate_theta2_res@min)
    optim_details <- list(theta1 = NA, theta2 = estimate_theta2_res@details)
    minuslogl <- function(p) {
      is.theta1 <- is.element(names(p), theta1)
      minuslogl <- list(theta1 = NA, theta2 = estimate_theta2_res@minuslogl(p = p[!is.theta1]))
      return(minuslogl)
    }
    optim_method <- c(theta1 = NA, theta2 = estimate_theta2_res@method)
  } else if (is.null(estimate_theta2_res)) {
    # no theta2 case
    coef <- estimate_theta1_res@coef
    vcov_names <- rownames(estimate_theta1_res@vcov)
    vcov <- estimate_theta1_res@vcov
    rownames(vcov) <- vcov_names
    colnames(vcov) <- vcov_names
    optim_min <- c(theta1 = estimate_theta1_res@min, theta2 = NA)
    optim_details <- list(theta1 = estimate_theta1_res@details, theta2 = NA)
    minuslogl <- function(p) {
      is.theta1 <- is.element(names(p), theta1)
      minuslogl <- list(theta1 = estimate_theta1_res@minuslogl(p = p[is.theta1]), theta2 = NA)
      return(minuslogl)
    }
    optim_method <- c(theta1 = estimate_theta1_res@method, theta2 = NA)
  } else {
    coef <- c(estimate_theta1_res@coef, estimate_theta2_res@coef)
    vcov_names <- c(rownames(estimate_theta1_res@vcov), rownames(estimate_theta2_res@vcov))
    vcov <- as.matrix(bdiag(list(estimate_theta1_res@vcov, estimate_theta2_res@vcov)))
    rownames(vcov) <- vcov_names
    colnames(vcov) <- vcov_names
    optim_min <- c(theta1 = estimate_theta1_res@min, theta2 = estimate_theta2_res@min)
    optim_details <- list(theta1 = estimate_theta1_res@details, theta2 = estimate_theta2_res@details)
    minuslogl <- function(p) {
      is.theta1 <- is.element(names(p), theta1)
      minuslogl <- list(theta1 = estimate_theta1_res@minuslogl(p = p[is.theta1]), theta2 = estimate_theta2_res@minuslogl(p = p[!is.theta1]))
      return(minuslogl)
    }
    optim_method <- c(theta1 = estimate_theta1_res@method, theta2 = estimate_theta2_res@method)
  }
  yuima.nobs <- as.integer(max(unlist(lapply(get.zoo.data(yuima), length)) - 1, na.rm = TRUE))
  res <- new("yuima.linear_state_space_qmle",
    call = call,
    coef = coef,
    fullcoef = coef,
    vcov = vcov,
    min = optim_min,
    details = optim_details,
    minuslogl = minuslogl,
    method = optim_method,
    nobs = yuima.nobs,
    model = yuima@model,
    drop_terms = drop_terms,
    explicit = explicit,
    mean_init = filter_mean_init
  )
  return(res)
}

# function to calculate minus quasi-log likelihood for theta1 estimation
# yuima: yuima objects
# theta1: values of parameters
# env: must contain objects below
# drop_terms: losses for first given terms will be ingnored on calculation
minuslogl.linear_state_space.theta1 <- function(yuima, theta1, env, drop_terms, delta.observed.variable, h) {
  # evaluate the diffusion of the observation equation
  tmp.env <- new.env(parent = env)
  for (i in 1:length(theta1)) {
    assign(names(theta1)[i], theta1[[i]], envir = tmp.env)
  }
  DIFFUSION <- yuima@model@diffusion[yuima@model@is.observed]
  nrow <- length(DIFFUSION)
  ncol <- length(DIFFUSION[[1]])
  observed.diffusion <- matrix(nrow = nrow, ncol = ncol)
  for (r in 1:nrow) {
    for (c in 1:ncol) {
      observed.diffusion[r, c] <- eval(DIFFUSION[[r]][c], envir = tmp.env)
    }
  }

  sq.observed.diffusion <- tcrossprod(as.matrix(observed.diffusion))
  inv.sq.observed.diffusion <- solve(sq.observed.diffusion)
  logdet.sq.observed.diffusion <- log(det(sq.observed.diffusion))

  # calculate likelihood
  QL <- minusloglcpp_linear_state_space_theta1(logdet.sq.observed.diffusion, inv.sq.observed.diffusion, delta.observed.variable, h)
  return(-drop(QL))
}

# estimate theta1
estimate.state_space.theta1 <- function(yuima, start, method, envir = globalenv(),
                                        lower, upper, theta1, delta.observed.variable, drop_terms, ...) {
  # validate arguments
  if (missing(yuima)) {
    yuima.stop("yuima object is missing.")
  }

  ## FIXME: maybe we should choose initial values at random within lower/upper
  ##        at present, qmle stops
  if (missing(start)) {
    yuima.stop("Starting values for the parameters are missing.")
  }

  yuima.nobs <- as.integer(max(unlist(lapply(get.zoo.data(yuima), length)) - 1, na.rm = TRUE))

  diff.par <- yuima@model@parameter@diffusion

  nm <- names(start)

  idx.diff <- match(diff.par, nm)

  new.start <- start[idx.diff]
  new.upper <- upper[nm[idx.diff]]
  new.lower <- lower[nm[idx.diff]]

  call <- match.call()

  ###################################
  # Solve theta1 explicitly
  # if observed equation is 1-dim
  #   and observed diffusion parameter is 1-dim
  #   and observed diffusion coef is parameter itself.
  ###################################
  # specify observed columns in diffusion matrix
  col_num <- length(yuima@model@diffusion[[1]])
  is.observed.column <- rep(TRUE, col_num)
  is.unobserved.column <- rep(TRUE, col_num)
  for (i in 1:col_num) {
    is.observed.column[i] <- all(yuima@model@is.observed | sapply(yuima@model@diffusion, function(row) as.character(row[i]) == "(0)"))
    is.unobserved.column[i] <- all(!yuima@model@is.observed | sapply(yuima@model@diffusion, function(row) as.character(row[i]) == "(0)"))
  }
  if (!all(is.observed.column | is.unobserved.column)) {
    yuima.stop("Invalid diffusion matrix. Cannnot divide columns to observed/unobserved.")
  }
  observed.diffusion <- lapply(yuima@model@diffusion[yuima@model@is.observed], function(x) x[is.observed.column])[[1]]
  h <- yuima@sampling@delta
  if (sum(yuima@model@is.observed) == 1 && sum(is.observed.column) == 1) {
    n <- yuima@sampling@n[[1]]
    param_name <- yuima@model@parameter@diffusion[!attr(yuima@model@parameter@diffusion, "unobserved")]
    if (as.character(observed.diffusion) == paste("(", param_name, ")", sep = "")) {
      dX <- delta.observed.variable
      coef <- sqrt(sum(dX^2) / (h * n))
      names(coef) <- param_name
      vcov <- matrix(coef^2 / (2 * n))
      rownames(vcov) <- param_name
      f <- function(p) {
        mycoef <- as.list(p)
        names(mycoef) <- c(param_name)

        return((-sum(dX^2) / (2 * h * mycoef[[param_name]])^2) + n * log(mycoef[[param_name]]))
      }

      res <- new("yuima.qmle",
        call = call,
        coef = coef,
        fullcoef = coef,
        vcov = vcov,
        min = f(coef),
        details = list(), # slot for return value of optim. no value because calculating explicitly.
        minuslogl = f,
        method = "explicit formula",
        nobs = yuima.nobs,
        model = yuima@model
      )
      return(res)
    }
  }

  # set args for optim
  ## define objective function
  f <- function(p) {
    mycoef <- as.list(p)
    names(mycoef) <- nm[idx.diff]
    return(minuslogl.linear_state_space.theta1(yuima, mycoef, envir, drop_terms = drop_terms, delta.observed.variable, h))
  }

  ## args for optim
  mydots <- as.list(call)[-(1:2)] # list of given args except `yuima` and `start`
  ### set proper values for `mydots`
  mydots$fn <- as.name("f")
  mydots$par <- unlist(new.start)
  mydots$hessian <- TRUE
  mydots$upper <- as.numeric(unlist(new.upper))
  mydots$lower <- as.numeric(unlist(new.lower))
  ### remove unnecessary params from `mydots`
  mydots$start <- NULL
  mydots$envir <- NULL
  mydots$theta1 <- NULL
  mydots$delta.observed.variable <- NULL
  mydots$drop_terms <- NULL

  # optimization
  if ((length(mydots$par) > 1) | any(is.infinite(c(mydots$upper, mydots$lower)))) {
    mydots$method <- method
    oout <- do.call(optim, args = mydots)
  } else {
    # use `optimize`
    ## set proper values for `mydots`
    mydots$f <- mydots$fn
    mydots$interval <- as.numeric(c(unlist(new.lower), unlist(new.upper)))

    ## remove unnecessary params from `mydots`
    mydots$fn <- NULL
    mydots$par <- NULL
    mydots$hessian <- NULL
    mydots$lower <- NULL
    mydots$upper <- NULL
    mydots$method <- NULL
    mydots$envir <- NULL

    ## optimization
    opt <- do.call(optimize, args = mydots)
    theta <- opt$minimum
    names(theta) <- diff.par

    ## calculate Hessian
    mydots$fn <- mydots$f
    mydots$par <- theta
    mydots$f <- NULL
    mydots$interval <- NULL
    HESS <- do.call(optimHess, args = mydots)

    oout <- list(
      par = theta, value = opt$objective,
      hessian = HESS
    )
  }

  coef <- oout$par[theta1] # estimation

  # calculate vcov
  hessian <- oout$hessian[theta1, theta1, drop = FALSE]
  vcov <- matrix(NA, length(coef), length(coef))
  if (length(coef)) {
    rrr <- try(solve(hessian), TRUE)
    if (class(rrr)[1] != "try-error") {
      vcov <- rrr
    }
  }

  dummycov <- matrix(0, length(coef), length(coef))
  rownames(dummycov) <- names(coef)
  colnames(dummycov) <- names(coef)
  dummycov[rownames(vcov), colnames(vcov)] <- vcov
  vcov <- dummycov

  # align return value and return
  final_res <- new("yuima.qmle",
    call = call,
    coef = coef,
    fullcoef = coef,
    vcov = vcov,
    min = oout$value,
    details = oout,
    minuslogl = f,
    method = method,
    nobs = yuima.nobs,
    model = yuima@model
  )

  return(final_res)
}

# function to calculate minus quasi-log likelihood for theta2 estimation
# yuima: yuima objects
# inv.squared.observed.diffusion: inverse of squared observed diffusion coef
# theta2: params in drift or unobserved diffusion
# filter_mean_init: estimated value for X_0
# env: env to evaluate drift and diffusion

minuslogl.linear_state_space.theta2 <- function(yuima, delta.observed.variable, inv.squared.observed.diffusion, theta2, filter_mean_init, env, explicit, drop_terms, h) {
  is.observed <- yuima@model@is.observed

  # define env for eval
  tmp.env <- new.env(parent = env)
  if (length(theta2) > 0) {
    for (i in 1:length(theta2)) {
      assign(names(theta2)[i], theta2[[i]], envir = tmp.env)
    }
  }

  # calculate `m` (estimation of `x`) using filter
  are <- TRUE # NOTE: Estimation using are=FALSE is not implemented yet.
  filter_res <- kalmanBucyFilter.inner(yuima, delta.observed.variable = delta.observed.variable, params = theta2, inv.squared.observed.diffusion = inv.squared.observed.diffusion, mean_init = filter_mean_init, are = are, explicit = explicit, minuslogl = TRUE, drop_terms = drop_terms, env = tmp.env)

  return(filter_res$minuslogl)
}

# estimate theta2
estimate.state_space.theta2 <- function(yuima, start, method = "L-BFGS-B", envir = globalenv(),
                                        lower, upper, theta2, theta1.est, delta.observed.variable, filter_mean_init, explicit, drop_terms, ...) {
  if (missing(yuima)) {
    yuima.stop("yuima object is missing.")
  }

  ## FIXME: maybe we should choose initial values at random within lower/upper
  ##        at present, qmle stops
  if (missing(start)) {
    yuima.stop("Starting values for the parameters are missing.")
  }

  yuima.nobs <- as.integer(max(unlist(lapply(get.zoo.data(yuima), length)) - 1, na.rm = TRUE))

  # env to evaluate drift/diffusion of given SDE
  env <- envir

  nm <- names(start)

  idx.theta2 <- match(theta2, nm)

  new.start <- start[idx.theta2]
  new.upper <- upper[nm[idx.theta2]]
  new.lower <- lower[nm[idx.theta2]]

  # calculate observed diffusion
  tmp.env <- new.env(parent = env)
  if (length(theta1.est) > 0) {
    for (i in 1:length(theta1.est)) {
      assign(names(theta1.est)[i], theta1.est[[i]], envir = tmp.env)
    }
  }
  DIFFUSION <- yuima@model@diffusion[yuima@model@is.observed]
  nrow <- length(DIFFUSION)
  ncol <- length(DIFFUSION[[1]])
  observed.diffusion.matrix <- matrix(nrow = nrow, ncol = ncol)
  for (r in 1:nrow) {
    for (c in 1:ncol) {
      observed.diffusion.matrix[r, c] <- eval(DIFFUSION[[r]][c], envir = tmp.env)
    }
  }
  inv.squared.observed.diffusion <- solve(tcrossprod(observed.diffusion.matrix))
  dim(inv.squared.observed.diffusion) <- c(dim(inv.squared.observed.diffusion), 1)
  h <- yuima@sampling@delta

  # set args for optim
  ## define objective function
  f <- function(p) {
    theta2.values <- as.list(p)
    names(theta2.values) <- nm[idx.theta2]
    return(minuslogl.linear_state_space.theta2(yuima, delta.observed.variable = delta.observed.variable, inv.squared.observed.diffusion = inv.squared.observed.diffusion, theta2 = theta2.values, filter_mean_init = filter_mean_init, env = tmp.env, explicit = explicit, drop_terms = drop_terms, h = h))
  }

  call <- match.call()

  ## args for optim
  mydots <- as.list(call)[-(1:2)] # list of given args except `yuima` and `start`
  ### set proper values for `mydots`
  mydots$fn <- as.name("f")
  mydots$par <- unlist(new.start)
  mydots$hessian <- TRUE
  if (method == "L-BFGS-B" | method == "Brent") {
    mydots$upper <- as.numeric(unlist(new.upper))
    mydots$lower <- as.numeric(unlist(new.lower))
  } else {
    mydots$upper <- NULL
    mydots$lower <- NULL
  }


  ### remove unnecessary params from `mydots`
  mydots$start <- NULL
  mydots$envir <- NULL
  mydots$theta1.est <- NULL
  mydots$theta2 <- NULL
  mydots$filter_mean_init <- NULL
  mydots$delta.observed.variable <- NULL
  mydots$explicit <- NULL
  mydots$drop_terms <- NULL

  # optimization
  if ((length(mydots$par) > 1) | any(is.infinite(c(mydots$upper, mydots$lower)))) {
    mydots$method <- method
    oout <- do.call(optim, args = mydots)
  } else {
    # use `optimize`
    ## set proper values for `mydots`
    mydots$f <- mydots$fn
    mydots$interval <- as.numeric(c(unlist(new.lower), unlist(new.upper)))

    ## remove unnecessary params from `mydots`
    mydots$fn <- NULL
    mydots$par <- NULL
    mydots$hessian <- NULL
    mydots$lower <- NULL
    mydots$upper <- NULL
    mydots$method <- NULL
    mydots$envir <- NULL

    ## optimization
    opt <- do.call(optimize, args = mydots)
    theta <- opt$minimum
    names(theta) <- theta2

    ## calculate Hessian
    mydots$fn <- mydots$f
    mydots$par <- theta
    mydots$f <- NULL
    mydots$interval <- NULL
    HESS <- do.call(optimHess, args = mydots)

    oout <- list(
      par = theta, value = opt$objective,
      hessian = HESS
    )
  }

  coef <- oout$par # estimation

  # calculate vcov
  vcov <- matrix(NA, length(coef), length(coef))
  if (length(coef)) {
    rrr <- try(solve(oout$hessian), TRUE)
    if (class(rrr)[1] != "try-error") {
      vcov <- rrr
    }
  }

  dummycov <- matrix(0, length(coef), length(coef))
  rownames(dummycov) <- names(coef)
  colnames(dummycov) <- names(coef)
  dummycov[rownames(vcov), colnames(vcov)] <- vcov
  vcov <- dummycov

  # align return value and return
  final_res <- new("yuima.qmle",
    call = call,
    coef = coef,
    fullcoef = coef,
    vcov = vcov,
    min = oout$value,
    details = oout,
    minuslogl = f,
    method = method,
    nobs = yuima.nobs,
    model = yuima@model
  )

  return(final_res)
}

split_parameters_into_theta1_and_theta2 <- function(yuima) {
  # args: yuima: yuima object
  # return: list(theta1="character", theta2="character")
  #   theta1: params in observed.diffusion
  #   theta2: params in drift, unobserved.diffusion
  # if there are at least 1 parameter included in both theta1 and 2, throw error

  theta1 <- yuima@model@parameter@diffusion[attr(yuima@model@parameter@diffusion, "observed")]
  theta2 <- unique(c(yuima@model@parameter@drift, yuima@model@parameter@diffusion[attr(yuima@model@parameter@diffusion, "unobserved")]))

  if (any(duplicated(c(theta1, theta2)))) {
    #yuima.stop("Cannot split parameters into theta1 and theta2.")
    theta2 <- theta2[!is.element(theta2,theta1)]
  }
  return(list(theta1 = theta1, theta2 = theta2))
}
