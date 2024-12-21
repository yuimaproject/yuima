kalmanBucyFilter <- function(yuima,
                             params,
                             mean_init,
                             vcov_init = NULL,
                             delta.vcov.solve = 0.001,
                             are = FALSE,
                             explicit = FALSE,
                             time_homogeneous = FALSE,
                             minuslogl = FALSE,
                             drop_terms = 0,
                             env = globalenv()) {
  # calculate diff of observed variables
  is.observed.equation <- yuima@model@is.observed
  observed.variables <- yuima@model@state.variable[is.observed.equation]
  delta.observed.variable <- array(
    dim = c(
      length(observed.variables),
      length(yuima@data@zoo.data[[1]]) - 1
    ),
    dimnames = list(observed.variables)
  )
  for (variable in observed.variables) {
    delta.observed.variable[variable, ] <- diff(matrix(yuima@data@zoo.data[[which(yuima@model@state.variable == variable)]]))
  }

  # create coefficient matrix of diffusion term of observed/unobserved variables
  col_num <- length(yuima@model@diffusion[[1]])
  is.observed.column <- rep(TRUE, col_num)
  is.unobserved.column <- rep(TRUE, col_num)
  for (i in 1:col_num) {
    is.observed.column[i] <- all(
      yuima@model@is.observed |
        sapply(yuima@model@diffusion, function(row) {
          as.character(row[i]) == "(0)"
        })
    )
    is.unobserved.column[i] <- all(
      !yuima@model@is.observed |
        sapply(yuima@model@diffusion, function(row) {
          as.character(row[i]) == "(0)"
        })
    )
  }
  if (!all(is.observed.column | is.unobserved.column)) {
    yuima.stop("Invalid diffusion matrix. Cannnot divide columns to observed/unobserved.")
  }
  observed.diffusion.expr <- lapply(yuima@model@diffusion[is.observed.equation], function(x) {
    x[is.observed.column]
  })
  unobserved.diffusion.expr <- lapply(yuima@model@diffusion[!is.observed.equation], function(x) {
    x[is.unobserved.column]
  })

  # evaluate observed.diffusion.expr
  tmp.env <- new.env(parent = env)
  for (i in 1:length(params)) {
    assign(names(params)[i], params[[i]], envir = tmp.env)
  }
  if (are | time_homogeneous) {
    observed.diffusion <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.diffusion.expr, tmp.env)[, , 1]
    inv.squared.observed.diffusion <- solve(tcrossprod(observed.diffusion))
    dim(inv.squared.observed.diffusion) <- c(dim(inv.squared.observed.diffusion), 1)
  } else {
    delta_smaller_than <- delta.vcov.solve
    K <- ceiling(yuima@sampling@delta / delta_smaller_than)
    delta <- yuima@sampling@delta / K
    n <- yuima@sampling@n[1] * K
    time.points <- as.matrix((0:n) * delta)
    observed.diffusion <- kalman_bucy_filter_eval_exp(
      observed.diffusion.expr,
      tmp.env,
      yuima@model@time.variable,
      time.points
    )
    inv.squared.observed.diffusion <- calc_inverse_square(observed.diffusion)
  }
  return(
    kalmanBucyFilter.inner(
      yuima,
      delta.observed.variable,
      params,
      inv.squared.observed.diffusion,
      mean_init,
      vcov_init,
      delta.vcov.solve,
      are,
      explicit,
      time_homogeneous,
      minuslogl,
      drop_terms,
      env
    )
  )
}

kalmanBucyFilter.inner <- function(yuima,
                                   delta.observed.variable,
                                   params,
                                   inv.squared.observed.diffusion,
                                   mean_init,
                                   vcov_init = NULL,
                                   delta.vcov.solve = 0.001,
                                   are = FALSE,
                                   explicit = FALSE,
                                   time_homogeneous = FALSE,
                                   minuslogl = FALSE,
                                   drop_terms = 0,
                                   env = globalenv()) {
  # are : flag if use algebraic Riccati equation or not
  # Calculation of `delta.observed.variable` is relatively slow and it can be a bottle neck in parameter estimation. So users can pass the values of delta.observed.variable.
  # Calculation of `inv_sq_ob_diff` is relatively slow and it can be a bottle neck in parameter estimation. So users can pass the values of inv.squared.observed.diffusion.
  # validate input
  if (!inherits(yuima@model, "yuima.linear_state_space_model")) {
    yuima.stop("model must be yuima.linear_state_space_model")
  }

  ## distinguish observed/unobserved variables
  is.observed.equation <- yuima@model@is.observed
  observed.variables <- yuima@model@state.variable[is.observed.equation]
  unobserved.variables <- yuima@model@state.variable[!is.observed.equation]

  if (are) {
    if (!is.null(vcov_init)) {
      yuima.warn("vcov init is given, but ignored because are = TRUE")
    }
    vcov_init <- matrix()
  } else {
    if (is.null(vcov_init)) {
      yuima.stop("'vcov_init' is required if are = FALSE")
    }
    ## if length of observed.variables is 1 and numeric of length 1 is given to vcov_init, convert vcov_init to matrix.
    if (inherits(vcov_init, "numeric")) {
      if (length(vcov_init) != 1) {
        yuima.stop(
          "Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variables) * length(unobserved.variables) or numeric of length 1."
        )
      }
      if (length(unobserved.variables) > 1) {
        yuima.stop(
          "Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variables) * length(unobserved.variables) when length(observed.variables) > 1."
        )
      }
      vcov_init <- matrix(vcov_init)
    } else if (inherits(vcov_init, "matrix")) {
      if (!all(dim(vcov_init) == length(unobserved.variables))) {
        yuima.stop(
          "Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variables) * length(unobserved.variables) or numeric of length 1."
        )
      }
    } else {
      yuima.stop(
        "Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variable) * length(unobserved.variable) or numeric of length 1."
      )
    }
  }

  # create coefficient matrix of diffusion term of observed/unobserved variables
  col_num <- length(yuima@model@diffusion[[1]])
  is.observed.column <- rep(TRUE, col_num)
  is.unobserved.column <- rep(TRUE, col_num)
  for (i in 1:col_num) {
    is.observed.column[i] <- all(
      yuima@model@is.observed |
        sapply(yuima@model@diffusion, function(row) {
          as.character(row[i]) == "(0)"
        })
    )
    is.unobserved.column[i] <- all(
      !yuima@model@is.observed |
        sapply(yuima@model@diffusion, function(row) {
          as.character(row[i]) == "(0)"
        })
    )
  }
  if (!all(is.observed.column | is.unobserved.column)) {
    yuima.stop("Invalid diffusion matrix. Cannnot divide columns to observed/unobserved.")
  }
  observed.diffusion.expr <- lapply(yuima@model@diffusion[is.observed.equation], function(x) {
    x[is.observed.column]
  })
  unobserved.diffusion.expr <- lapply(yuima@model@diffusion[!is.observed.equation], function(x) {
    x[is.unobserved.column]
  })

  # get coefficient matrix of drift term of observed/unobserved variables
  observed.drift.slope.expr <- yuima@model@drift_slope[is.observed.equation]
  observed.drift.intercept.expr <- yuima@model@drift_intercept[is.observed.equation]
  unobserved.drift.slope.expr <- yuima@model@drift_slope[!is.observed.equation]
  unobserved.drift.intercept.expr <- yuima@model@drift_intercept[!is.observed.equation]


  tmp.env <- new.env(parent = env)
  for (i in 1:length(params)) {
    assign(names(params)[i], params[[i]], envir = tmp.env)
  }

  # evaluate each coefficients
  if (are) {
    # coefficients are independent of t.
    unobserved.drift.slope <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.drift.slope.expr, tmp.env)
    unobserved.drift.intercept <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.drift.intercept.expr, tmp.env)
    unobserved.diffusion <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.diffusion.expr, tmp.env)
    observed.drift.slope <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.drift.slope.expr, tmp.env)
    observed.drift.intercept <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.drift.intercept.expr, tmp.env)

    subsamp_rate <- 1
  } else {
    # evaluate in higher frequency for vcov
    delta_smaller_than <- delta.vcov.solve
    subsamp_rate <- ceiling(yuima@sampling@delta / delta_smaller_than)

    if (time_homogeneous) {
      # coefficients are independent of t.
      unobserved.drift.slope <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.drift.slope.expr, tmp.env)
      unobserved.drift.intercept <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.drift.intercept.expr, tmp.env)
      unobserved.diffusion <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.diffusion.expr, tmp.env)
      observed.drift.slope <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.drift.slope.expr, tmp.env)
      observed.drift.intercept <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.drift.intercept.expr, tmp.env)
    } else {
      # coefficients are dependent of t.
      subsamp_delta <- yuima@sampling@delta / subsamp_rate
      subsamp_n <- yuima@sampling@n[1] * subsamp_rate + 1
      time.points <- as.matrix((0:(subsamp_n - 1)) * subsamp_delta)

      unobserved.drift.slope <- kalman_bucy_filter_eval_exp(
        unobserved.drift.slope.expr,
        tmp.env,
        yuima@model@time.variable,
        time.points
      )
      unobserved.drift.intercept <- kalman_bucy_filter_eval_exp(
        unobserved.drift.intercept.expr,
        tmp.env,
        yuima@model@time.variable,
        time.points
      )
      unobserved.diffusion <- kalman_bucy_filter_eval_exp(
        unobserved.diffusion.expr,
        tmp.env,
        yuima@model@time.variable,
        time.points
      )
      observed.drift.slope <- kalman_bucy_filter_eval_exp(
        observed.drift.slope.expr,
        tmp.env,
        yuima@model@time.variable,
        time.points
      )
      observed.drift.intercept <- kalman_bucy_filter_eval_exp(
        observed.drift.intercept.expr,
        tmp.env,
        yuima@model@time.variable,
        time.points
      )
    }
  }
  dim(unobserved.drift.intercept) <- dim(unobserved.drift.intercept)[c(1, 3)]
  dim(observed.drift.intercept) <- dim(observed.drift.intercept)[c(1, 3)]
  filter_res <- calc_kalman_bucy_filter_cpp(
    unobserved.drift.slope,
    unobserved.drift.intercept,
    unobserved.diffusion,
    observed.drift.slope,
    observed.drift.intercept,
    inv.squared.observed.diffusion,
    vcov_init,
    mean_init,
    yuima@sampling@delta,
    delta.observed.variable,
    are,
    explicit,
    time_homogeneous,
    minuslogl,
    drop_terms,
    subsamp_rate
  )
  vcov <- filter_res$vcov
  mean <- filter_res$mean
  minuslogl <- filter_res$minuslogl
  rownames(mean) <- unobserved.variables
  ts.mean <- ts(
    t(mean),
    start = start(yuima@data@zoo.data[[1]]),
    frequency = frequency(yuima@data@zoo.data[[1]])
  )
  res <- new(
    "yuima.kalmanBucyFilter",
    model = yuima@model,
    mean = ts.mean,
    vcov = vcov,
    mean.init = mean_init,
    vcov.init = vcov_init,
    delta = yuima@sampling@delta,
    data = yuima@data,
    minuslogl = minuslogl
  )
  return(res)
}

# evaluate coefficients
# When expr is dependent of t, evaluated values are 3-dim matrix.
kalman_bucy_filter_eval_exp <- function(expr, env, time.variable, time.points) {
  vec <- diffusionTermCpp(expr, time.variable, time.points, env)
  res <- array(vec, dim = c(length(expr), length(expr[[1]]), length(time.points)))
  return(res)
}

# When expr is independent of t, evaluated values are 2-dim matrix.
# Returns 3-dim array with dim(res)[3] == 1.
kalman_bucy_filter_eval_exp_time_homogeneous <- function(expr, env) {
  n_row <- length(expr)
  n_col <- length(expr[[1]])
  res <- array(dim = c(n_row, n_col, 1))
  for (r in 1:n_row) {
    for (c in 1:n_col) {
      res[r, c, 1] <- eval(expr[[r]][c], envir = env)
    }
  }

  return(res)
}
