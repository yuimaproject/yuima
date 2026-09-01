#' Continuous-time particle filter for state-space models
#'
#' Implements the continuous-time particle approximation and minimal-variance
#' branching algorithm described in Sections 9.2 and 9.6 of Bain and Crisan
#' (2009). The signal is propagated with the Euler-Maruyama method and the
#' particle weights are updated from all observed dimensions.
#'
#' @param yuima A `yuima` object containing a `yuima.state_space_model` and
#'   observations.
#' @param xinits A numeric matrix whose rows are the initial values of the
#'   unobserved state for each particle.
#' @param params A named list containing all model parameters.
#' @param init A positive integer giving the index of the initial sampling
#'   time. Defaults to 1.
#' @param steps A positive integer giving the number of observation intervals
#'   to process. Defaults to all intervals following `init`.
#' @param simulations_per_weight_update A positive integer giving the number
#'   of Euler-Maruyama steps per observation interval. Defaults to 1.
#' @param weight_updates_per_branching A positive integer giving the number of
#'   weight updates between branching corrections. Defaults to 100.
#' @param store_particles Particle storage mode. `"last"` stores the final
#'   filtering cloud, `"observation"` stores clouds at observation times,
#'   `"all"` additionally stores every Euler-Maruyama state, and `"none"`
#'   stores no individual particles.
#' @param keep_ancestors Logical; whether to store particle-parent indices at
#'   observation times. This requires `store_particles = "observation"` or
#'   `"all"`.
#' @param interval_levels Numeric vector of central weighted particle interval
#'   levels strictly between zero and one. Use `NULL` to omit intervals.
#' @param seed An optional random seed.
#'
#' @return An object of class `yuima.particleFilter`. Filtered means,
#'   variance-covariance matrices, effective sample sizes, branching
#'   indicators, likelihood diagnostics, and requested weighted intervals are
#'   always stored. Individual particles are controlled by `store_particles`.
#'
#' @details
#' The observation covariance is calculated from the observed rows of the
#' model diffusion matrix. It must be positive definite. As in Chapter 9 of
#' Bain and Crisan (2009), the Brownian noises driving the observed and
#' unobserved equations must be independent.
#'
#' For validation with synthetic Euler-Maruyama data, the propagation step
#' used to generate the signal should match the filter propagation step. If
#' `simulations_per_weight_update` is greater than one, generate the synthetic
#' signal on the same finer grid before subsampling its observations. This is
#' especially important for chaotic dynamics such as Lorenz-63.
#'
#' @references
#' Bain, A. and Crisan, D. (2009). Fundamentals of Stochastic Filtering.
#' Springer, Sections 9.2 and 9.6.
#'
#' @examples
#' \dontrun{
#' ## One-dimensional nonlinear double-well signal
#' nonlinear_model <- setModel(
#'   drift = c("a * X - b * X^3", "X + 0.2 * X^3"),
#'   diffusion = matrix(
#'     c("state_sd", "0", "0", "obs_sd"),
#'     nrow = 2, byrow = TRUE
#'   ),
#'   solve.variable = c("X", "Y"),
#'   state.variable = c("X", "Y"),
#'   observed.variable = "Y"
#' )
#' nonlinear_sampling <- setSampling(delta = 0.01, n = 1000)
#' nonlinear_parameters <- list(
#'   a = 1, b = 1, state_sd = 0.5, obs_sd = 0.5
#' )
#' set.seed(1001)
#' nonlinear_object <- simulate(
#'   nonlinear_model,
#'   sampling = nonlinear_sampling,
#'   xinit = c(-1, 0),
#'   true.parameter = nonlinear_parameters
#' )
#' nonlinear_initial_particles <- matrix(
#'   rnorm(1000, mean = -1, sd = 0.4),
#'   ncol = 1
#' )
#' nonlinear_filter <- particleFilter(
#'   nonlinear_object,
#'   xinits = nonlinear_initial_particles,
#'   params = nonlinear_parameters,
#'   simulations_per_weight_update = 1,
#'   weight_updates_per_branching = 10,
#'   store_particles = "observation",
#'   interval_levels = 0.9,
#'   seed = 1002
#' )
#' plot(nonlinear_filter, plot_truth = TRUE, level = 0.9)
#'
#' ## Three-dimensional Lorenz-63 signal with two observations
#' lorenz_model <- setModel(
#'   drift = c(
#'     "sigma * (X2 - X1)",
#'     "X1 * (rho - X3) - X2",
#'     "X1 * X2 - beta * X3",
#'     "X1",
#'     "X3"
#'   ),
#'   diffusion = diag(c(0.10, 0.10, 0.10, 0.60, 0.60)),
#'   solve.variable = c("X1", "X2", "X3", "Y1", "Y2"),
#'   state.variable = c("X1", "X2", "X3", "Y1", "Y2"),
#'   observed.variable = c("Y1", "Y2")
#' )
#' lorenz_sampling <- setSampling(delta = 0.01, n = 1000)
#' lorenz_parameters <- list(sigma = 10, rho = 28, beta = 8 / 3)
#' set.seed(3001)
#' lorenz_object <- simulate(
#'   lorenz_model,
#'   sampling = lorenz_sampling,
#'   xinit = c(1, 1, 1, 0, 0),
#'   true.parameter = lorenz_parameters
#' )
#' lorenz_initial_particles <- cbind(
#'   rnorm(1000, 1, 0.2),
#'   rnorm(1000, 1, 0.2),
#'   rnorm(1000, 1, 0.2)
#' )
#' lorenz_filter <- particleFilter(
#'   lorenz_object,
#'   xinits = lorenz_initial_particles,
#'   params = lorenz_parameters,
#'   simulations_per_weight_update = 1,
#'   weight_updates_per_branching = 8,
#'   store_particles = "observation",
#'   interval_levels = 0.9,
#'   seed = 3002
#' )
#' plot(lorenz_filter, plot_truth = TRUE, level = 0.9)
#' }
#'
#' @export
particleFilter <- function(
    yuima,
    xinits,
    params = list(),
    init = 1L,
    steps = NULL,
    simulations_per_weight_update = 1L,
    weight_updates_per_branching = 100,
    store_particles = c("last", "observation", "all", "none"),
    keep_ancestors = FALSE,
    interval_levels = 0.95,
    seed) {
  filter_call <- match.call()
  if (missing(yuima)) {
    yuima.stop("yuima object is missing.")
  }
  if (!inherits(yuima, "yuima")) {
    yuima.stop("yuima object is not of class yuima.")
  }

  model <- yuima@model
  if (!inherits(model, "yuima.state_space_model")) {
    yuima.stop("yuima@model must be a state space model.")
  }
  if (length(model@hurst) != 1L || model@hurst != 0.5) {
    yuima.stop("Hurst parameter must be 0.5 for particle filtering.")
  }
  if (length(model@jump.coeff) != 0L) {
    yuima.stop("Models with jump coefficients are not supported.")
  }

  is_observed <- model@is.observed
  d_observed <- sum(is_observed)
  d_unobserved <- model@equation.number - d_observed
  if (d_observed < 1L || d_unobserved < 1L) {
    yuima.stop("The model must have at least one observed and one unobserved variable.")
  }

  if (missing(xinits) || is.null(xinits)) {
    yuima.stop("Initial values are missing.")
  }
  if (!is.matrix(xinits) || !is.numeric(xinits)) {
    yuima.stop("xinits must be a numeric matrix.")
  }
  if (nrow(xinits) < 1L || ncol(xinits) != d_unobserved) {
    yuima.stop("xinits must have one row per particle and one column per unobserved variable.")
  }
  if (any(!is.finite(xinits))) {
    yuima.stop("xinits must contain only finite values.")
  }

  grid <- as.numeric(yuima@sampling@grid[[1]])
  if (length(grid) < 2L || any(!is.finite(grid)) || any(diff(grid) <= 0)) {
    yuima.stop("yuima@sampling must contain a strictly increasing finite grid.")
  }

  if (missing(init) || is.null(init)) {
    init <- 1L
  }
  if (length(init) != 1L || !is.numeric(init) || !is.finite(init) ||
      init <= 0 || init %% 1 != 0) {
    yuima.stop("init must be a positive integer.")
  }
  init <- as.integer(init)
  if (init >= length(grid)) {
    yuima.stop("init must identify a sampling time with at least one following interval.")
  }

  if (missing(steps) || is.null(steps)) {
    steps <- length(grid) - init
  }
  if (length(steps) != 1L || !is.numeric(steps) || !is.finite(steps) ||
      steps <= 0 || steps %% 1 != 0) {
    yuima.stop("steps must be a positive integer.")
  }
  steps <- as.integer(steps)
  if (init + steps > length(grid)) {
    yuima.stop("init + steps must not exceed the sampling grid length.")
  }

  interval_deltas <- diff(grid)[seq.int(init, length.out = steps)]
  delta <- interval_deltas[[1L]]
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(delta))
  if (any(abs(interval_deltas - delta) > tolerance)) {
    yuima.stop("Particle filtering currently requires an equidistant sampling grid.")
  }

  all_parameters <- model@parameter@all
  if (missing(params) || is.null(params)) {
    params <- list()
  }
  if (!is.list(params)) {
    yuima.stop("params must be a named list.")
  }
  if (length(all_parameters) != length(params) ||
      !setequal(all_parameters, names(params))) {
    yuima.stop("params must contain exactly all model parameters.")
  }
  if (length(params) > 0L && any(!vapply(params, function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x)
  }, logical(1)))) {
    yuima.stop("Every model parameter must be a finite numeric scalar.")
  }

  validate_positive_integer <- function(x, name, default) {
    if (missing(x) || is.null(x)) {
      return(default)
    }
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
        x <= 0 || x %% 1 != 0) {
      yuima.stop(paste(name, "must be a positive integer."))
    }
    as.integer(x)
  }
  simulations_per_weight_update <- validate_positive_integer(
    simulations_per_weight_update,
    "simulations_per_weight_update",
    1L
  )
  weight_updates_per_branching <- validate_positive_integer(
    weight_updates_per_branching,
    "weight_updates_per_branching",
    100L
  )

  store_particles <- match.arg(store_particles)
  if (length(keep_ancestors) != 1L || !is.logical(keep_ancestors) ||
      is.na(keep_ancestors)) {
    yuima.stop("keep_ancestors must be TRUE or FALSE.")
  }
  if (keep_ancestors &&
      !store_particles %in% c("observation", "all")) {
    yuima.stop(paste(
      "keep_ancestors = TRUE requires store_particles =",
      "\"observation\" or \"all\"."
    ))
  }
  if (is.null(interval_levels)) {
    interval_levels <- numeric()
  }
  if (!is.numeric(interval_levels) || any(!is.finite(interval_levels)) ||
      any(interval_levels <= 0 | interval_levels >= 1)) {
    yuima.stop("interval_levels must lie strictly between zero and one.")
  }
  interval_levels <- sort(unique(as.numeric(interval_levels)))

  diffusion <- model@diffusion
  r_size <- model@noise.number
  diffusion_is_zero <- vapply(diffusion, function(row) {
    vapply(seq_len(r_size), function(column) {
      as.character(row[column]) %in% c("0", "(0)")
    }, logical(1))
  }, logical(r_size))
  if (r_size == 1L) {
    diffusion_is_zero <- matrix(diffusion_is_zero, nrow = 1L)
  }
  nonzero_observed <- rowSums(!diffusion_is_zero[, is_observed, drop = FALSE]) > 0L
  nonzero_unobserved <- rowSums(!diffusion_is_zero[, !is_observed, drop = FALSE]) > 0L
  if (any(nonzero_observed & nonzero_unobserved)) {
    yuima.stop(paste(
      "Observed and unobserved equations must be driven by independent",
      "Brownian-noise columns, as assumed by the Section 9 particle filter."
    ))
  }

  state_variables <- model@state.variable
  observed_variables <- state_variables[is_observed]
  unobserved_variables <- state_variables[!is_observed]
  zoo_data <- yuima@data@zoo.data
  observed_data <- matrix(
    NA_real_,
    nrow = d_observed,
    ncol = steps + 1L,
    dimnames = list(observed_variables, NULL)
  )
  data_indices <- seq.int(init, init + steps)
  for (j in seq_along(observed_variables)) {
    state_index <- match(observed_variables[[j]], state_variables)
    if (state_index > length(zoo_data)) {
      yuima.stop(paste("Observation data are missing for", observed_variables[[j]]))
    }
    series <- as.numeric(zoo_data[[state_index]])
    if (length(series) < max(data_indices)) {
      yuima.stop(paste("Observation data are too short for", observed_variables[[j]]))
    }
    observed_data[j, ] <- series[data_indices]
  }
  if (any(!is.finite(observed_data))) {
    yuima.stop("Observation data must contain only finite values in the selected interval.")
  }
  delta_y <- observed_data[, -1L, drop = FALSE] -
    observed_data[, -ncol(observed_data), drop = FALSE]
  observed_left <- observed_data[, -ncol(observed_data), drop = FALSE]

  env <- list2env(params, parent = baseenv())
  partial_drift <- partial.eval(model@drift, env)
  observed_drift <- partial_drift[is_observed]
  unobserved_drift <- partial_drift[!is_observed]
  observed_diffusion <- partial.eval(unlist(diffusion[is_observed]), env)
  unobserved_diffusion <- partial.eval(unlist(diffusion[!is_observed]), env)

  collapse_expressions <- function(expressions) {
    parse(text = paste0(
      "c(", paste(as.character(expressions), collapse = ","), ")"
    ))
  }
  observed_drift <- collapse_expressions(observed_drift)
  unobserved_drift <- collapse_expressions(unobserved_drift)
  observed_diffusion <- collapse_expressions(observed_diffusion)
  unobserved_diffusion <- collapse_expressions(unobserved_diffusion)

  if (!missing(seed)) {
    if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed)) {
      yuima.stop("seed must be a finite numeric scalar.")
    }
    set.seed(seed)
  }

  storage_code <- match(
    store_particles,
    c("none", "last", "observation", "all")
  ) - 1L
  result <- .Call(
    "_yuima_euler_multi_particles_with_weights_and_branching",
    xinits,
    rep(1 / nrow(xinits), nrow(xinits)),
    grid[[init]],
    as.integer(r_size),
    as.integer(d_observed),
    delta,
    steps,
    model@time.variable,
    unobserved_variables,
    observed_variables,
    simulations_per_weight_update,
    weight_updates_per_branching,
    observed_drift,
    unobserved_drift,
    observed_diffusion,
    unobserved_diffusion,
    delta_y,
    observed_left,
    env,
    as.integer(storage_code),
    keep_ancestors,
    interval_levels,
    PACKAGE = "yuima"
  )

  filter_time <- grid[data_indices]
  particle_time <- if (store_particles == "last") {
    tail(filter_time, 1L)
  } else {
    filter_time
  }
  path_time <- seq(
    from = filter_time[[1L]],
    by = delta / simulations_per_weight_update,
    length.out = steps * simulations_per_weight_update + 1L
  )

  empty_array <- array(numeric(), dim = c(0L, 0L, 0L))
  particles <- if (is.null(result$particles)) empty_array else result$particles
  weights <- if (is.null(result$weights)) {
    matrix(numeric(), nrow = 0L, ncol = 0L)
  } else {
    result$weights
  }
  paths <- if (is.null(result$paths)) empty_array else result$paths
  ancestors <- if (is.null(result$ancestors)) {
    matrix(integer(), nrow = 0L, ncol = 0L)
  } else {
    matrix(as.integer(result$ancestors), nrow = nrow(result$ancestors))
  }

  if (length(particles) > 0L) {
    dimnames(particles) <- list(
      particle = seq_len(nrow(xinits)),
      state = unobserved_variables,
      time = format(particle_time, trim = TRUE)
    )
    dimnames(weights) <- list(
      particle = seq_len(nrow(xinits)),
      time = format(particle_time, trim = TRUE)
    )
  }
  if (length(paths) > 0L) {
    dimnames(paths) <- list(
      particle = seq_len(nrow(xinits)),
      state = unobserved_variables,
      time = format(path_time, trim = TRUE)
    )
  }
  if (length(ancestors) > 0L) {
    dimnames(ancestors) <- list(
      particle = seq_len(nrow(xinits)),
      time = format(filter_time, trim = TRUE)
    )
  }

  filtered_mean <- stats::ts(
    t(result$mean),
    start = filter_time[[1L]],
    frequency = 1 / delta
  )
  colnames(filtered_mean) <- unobserved_variables
  filtered_vcov <- result$vcov
  dimnames(filtered_vcov) <- list(
    state = unobserved_variables,
    state = unobserved_variables,
    time = format(filter_time, trim = TRUE)
  )
  intervals <- if (length(interval_levels) == 0L) {
    array(
      numeric(),
      dim = c(d_unobserved, 2L, length(filter_time), 0L)
    )
  } else {
    result$intervals
  }
  dimnames(intervals) <- list(
    state = unobserved_variables,
    bound = c("lower", "upper"),
    time = format(filter_time, trim = TRUE),
    level = if (length(interval_levels) == 0L) {
      NULL
    } else {
      paste0(format(100 * interval_levels, trim = TRUE), "%")
    }
  )

  new(
    "yuima.particleFilter",
    model = model,
    data = yuima@data,
    mean = filtered_mean,
    vcov = filtered_vcov,
    time = filter_time,
    particles = particles,
    weights = weights,
    paths = paths,
    ancestors = ancestors,
    ess = as.numeric(result$ess),
    branched = as.logical(result$branched),
    intervals = intervals,
    interval.levels = interval_levels,
    logLik = as.numeric(result$logLik),
    logLik.increment = as.numeric(result$logLik_increment),
    call = filter_call,
    settings = list(
      n_particles = nrow(xinits),
      init = init,
      steps = steps,
      simulations_per_weight_update = simulations_per_weight_update,
      weight_updates_per_branching = weight_updates_per_branching,
      store_particles = store_particles,
      keep_ancestors = keep_ancestors
    )
  )
}

#' @rdname particleFilter
#' @export
simulate_multi_particles_with_weights <- function(
    yuima,
    xinits,
    init = 1L,
    steps = NULL,
    params = list(),
    simulations_per_weight_update = 1L,
    weight_updates_per_branching = 100L,
    seed,
    store_particles = c("last", "observation", "all", "none"),
    keep_ancestors = FALSE,
    interval_levels = 0.95) {
  arguments <- list(
    yuima = yuima,
    xinits = xinits,
    params = params,
    init = init,
    steps = steps,
    simulations_per_weight_update = simulations_per_weight_update,
    weight_updates_per_branching = weight_updates_per_branching,
    store_particles = store_particles,
    keep_ancestors = keep_ancestors,
    interval_levels = interval_levels
  )
  if (!missing(seed)) {
    arguments$seed <- seed
  }
  do.call(particleFilter, arguments)
}

setMethod("mean", "yuima.particleFilter", function(x) x@mean)
setMethod("vcov", "yuima.particleFilter", function(object) object@vcov)

#' Plotting Method for Particle Filter
#'
#' Plotting method for objects of class \code{yuima.particleFilter}. The
#' graphical layout, colors, line types, labels, and optional truth overlay
#' follow the plotting method for \code{yuima.kalmanBucyFilter}.
#'
#' @details
#' The shaded interval is the weighted empirical particle interval calculated
#' by \code{particleFilter()}. Therefore, a requested \code{level} must match
#' one of the values supplied through \code{interval_levels} when the filter
#' was run. Otherwise the method raises an error. Set \code{level = 0} to draw
#' no interval.
#'
#' All writable graphical parameters are restored on exit, including when
#' plotting stops with an error.
#'
#' @param x A \code{\link{yuima.particleFilter-class}} object.
#' @param plot_truth Logical. If \code{TRUE}, plot true values of state
#'   variables when they are present in the data.
#' @param level Numeric. If \code{0 < level < 1}, plot the stored weighted
#'   particle interval of that level.
#'
#' @return \code{NULL}, invisibly. A plot is drawn on the active device.
#'
#' @author The YUIMA Project Team
#'
#' @examples
#' \dontrun{
#' model <- setModel(
#'   drift = c("-a * X", "X"),
#'   diffusion = matrix(c("b", "0", "0", "sigma"), 2, 2,
#'                      byrow = TRUE),
#'   solve.variable = c("X", "Y"),
#'   state.variable = c("X", "Y"),
#'   observed.variable = "Y"
#' )
#' sampling <- setSampling(delta = 0.01, n = 100)
#' parameters <- list(a = 0.7, b = 0.3, sigma = 0.4)
#' object <- simulate(
#'   model, sampling = sampling, xinit = c(0, 0),
#'   true.parameter = parameters
#' )
#' result <- particleFilter(
#'   object,
#'   xinits = matrix(rnorm(1000, sd = 0.5), ncol = 1),
#'   params = parameters,
#'   interval_levels = 0.95
#' )
#' plot(result, plot_truth = TRUE, level = 0.95)
#' }
setMethod(
  "plot", "yuima.particleFilter",
  function(x, plot_truth = FALSE, level = 0) {
    if (length(plot_truth) != 1L || !is.logical(plot_truth) ||
        is.na(plot_truth)) {
      yuima.stop("plot_truth must be TRUE or FALSE.")
    }
    if (length(level) != 1L || !is.numeric(level) || !is.finite(level)) {
      yuima.stop("level must be a finite numeric scalar.")
    }

    print_level_interval <- 0 < level && level < 1
    interval_index <- integer()
    if (print_level_interval) {
      tolerance <- sqrt(.Machine$double.eps) * max(1, abs(level))
      interval_index <- which(abs(x@interval.levels - level) <= tolerance)
      if (length(interval_index) == 0L) {
        available <- if (length(x@interval.levels) == 0L) {
          "none"
        } else {
          paste(x@interval.levels, collapse = ", ")
        }
        yuima.stop(paste0(
          "The requested interval level ", level,
          " was not stored. Available interval levels: ", available, "."
        ))
      }
      interval_index <- interval_index[[1L]]
    }

    original_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(original_par), add = TRUE)

    mar <- c(4, 4, 0, 2)
    oma <- c(1, 1, 1, 1)
    mgp <- c(2.5, 1, 0)
    cols <- c(
      "black",
      "blue",
      grDevices::rgb(173, 216, 230, maxColorValue = 255, alpha = 100)
    )
    ltys <- c(1, 2, 1)
    lwd <- c(1, 1, 8)
    upper_margin_coef <- if (print_level_interval || plot_truth) 1.5 else 1.2

    unobserved_variables <-
      x@model@state.variable[!x@model@is.observed]
    n_screen <- length(unobserved_variables)
    graphics::par(
      mfrow = c(n_screen, 1L),
      mar = mar,
      oma = oma,
      mgp = mgp,
      cex = original_par$cex,
      mex = original_par$mex
    )

    time_data <- as.vector(stats::time(x@mean))
    truth_indices <- seq.int(
      x@settings$init,
      length.out = length(time_data)
    )

    for (i in seq_len(n_screen)) {
      var_name <- unobserved_variables[[i]]
      est_x_data <- as.numeric(x@mean[, var_name])
      mean_value <- mean(est_x_data)

      if (print_level_interval) {
        lower_bound <- as.numeric(x@intervals[i, 1L, , interval_index])
        upper_bound <- as.numeric(x@intervals[i, 2L, , interval_index])
        ylim <- c(
          (min(lower_bound) - mean_value) * 1.2 + mean_value,
          (max(upper_bound) - mean_value) * upper_margin_coef + mean_value
        )
      } else {
        ylim <- c(
          (min(est_x_data) - mean_value) * 1.2 + mean_value,
          (max(est_x_data) - mean_value) * upper_margin_coef + mean_value
        )
      }

      if (plot_truth) {
        original_index <- match(var_name, x@model@state.variable)
        true_series <- as.numeric(x@data@zoo.data[[original_index]])
        if (length(true_series) < max(truth_indices)) {
          yuima.stop(paste("Data for", var_name, "are shorter than the filter interval."))
        }
        true_x_data <- true_series[truth_indices]
        if (all(is.na(true_x_data))) {
          yuima.stop(paste("No data for", var_name, "is found."))
        }
        if (any(!is.finite(true_x_data))) {
          yuima.stop(paste("Non-finite data for", var_name, "were found."))
        }
        ylim <- c(
          min(ylim[[1L]], true_x_data * 1.2),
          max(ylim[[2L]], true_x_data * upper_margin_coef)
        )
      }

      graphics::plot(
        0, 0,
        type = "n",
        xlim = range(time_data),
        xlab = "Time",
        ylim = ylim,
        ylab = var_name,
        cex.lab = 1.5,
        cex.axis = 1.2
      )

      if (print_level_interval) {
        graphics::polygon(
          c(time_data, rev(time_data)),
          c(lower_bound, rev(upper_bound)),
          col = cols[[3L]],
          border = NA
        )
      }

      if (plot_truth) {
        graphics::lines(
          time_data, true_x_data,
          col = cols[[1L]], lty = ltys[[1L]]
        )
        estimation_line_style_index <- 2L
      } else {
        estimation_line_style_index <- 1L
      }
      graphics::lines(
        time_data, est_x_data,
        col = cols[[estimation_line_style_index]],
        lty = ltys[[estimation_line_style_index]]
      )

      if (plot_truth) {
        legends <- if (print_level_interval) {
          c(
            paste("ture", var_name),
            "Particle filter",
            paste0(100 * level, "% confidence interval")
          )
        } else {
          c(paste("ture", var_name), "Particle filter")
        }
        graphics::legend(
          "top", legend = legends,
          col = cols, lty = ltys, lwd = lwd
        )
      } else if (print_level_interval) {
        legends <- c(
          "Particle filter",
          paste0(100 * level, "% confidence interval")
        )
        graphics::legend(
          "top", legend = legends,
          col = cols[c(1L, 3L)], lty = ltys[c(1L, 3L)]
        )
      }
    }

    invisible(NULL)
  }
)
