kalmanBucyFilter <- function(yuima, params, mean_init, vcov_init = NULL, delta.vcov.solve = 0.001, are = FALSE,
    explicit = FALSE, time_homogeneous = FALSE, env = globalenv()) {
    # calculate diff of observed variables
    is.observed.equation <- yuima@model@is.observed
    observed.variables <- yuima@model@state.variable[is.observed.equation]
    delta.observed.variable <- array(dim = c(length(observed.variables), length(yuima@data@zoo.data[[1]]) -
        1), dimnames = list(observed.variables))
    for (variable in observed.variables) {
        delta.observed.variable[variable, ] <- diff(matrix(yuima@data@zoo.data[[which(yuima@model@state.variable ==
            variable)]]))
    }

    # create coefficient matrix of diffusion term of observed/unobserved variables
    col_num <- length(yuima@model@diffusion[[1]])
    is.observed.column <- rep(TRUE, col_num)
    is.unobserved.column <- rep(TRUE, col_num)
    for (i in 1:col_num) {
        is.observed.column[i] <- all(yuima@model@is.observed | sapply(yuima@model@diffusion, function(row) {
            as.character(row[i]) == "(0)"
        }))
        is.unobserved.column[i] <- all(!yuima@model@is.observed | sapply(yuima@model@diffusion, function(row) {
            as.character(row[i]) == "(0)"
        }))
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
        observed.diffusion <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.diffusion.expr,
            tmp.env)[, , 1]
        inv.squared.observed.diffusion <- solve(tcrossprod(observed.diffusion))
        dim(inv.squared.observed.diffusion) <- c(dim(inv.squared.observed.diffusion), 1)
    } else {
        delta_smaller_than <- delta.vcov.solve
        K <- ceiling(yuima@sampling@delta/delta_smaller_than)
        delta <- yuima@sampling@delta/K
        n <- yuima@sampling@n[1] * K
        time.points <- as.matrix((0:n) * delta)
        observed.diffusion <- kalman_bucy_filter_eval_exp(observed.diffusion.expr, tmp.env, yuima@model@time.variable,
            time.points)
        inv.squared.observed.diffusion <- calc_inverse_square(observed.diffusion)
    }
    return(kalmanBucyFilter.inner(yuima, delta.observed.variable, params, inv.squared.observed.diffusion,
        mean_init, vcov_init, delta.vcov.solve, are, explicit, time_homogeneous, minuslogl = FALSE,
        drop_terms = 0, env)$filter_res)
}

kalmanBucyFilter.inner <- function(yuima, delta.observed.variable, params, inv.squared.observed.diffusion,
    mean_init, vcov_init = NULL, delta.vcov.solve = 0.001, are = FALSE, explicit = FALSE, time_homogeneous = FALSE,
    minuslogl = FALSE, drop_terms = 0, env = globalenv()) {
    # are : flag if use algebraic Riccati equation or not Calculation of `delta.observed.variable`
    # is relatively slow and it can be a bottle neck in parameter estimation. So users can pass
    # the values of delta.observed.variable.  Calculation of `inv_sq_ob_diff` is relatively slow
    # and it can be a bottle neck in parameter estimation. So users can pass the values of
    # inv.squared.observed.diffusion.  validate input
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
        ## if length of observed.variables is 1 and numeric of length 1 is given to vcov_init,
        ## convert vcov_init to matrix.
        if (inherits(vcov_init, "numeric")) {
            if (length(vcov_init) != 1) {
                yuima.stop("Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variables) * length(unobserved.variables) or numeric of length 1.")
            }
            if (length(unobserved.variables) > 1) {
                yuima.stop("Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variables) * length(unobserved.variables) when length(observed.variables) > 1.")
            }
            vcov_init <- matrix(vcov_init)
        } else if (inherits(vcov_init, "matrix")) {
            if (!all(dim(vcov_init) == length(unobserved.variables))) {
                yuima.stop("Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variables) * length(unobserved.variables) or numeric of length 1.")
            }
        } else {
            yuima.stop("Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variable) * length(unobserved.variable) or numeric of length 1.")
        }
    }

    # create coefficient matrix of diffusion term of observed/unobserved variables
    col_num <- length(yuima@model@diffusion[[1]])
    is.observed.column <- rep(TRUE, col_num)
    is.unobserved.column <- rep(TRUE, col_num)
    for (i in 1:col_num) {
        is.observed.column[i] <- all(yuima@model@is.observed | sapply(yuima@model@diffusion, function(row) {
            as.character(row[i]) == "(0)"
        }))
        is.unobserved.column[i] <- all(!yuima@model@is.observed | sapply(yuima@model@diffusion, function(row) {
            as.character(row[i]) == "(0)"
        }))
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
        unobserved.drift.slope <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.drift.slope.expr,
            tmp.env)
        unobserved.drift.intercept <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.drift.intercept.expr,
            tmp.env)
        unobserved.diffusion <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.diffusion.expr,
            tmp.env)
        observed.drift.slope <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.drift.slope.expr,
            tmp.env)
        observed.drift.intercept <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.drift.intercept.expr,
            tmp.env)

        subsamp_rate <- 1
    } else {
        # evaluate in higher frequency for vcov
        delta_smaller_than <- delta.vcov.solve
        subsamp_rate <- ceiling(yuima@sampling@delta/delta_smaller_than)

        if (time_homogeneous) {
            # coefficients are independent of t.
            unobserved.drift.slope <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.drift.slope.expr,
                tmp.env)
            unobserved.drift.intercept <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.drift.intercept.expr,
                tmp.env)
            unobserved.diffusion <- kalman_bucy_filter_eval_exp_time_homogeneous(unobserved.diffusion.expr,
                tmp.env)
            observed.drift.slope <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.drift.slope.expr,
                tmp.env)
            observed.drift.intercept <- kalman_bucy_filter_eval_exp_time_homogeneous(observed.drift.intercept.expr,
                tmp.env)
        } else {
            # coefficients are dependent of t.
            subsamp_delta <- yuima@sampling@delta/subsamp_rate
            subsamp_n <- yuima@sampling@n[1] * subsamp_rate + 1
            time.points <- as.matrix((0:(subsamp_n - 1)) * subsamp_delta)

            unobserved.drift.slope <- kalman_bucy_filter_eval_exp(unobserved.drift.slope.expr, tmp.env,
                yuima@model@time.variable, time.points)
            unobserved.drift.intercept <- kalman_bucy_filter_eval_exp(unobserved.drift.intercept.expr,
                tmp.env, yuima@model@time.variable, time.points)
            unobserved.diffusion <- kalman_bucy_filter_eval_exp(unobserved.diffusion.expr, tmp.env,
                yuima@model@time.variable, time.points)
            observed.drift.slope <- kalman_bucy_filter_eval_exp(observed.drift.slope.expr, tmp.env,
                yuima@model@time.variable, time.points)
            observed.drift.intercept <- kalman_bucy_filter_eval_exp(observed.drift.intercept.expr, tmp.env,
                yuima@model@time.variable, time.points)
        }
    }
    dim(unobserved.drift.intercept) <- dim(unobserved.drift.intercept)[c(1, 3)]
    dim(observed.drift.intercept) <- dim(observed.drift.intercept)[c(1, 3)]
    filter_res <- calc_kalman_bucy_filter_cpp(unobserved.drift.slope, unobserved.drift.intercept, unobserved.diffusion,
        observed.drift.slope, observed.drift.intercept, inv.squared.observed.diffusion, vcov_init, mean_init,
        yuima@sampling@delta, delta.observed.variable, are, explicit, time_homogeneous, minuslogl, drop_terms,
        subsamp_rate)
    vcov <- filter_res$vcov
    mean <- filter_res$mean
    minuslogl <- filter_res$minuslogl
    rownames(mean) <- unobserved.variables
    ts.mean <- ts(t(mean), start = start(yuima@data@zoo.data[[1]]), frequency = frequency(yuima@data@zoo.data[[1]]))
    filter_res <- new("yuima.kalmanBucyFilter", model = yuima@model, mean = ts.mean, vcov = vcov, mean.init = mean_init,
        vcov.init = vcov_init, delta = yuima@sampling@delta, data = yuima@data)
    return(list(filter_res = filter_res, minuslogl = minuslogl))
}

# evaluate coefficients When expr is dependent of t, evaluated values are 3-dim matrix.
kalman_bucy_filter_eval_exp <- function(expr, env, time.variable, time.points) {
    vec <- diffusionTermCpp(expr, time.variable, time.points, env)
    res <- array(vec, dim = c(length(expr), length(expr[[1]]), length(time.points)))
    return(res)
}

# When expr is independent of t, evaluated values are 2-dim matrix.  Returns 3-dim array with
# dim(res)[3] == 1.
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

setGeneric("mean")
setMethod("mean", "yuima.kalmanBucyFilter", function(x) x@mean)

setGeneric("vcov")
setMethod("vcov", "yuima.kalmanBucyFilter", function(object) object@vcov)

#' Plotting Method for Kalman-Bucy Filter
#'
#' Plotting method for objects of class \code{yuima.kalmanBucyFilter}.
#' 
#' @details This method plots the estimated values of state variables by Kalman-Bucy filter. 
#' Optionally, it can plot true values of state variables which may exist when using simulated data.
#' Also, it can plot confidence interval of \code{level} if 0 < \code{level} < 1.
#' 
#' @param x A \code{\link{yuima.kalmanBucyFilter-class}} object.
#' @param plot_truth Logical. If \code{TRUE}, plot true values of state variables.
#' @param level Numeric. If 0 < \code{level} < 1, plot confidence interval of \code{level}.
#' 
#' @return NULL (plot is drawn)
#' 
#' @author The YUIMA Project Team
#' 
#' @examples
#' \dontrun{
#' # create Kalman-Bucy filter object
#' drift <- c('a*X', 'c*X')
#' diffusion <- matrix(
#'   c('b', '0', '0', 'sigma'),
#'   2, 2
#' )
#' vars <- c('X', 'Y')
#' mod <- setModel(
#'   drift = drift, diffusion = diffusion, solve.variable = vars,
#'   state.variable = vars, observed.variable = 'Y', xinit = c(0, 0)
#' )
#' samp <- setSampling(delta = 0.01, n = 10^3)
#' trueparam <- list(a = -1.5, b = 0.3, c = 1, sigma = 0.02)
#' 
#' ### simulate
#' yuima <- simulate(mod, sampling = samp, true.parameter = trueparam)
#' res <- kalmanBucyFilter(
#'   yuima,
#'   params = trueparam, mean_init = 0, vcov_init = 0.1,
#'   delta.vcov.solve = 0.001, are = FALSE
#' )
#' 
#' ### visualize
#' plot(res, plot_truth = TRUE, level = 0.95)
#' }
setMethod("plot", "yuima.kalmanBucyFilter", function(x, plot_truth = FALSE, level = 0) {
    orig_par <- par(no.readonly = TRUE)
    # config
    mar = c(4, 4, 0, 2)
    oma = c(1, 1, 1, 1)
    mgp = c(2.5, 1, 0)
    cols <- c("black", "blue", rgb(173, 216, 230, maxColorValue = 255, alpha=100))
    ltys <- c(1, 2, 1)
    lwd <- c(1, 1, 8)
    
    print_level_interval = 0 < level && level < 1
    upper_margin_coef = 1.2
    if (print_level_interval || plot_truth) {
        upper_margin_coef = 1.5
    }

    # spitting screen
    unobserved_variables = x@model@state.variable[!x@model@is.observed]

    n_screen <- length(unobserved_variables)
    split.screen(figs = c(n_screen, 1))


    for (i in 1:n_screen) {
        screen(i)
        par(mar = mar, oma = oma, mgp = mgp)
        var_name = unobserved_variables[i]
        # create frame
        est.x.data <- x@mean[, var_name]
        mean_value = mean(est.x.data)
        if (print_level_interval) {
            lower_coef = qnorm((1 - level)/2)
            upper_coef = qnorm(1 - (1 - level)/2)
            lower_bound = est.x.data + lower_coef * sqrt(x@vcov[i, i, ])
            upper_bound = est.x.data + upper_coef * sqrt(x@vcov[i, i, ])
            ylim <- c((min(lower_bound) - mean_value) * 1.2 + mean_value, (max(upper_bound) - mean_value) *
                upper_margin_coef + mean_value)
        } else {
            ylim <- c((min(est.x.data) - mean_value) * 1.2 + mean_value, (max(est.x.data) - mean_value) *
                upper_margin_coef + mean_value)
        }
        if (plot_truth) {
            orig_index = match(var_name, x@model@state.variable)
            true.x.data <- as.vector(x@data@zoo.data[[orig_index]])
            if (all(is.na(true.x.data))) {
                yuima.stop(paste("No data for", var_name, "is found."))
            }
            ylim <- c(min(ylim[1], (true.x.data * 1.2)), max(ylim[2], (true.x.data * upper_margin_coef)))
        }
        time.data <- as.vector(time(x@mean))
        xlim <- c(min(time.data), max(time.data))
        plot(0, 0, type = "n", xlim = xlim, xlab = "Time", ylim = ylim, ylab = paste(var_name), cex.lab = 1.5,
            cex.axis = 1.2)
        
        if (print_level_interval) {
          #lines(time.data, lower_bound, col = cols[3], lty = ltys[3])
          #lines(time.data, upper_bound, col = cols[3], lty = ltys[3])
          polygon(c(time.data, rev(time.data)), c(lower_bound, rev(upper_bound)), col = cols[3], border = NA)
        }

        if (plot_truth) {
            # draw true X line
            lines(as.vector(time(x@mean)), true.x.data, col = cols[1], lty = ltys[1])
            estimation_line_style_index = 2
        } else {
            estimation_line_style_index = 1
        }

        # draw estimation line
        lines(as.vector(time(x@mean)), est.x.data, col = cols[estimation_line_style_index], lty = ltys[estimation_line_style_index])

        # legend
        if (plot_truth) {
            if (print_level_interval) {
                legends <- c(paste("ture", var_name), "Kalman-Bucy filter", paste0(100 * level, "% confidence interval"))
            } else {
                legends <- c(paste("ture", var_name), "Kalman-Bucy filter")
            }
            legend("top", legend = legends, col = cols, lty = ltys, lwd = lwd)
        } else if (print_level_interval) {
            legends <- c("Kalman-Bucy filter", paste0(100 * level, "% confidence interval"))
            legend("top", legend = legends, col = cols[c(1, 3)], lty = ltys[c(1, 3)])
        }
    }
    close.screen(all.screens = TRUE)
    par(orig_par)
})
