simulate_multi_particles_with_weights <- function(yuima, nsim, seed, xinits, true_parameter,
    method, sampling) {
    ## :: errors checks

    ## :1: error on yuima model
    if (missing(yuima)) {
        yuima.warn("yuima object is missing.")
        return(NULL)
    }
    if (!inherits(yuima@model, "yuima.state_space_model")) {
        yuima.warn("yuima object must be a stateSpaceModel.")
        return(NULL)
    }
    if (yuima@model@hurst != 0.5) {
        yuima.warn("Hurst parameter must be 0.5 for multi-particle simulation.")
        return(NULL)
    }
    if (length(yuima@model@jump.coeff) != 0) {
        yuima.warn("Model with jump coefficients is not supported.")
        return(NULL)
    }

    if (missing(sampling) || is.null(sampling)) {
        if (is.null(yuima@sampling)) {
            yuima.warn("Sampling object is missing.")
        } else {
            sampling <- yuima@sampling
        }
    }

    sdeModel <- yuima@model
    Terminal <- sampling@Terminal[1]
    Initial <- sampling@Initial[1]

    n <- sampling@n[1]
    r.size <- sdeModel@noise.number
    d.size <- sdeModel@equation.number

    ## :2: error on xinit
    if (missing(xinits) || is.null(xinits)) {
        yuima.warn("Initial values are missing.")
        return(NULL)
    }
    # TODO: check if xinits is a 2-dim array of numerics TODO: consider the
    # type of xinits TODO: maybe nsim is not needed if xinits is not a matrix,
    # error
    if (!is.matrix(xinits)) {
        yuima.warn("Initial values must be a matrix.")
        return(NULL)
    }
    if (missing(nsim) || is.null(nsim)) {
        nsim <- nrow(xinits)
    }

    if (nrow(xinits) != nsim) {
        yuima.warn("Length of xinits missmatches with nsim.")
        return(NULL)
    }
    if (ncol(xinits) != d.size) {
        yuima.warn("Length of items in xinits missmatches with the dimension of the model.")
        return(NULL)
    }


    par.len <- length(sdeModel@parameter@all)

    if (missing(true_parameter) & par.len > 0) {
        true_parameter <- vector(par.len, mode = "list")
        for (i in 1:par.len) {
            true_parameter[[i]] <- 0
        }
        names(true_parameter) <- sdeModel@parameter@all
    }

    yuimaEnv <- new.env()

    if (par.len > 0) {
        for (i in 1:par.len) {
            pars <- sdeModel@parameter@all[i]

            for (j in 1:length(true_parameter)) {
                if (is.na(match(pars, names(true_parameter)[j])) != TRUE) {
                  assign(sdeModel@parameter@all[i], true_parameter[[j]], yuimaEnv)
                }
            }
        }
    }

    ## :: using Euler-Maruyama method
    delta <- (Terminal - Initial)/n
    dW <- rnorm(nsim * n * r.size, 0, sqrt(delta))

    data <- euler_multi_particles_with_weights(xinits, yuima, dW, yuimaEnv)

    return(data)
}
