simulate_multi_particles_with_weights <- function(yuima, xinits, true_parameter,
    method, sampling, seed) {
    ## :: errors checks

    ## :1: error on yuima model
    if (missing(yuima)) {
        yuima.warn("yuima object is missing.")
        return(NULL)
    }

    model <- yuima@model
    if (!inherits(model, "yuima.state_space_model")) {
        yuima.warn("yuima object must be a stateSpaceModel.")
        return(NULL)
    }
    if (model@hurst != 0.5) {
        yuima.warn("Hurst parameter must be 0.5 for multi-particle simulation.")
        return(NULL)
    }
    if (length(model@jump.coeff) != 0) {
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


    ## :2: error on xinit
    if (missing(xinits) || is.null(xinits)) {
        yuima.warn("Initial values are missing.")
        return(NULL)
    }
    # TODO: check if xinits is a 2-dim array of numerics
    # TODO: consider the type of xinits
    if (!is.matrix(xinits)) {
        yuima.warn("Initial values must be a matrix.")
        return(NULL)
    }

    d.size <- model@equation.number
    if (ncol(xinits) != d.size) {
        yuima.warn("ncol(xinits) missmatches with the dimension of the model.")
        return(NULL)
    }

    all_parameters <- model@parameter@all
    npar <- length(all_parameters)
    if (missing(true_parameter) & npar > 0) {
        true_parameter <- vector(npar, mode = "list")
        for (i in 1:npar) {
            true_parameter[[i]] <- 0
        }
        names(true_parameter) <- all_parameters
    }

    yuimaEnv <- new.env()

    if (npar > 0) {
        for (i in 1:npar) {
            par_name <- all_parameters[i]

            for (j in 1:length(true_parameter)) {
                if (par_name == names(true_parameter)[j]) {
                  assign(par_name, true_parameter[[j]], yuimaEnv)
                }
            }
        }
    }

    # prepare values of random noise
    delta <- sampling@delta
    nsim <- nrow(xinits)
    n <- sampling@n[1]
    r.size <- model@noise.number

    if (!missing(seed)) {
        set.seed(seed)
    }
    dW <- rnorm(nsim * n * r.size, 0, sqrt(delta))

    ## simulate using Euler-Maruyama method
    data <- euler_multi_particles_with_weights(xinits, yuima, dW, yuimaEnv)

    return(data)
}
