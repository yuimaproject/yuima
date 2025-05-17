# yuima: yuima object
# xinits: initial values of unobserved variables, a matrix of size (nsim, d_unob_size)
# init: initial time. interger indicating the index of the initial time in yuima@sampling@grid
# steps: number of steps to simulate, positive integer
# params: a list of parameters
# simulations_per_weight_update: number of simulations per weight update
# seed: seed for random number generation
simulate_multi_particles_with_weights <- function(yuima, xinits, init, steps, params,
    simulations_per_weight_update, seed) {
    # errors checks

    ## :1: error on yuima model
    if (missing(yuima)) {
        yuima.warn("yuima object is missing.")
        return(NULL)
    }
    if (!inherits(yuima, "yuima")) {
        yuima.warn("yuima object is not of class yuima.")
        return(NULL)
    }

    model <- yuima@model
    if (!inherits(model, "yuima.state_space_model")) {
        yuima.warn("yuima@model must be a state space model.")
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
    d_ob_size <- sum(model@is.observed)
    d_unob_size <- d.size - d_ob_size
    if (ncol(xinits) != d_unob_size) {
        yuima.warn("ncol(xinits) missmatches with the number of unobserved variables.")
        return(NULL)
    }
    nsim <- nrow(xinits)
    
    ## :3: error on init
    if (missing(init) || is.null(init)) {
        init <- 1
    }
    if (!is.numeric(init) || init <= 0 || init %% 1 != 0) {
        yuima.warn("init must be a positive integer.")
        return(NULL)
    }
    if (init > length(yuima@sampling@grid[[1]])) {
        yuima.warn("init must be less than or equal to length(yuima@sampling@grid[[1]]).")
        return(NULL)
    }
    
    ## :4: error on steps
    if (missing(steps) || is.null(steps)) {
        steps <- length(yuima@sampling@grid[[1]]) - init 
    }
    if (!is.numeric(steps) || steps <= 0 || steps %% 1 != 0) {
        yuima.warn("steps must be a positive integer.")
        return(NULL)
    }
    if (init + steps > length(yuima@sampling@grid[[1]])) {
        yuima.warn("init + steps must be less than or equal to length(yuima@sampling@grid[[1]]).")
        return(NULL)
    }
    
    ## :5: error on params
    all_parameters <- model@parameter@all
    npar <- length(all_parameters)
    if (npar > 0) {
        if (missing(params) || is.null(params)) {
            yuima.warn("Parameters are missing. Using 0 as default.")
            params <- vector(npar, mode = "list")
            for (i in 1:npar) {
                params[[i]] <- 0
            }
            names(params) <- all_parameters
        }
        if (!is.list(params)) {
            yuima.warn("Parameters must be a list.")
            return(NULL)
        }
        if (length(params) != npar) {
            yuima.warn("Length of params does not match the number of parameters.")
            return(NULL)
        }
        for (i in 1:npar) {
            if (!all_parameters[i] %in% names(params)) {
                yuima.warn(paste("Parameter", all_parameters[i], "is missing."))
                return(NULL)
            }
        }
    }
    ## :6: error on simulations_per_weight_update
    if (missing(simulations_per_weight_update) || is.null(simulations_per_weight_update)) {
        simulations_per_weight_update <- 1
    }
    if (!is.numeric(simulations_per_weight_update) || simulations_per_weight_update <= 0 || simulations_per_weight_update %% 1 != 0) {
        yuima.warn("simulations_per_weight_update must be a positive integer.")
        return(NULL)
    }
    
    # create environment for parameters
    env <- list2env(params)

    sampling <- yuima@sampling
    delta <- sampling@delta
    r_size <- model@noise.number

    if (!missing(seed)) {
        set.seed(seed)
    }

    modelstate <- model@state.variable
    observed_variables <- modelstate[model@is.observed]
    zoodata <- yuima@data@zoo.data
    deltaY <- array(dim = c(length(observed_variables), length(zoodata[[1]]) - 1),
        dimnames = list(observed_variables))
    for (variable in observed_variables) {
        deltaY[variable, ] <- diff(matrix(zoodata[[which(modelstate == variable)]]))
    }
    
    ## simulate using Euler-Maruyama method
    # args for euler method
    weight_init <- rep(1/nsim, nsim)
    is_observed <- model@is.observed
    unobserved_variables <- model@solve.variable[!is_observed] # TODO: rename
    time_variable <- model@time.variable
    drift <- model@drift
    diffusion <- model@diffusion
    #r_size <- model@noise.number
    initial_time <- sampling@grid[[1]][init]
    #n <- sampling@n[1]
    #delta <- sampling@delta

    # partially evaluate drift and diffusion with parameter values
    # TODO: is it possible env contains variable values, not just parameters?
    partial_evaled_drift <- partial.eval(drift, env)
    observed_partial_evaled_drift <- partial_evaled_drift[is_observed]
    unobserved_partial_evaled_drift <- partial_evaled_drift[!is_observed]
    observed_partial_evaled_diffusion <- partial.eval(unlist(diffusion[is_observed]),
                                                      env)
    unobserved_partial_evaled_diffusion <- partial.eval(unlist(diffusion[!is_observed]),
                                                        env)
    
    # convert expressions to length 1
    observed_partial_evaled_drift <- parse(text = paste("c(", paste(as.character(observed_partial_evaled_drift),
                                                                    collapse = ","), ")"))
    unobserved_partial_evaled_drift <- parse(text = paste("c(", paste(as.character(unobserved_partial_evaled_drift),
                                                                      collapse = ","), ")"))
    # observed_partial_evaled_diffusion <- parse(text = paste("matrix(c(", paste(as.character(observed_partial_evaled_diffusion),
    #     collapse = ","), "), ncol = ", r_size, ", byrow = TRUE)"))
    # unobserved_partial_evaled_diffusion <- parse(text = paste("matrix(c(", paste(as.character(unobserved_partial_evaled_diffusion),
    #     collapse = ","), "), ncol = ", r_size, ", byrow = TRUE)"))
    observed_partial_evaled_diffusion <- parse(text = paste("c(", paste(as.character(observed_partial_evaled_diffusion),
                                                                        collapse = ","), ")"))
    unobserved_partial_evaled_diffusion <- parse(text = paste("c(", paste(as.character(unobserved_partial_evaled_diffusion),
                                                                          collapse = ","), ")"))
    data <- .Call("_yuima_euler_multi_particles_with_weights", xinits, weight_init,
                    initial_time, r_size, delta, steps, time_variable, unobserved_variables, simulations_per_weight_update,
                    observed_partial_evaled_drift,
                    unobserved_partial_evaled_drift, observed_partial_evaled_diffusion, unobserved_partial_evaled_diffusion,
                    deltaY, env, PACKAGE = "yuima")
    return(data)
}
