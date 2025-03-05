euler_multi_particles_with_weights <- function(xinits, weight_init, model, sampling,
    dW, deltaY, env) {
    # args for euler method
    is_observed <- model@is.observed
    modelstate <- model@solve.variable[!is_observed]
    modeltime <- model@time.variable
    drift <- model@drift
    diffusion <- model@diffusion
    r_size <- model@noise.number
    initial_time <- sampling@Initial[1]
    n <- sampling@n[1]
    delta <- sampling@delta

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
    observed_partial_evaled_diffusion <- parse(text = paste("matrix(c(", paste(as.character(observed_partial_evaled_diffusion),
        collapse = ","), "), ncol = ", r_size, ", byrow = TRUE)"))
    unobserved_partial_evaled_diffusion <- parse(text = paste("matrix(c(", paste(as.character(unobserved_partial_evaled_diffusion),
        collapse = ","), "), ncol = ", r_size, ", byrow = TRUE)"))

    X_cube <- .Call("_yuima_euler_multi_particles_with_weights", xinits, weight_init,
        initial_time, r_size, delta, n, dW, modeltime, modelstate, observed_partial_evaled_drift,
        unobserved_partial_evaled_drift, observed_partial_evaled_diffusion, unobserved_partial_evaled_diffusion,
        deltaY, env, new.env(), PACKAGE = "yuima")
    return(X_cube)  # TODO: consider the type of return value
}
