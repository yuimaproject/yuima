euler_multi_particles_with_weights <- function(xinits, model, sampling,
    dW, env) {
    # args for euler method
    modelstate <- model@solve.variable
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
    partial_evaled_diffusion <- partial.eval(unlist(diffusion), env)

    X_cube <- .Call("_yuima_euler_multi_particles_with_weights", xinits,
        initial_time, r_size, delta, n, dW, modeltime, modelstate,
        partial_evaled_drift, partial_evaled_diffusion, env, new.env(),
        PACKAGE = "yuima")

    return(X_cube)  # TODO: consider the type of return value
}
