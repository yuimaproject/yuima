euler_multi_particles_with_weights <- function(xinits, model, sampling,
    dW, env) {

    modelstate <- model@solve.variable
    modeltime <- model@time.variable
    V0 <- model@drift
    V <- model@diffusion
    r.size <- model@noise.number
    Initial <- sampling@Initial[1]
    n <- sampling@n[1]

    ## :: set time step
    delta <- sampling@delta

    ## :: using Euler-Maruyama method
    partial_evaled_drift <- partial.eval(V0, env)
    # TODO: is it possible env contains variable values, not just parameters?
    partial_evaled_diffusion <- partial.eval(unlist(V), env)

    X_cube <- .Call("_yuima_euler_multi_particles_with_weights",
        xinits, Initial, r.size, delta, n, dW, modeltime, modelstate,
        partial_evaled_drift, partial_evaled_diffusion, env, new.env(),
        PACKAGE = "yuima")

    return(X_cube)  # TODO: consider the type of return value
}
