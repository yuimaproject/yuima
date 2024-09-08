setMethod(
  "initialize", "yuima.state_space_model",
  function(.Object,
           drift = expression(),
           diffusion = list(),
           hurst = 0.5,
           jump.coeff = list(),
           measure = list(),
           measure.type = character(),
           parameter = new("model.parameter"),
           state.variable = "x",
           jump.variable = "z",
           time.variable = "t",
           noise.number = numeric(),
           equation.number = numeric(),
           dimension = numeric(),
           solve.variable = character(),
           xinit = expression(),
           J.flag = logical(),
           is.observed = logical()) {
    .Object@drift <- drift
    .Object@diffusion <- diffusion
    .Object@hurst <- hurst
    .Object@jump.coeff <- jump.coeff
    .Object@measure <- measure
    .Object@measure.type <- measure.type
    .Object@parameter <- parameter
    .Object@state.variable <- state.variable
    .Object@jump.variable <- jump.variable
    .Object@time.variable <- time.variable
    .Object@noise.number <- noise.number
    .Object@equation.number <- equation.number
    .Object@dimension <- dimension
    .Object@solve.variable <- solve.variable
    .Object@xinit <- xinit
    .Object@J.flag <- J.flag
    .Object@is.observed <- is.observed
    return(.Object)
  }
)

## setStateSpaceModel
## setter of class 'yuima.model'
## set yuima model from SDE
setStateSpaceModel <- function(drift = NULL,
                               diffusion = NULL,
                               hurst = 0.5,
                               jump.coeff = NULL,
                               measure = list(),
                               measure.type = character(),
                               state.variable = "x",
                               jump.variable = "z",
                               time.variable = "t",
                               solve.variable,
                               xinit = NULL,
                               observed.variable = NULL,
                               unobserved.variable = NULL) {
  # validation
  if (is.null(observed.variable) && is.null(unobserved.variable)) {
    yuima.stop("Either observed or unobserved variable must be given.")
  }

  if (!is.null(observed.variable) && !is.null(unobserved.variable)) {
    all.variable <- c(observed.variable, unobserved.variable)
    if (!all(is.element(all.variable, state.variable)) || !all(is.element(state.variable, all.variable))) {
      yuima.stop("observed and unobserved variables are not consistent with the state variables.")
    }
    if (length(all.variable) != length(state.variable)) {
      yuima.stop("observed and unobserved variables may have duplicated elements.")
    }
  }

  if (is.null(observed.variable) && !is.null(unobserved.variable)) {
    if (!all(is.element(unobserved.variable, state.variable))) {
      yuima.stop("invalid variable name is given in unobserved.varialbe.")
    } else {
      observed.variable <- state.variable[!is.element(state.variable, unobserved.variable)]
    }
  }
  if (!all(is.element(observed.variable, state.variable))) {
    yuima.stop("invalid variable name is given in observed.varialbe.")
  }
  model <- setModel(
    drift = drift,
    diffusion = diffusion,
    hurst = hurst,
    jump.coeff = jump.coeff,
    measure = measure,
    measure.type = measure.type,
    state.variable = state.variable,
    jump.variable = jump.variable,
    time.variable = time.variable,
    solve.variable = solve.variable,
    xinit = xinit
  )

  state_space_model <- new("yuima.state_space_model",
    drift = model@drift,
    diffusion = model@diffusion,
    hurst = model@hurst,
    jump.coeff = model@jump.coeff,
    measure = model@measure,
    measure.type = model@measure.type,
    parameter = model@parameter,
    state.variable = model@state.variable,
    jump.variable = model@jump.variable,
    time.variable = model@time.variable,
    noise.number = model@noise.number,
    equation.number = model@equation.number,
    dimension = model@dimension,
    solve.variable = model@solve.variable,
    xinit = model@xinit,
    J.flag = model@J.flag,
    is.observed = is.element(model@state.variable, observed.variable)
  )

  return(state_space_model)
}

setMethod(
  "initialize", "yuima.linear_state_space_model",
  function(.Object,
           drift = expression(),
           diffusion = list(),
           hurst = 0.5,
           jump.coeff = list(),
           measure = list(),
           measure.type = character(),
           parameter = new("model.parameter"),
           state.variable = "x",
           jump.variable = "z",
           time.variable = "t",
           noise.number = numeric(),
           equation.number = numeric(),
           dimension = numeric(),
           solve.variable = character(),
           xinit = expression(),
           J.flag = logical(),
           is.observed = logical(),
           drift_slope,
           drift_intercept) {
    .Object@drift <- drift
    .Object@diffusion <- diffusion
    .Object@hurst <- hurst
    .Object@jump.coeff <- jump.coeff
    .Object@measure <- measure
    .Object@measure.type <- measure.type
    .Object@parameter <- parameter
    .Object@state.variable <- state.variable
    .Object@jump.variable <- jump.variable
    .Object@time.variable <- time.variable
    .Object@noise.number <- noise.number
    .Object@equation.number <- equation.number
    .Object@dimension <- dimension
    .Object@solve.variable <- solve.variable
    .Object@xinit <- xinit
    .Object@J.flag <- J.flag
    .Object@is.observed <- is.observed
    .Object@drift_slope <- drift_slope
    .Object@drift_intercept <- drift_intercept
    return(.Object)
  }
)


setLinearStateSpaceModel <- function(drift = NULL,
                                     diffusion = NULL,
                                     hurst = 0.5,
                                     jump.coeff = NULL,
                                     measure = list(),
                                     measure.type = character(),
                                     state.variable = "x",
                                     jump.variable = "z",
                                     time.variable = "t",
                                     solve.variable,
                                     xinit = NULL,
                                     observed.variable = NULL,
                                     unobserved.variable = NULL) {
  model <- setStateSpaceModel(
    drift = drift,
    diffusion = diffusion,
    hurst = hurst,
    jump.coeff = jump.coeff,
    measure = measure,
    measure.type = measure.type,
    state.variable = state.variable,
    jump.variable = jump.variable,
    time.variable = time.variable,
    solve.variable = solve.variable,
    xinit = xinit,
    observed.variable = observed.variable,
    unobserved.variable = unobserved.variable
  )

  ## distinguish observed/unobserved.variable
  unobserved.variable <- model@state.variable[!model@is.observed]

  ## distingish each params in observed/unobserved
  attr(model@parameter@drift, "observed") <- params_in_expr(model@parameter@drift, model@drift[model@is.observed])
  attr(model@parameter@drift, "unobserved") <- params_in_expr(model@parameter@drift, model@drift[!model@is.observed])
  attr(model@parameter@diffusion, "observed") <- params_in_exprs(model@parameter@diffusion, model@diffusion[model@is.observed])
  attr(model@parameter@diffusion, "unobserved") <- params_in_exprs(model@parameter@diffusion, model@diffusion[!model@is.observed])

  # set coefficient matrix of drift term
  tmp.env <- new.env()
  eqnum <- model@equation.number

  ## get intercept
  drift.intercept <- list()
  for (i in 1:eqnum) {
    drift.intercept[[i]] <- expression(0)
  }
  for (k in 1:length(unobserved.variable)) {
    assign(unobserved.variable[k], 0, envir = tmp.env)
  }
  for (i in 1:eqnum) {
    drift.intercept[[i]] <- partial.eval(model@drift[i], tmp.env)
  }

  ## get slope
  drift.slope <- list()
  for (i in 1:eqnum) {
    drift.slope[[i]] <- rep(expression(0), length(unobserved.variable))
  }
  for (j in 1:length(unobserved.variable)) {
    for (k in 1:length(unobserved.variable)) {
      assign(unobserved.variable[k], as.numeric(k == j), envir = tmp.env)
    }
    for (i in 1:eqnum) {
      slope_and_intercept <- partial.eval(model@drift[i], tmp.env)
      slope.call <- call("-", as.call(slope_and_intercept)[[1]], as.call(drift.intercept[[i]])[[1]])
      slope <- as.expression(slope.call)
      drift.slope[[i]][j] <- slope
    }
  }
  rm(tmp.env)

  linear_state_space_model <- new("yuima.linear_state_space_model",
    drift = model@drift,
    diffusion = model@diffusion,
    hurst = model@hurst,
    jump.coeff = model@jump.coeff,
    measure = model@measure,
    measure.type = model@measure.type,
    parameter = model@parameter,
    state.variable = model@state.variable,
    jump.variable = model@jump.variable,
    time.variable = model@time.variable,
    noise.number = model@noise.number,
    equation.number = model@equation.number,
    dimension = model@dimension,
    solve.variable = model@solve.variable,
    xinit = model@xinit,
    J.flag = model@J.flag,
    is.observed = is.element(model@state.variable, observed.variable),
    drift_slope = drift.slope,
    drift_intercept = drift.intercept
  )
  return(linear_state_space_model)
}
