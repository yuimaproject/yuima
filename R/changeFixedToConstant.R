## 7/8/2021 Kito
## Some functions didn't work when 'fixed' parameters are given. (at present, qmle and adaBayes)
## We added functions below to allow functions above to work with fixed params by rewriting yuima model. 
transform.drift <- function(drift, fixed, unfixed.nms, state.nms ,index, unfixed.vals, state.vals){## Kurisaki 4/10/2021
  drift <- drift
  fixed <- fixed
  unfixed.nms <- unfixed.nms
  state.nms <- state.nms
  unfixed.vals <-  unfixed.vals
  state.vals <- state.vals
  
  ## substitute for the fixed parameters
  for(var in names(fixed)) {
      assign(var, fixed[[var]])
  }
  
  ## substitute for the unfixed parameters
  i <- 1
  for(nm in unfixed.nms) {
      assign(nm, unfixed.vals[i])
      i <- i+1
  }

  ## substitute for the state parameters
  i <- 1
  for(nm in state.nms) {
    assign(nm, state.vals[,i])
    i <- i+1
  }
  
  eval(drift[index])
}

transform.diffusion <- function(diffusion, fixed, unfixed.nms, state.nms, row, index, unfixed.vals, state.vals){
  diffusion <- diffusion
  fixed <- fixed
  unfixed.nms <- unfixed.nms
  state.nms <- state.nms
  unfixed.vals <-  unfixed.vals
  state.vals <- state.vals
  
  ## substitute for the fixed parameters
  for(var in names(fixed)) {
      assign(var, fixed[[var]])
  }
  
  ## substitute x for the unfixed parameters
  i <- 1
  for(nm in unfixed.nms) {
      assign(nm, unfixed.vals[i])
      i <- i+1
  }
  ## substitute for the state parameters
  i <- 1
  for(nm in state.nms) {
    assign(nm, state.vals[,i])
    i <- i+1
  }
  
  eval(diffusion[[row]][index])
}

transform.jump <- function(jump, fixed, unfixed.nms, state.nms, row, index, unfixed.vals, state.vals){
  jump <- jump
  fixed <- fixed
  unfixed.nms <- unfixed.nms
  state.nms <- state.nms
  unfixed.vals <-  unfixed.vals
  state.vals <- state.vals
  
  ## substitute for the fixed parameters
  for(var in names(fixed)) {
      assign(var, fixed[[var]])
  }
  
  ## substitute for the unfixed parameters
  i <- 1
  for(nm in unfixed.nms) {
      assign(nm, unfixed.vals[i])
      i <- i+1
  }
  ## substitute for the state parameters
  i <- 1
  for(nm in state.nms) {
    assign(nm, state.vals[,i])
    i <- i+1
  }
  
  eval(jump[[row]][index])
}

changeFixedParametersToConstant <- function(yuima, fixed) {
    env <- new.env() # environment to calculate estimation
    
    yuima = yuima
    fixed = fixed
    
    # list of names of unfixed parameters
    drift.unfixed.nms = yuima@model@parameter@drift[!is.element(yuima@model@parameter@drift, names(fixed))]
    diffusion.unfixed.nms = yuima@model@parameter@diffusion[!is.element(yuima@model@parameter@diffusion, names(fixed))]

    # list of names of state variables
    state.nms = yuima@model@state.variable
    
    # arguments for new.drift.func & new.diffusion.func
    state.vals = paste("matrix(c(", paste(state.nms, collapse=", "), "), ncol=", length(state.nms),")")
    drift.unfixed.vals = paste("c(", paste(drift.unfixed.nms, collapse=", "), ")")
    diffusion.unfixed.vals = paste("c(", paste(diffusion.unfixed.nms, collapse=", "), ")")
    
    new.drift.func <- function(index, unfixed.vals, state.vals){transform.drift(yuima@model@drift, fixed, drift.unfixed.nms, state.nms, index, unfixed.vals, state.vals)}
    new.diffusion.func <- function(row, index, unfixed.vals, state.vals){transform.diffusion(yuima@model@diffusion, fixed, diffusion.unfixed.nms, state.nms, row, index, unfixed.vals, state.vals)}

    # new drift and diffusion term
    transformed.drift <- paste("new.drift.func(", 1:length(yuima@model@drift), ", ", drift.unfixed.vals, ", ", state.vals, ")", sep="")

    vector.diffusion <- c()
    for(row in 1:length(yuima@model@diffusion)) {
      vector.diffusion <- c(vector.diffusion, paste("new.diffusion.func(", row, ", ", 1:length(yuima@model@diffusion[[row]]), ", ", diffusion.unfixed.vals, ", ", state.vals, ")", sep=""))
    }
    transformed.diffusion <- matrix(vector.diffusion, nrow = length(yuima@model@diffusion),byrow=T)

    # new jump term
    if(length(yuima@model@jump.coeff) > 0) {
      jump.unfixed.nms <- yuima@model@parameter@jump[!is.element(yuima@model@parameter@jump, names(fixed))]
      jump.unfixed.vals = paste("c(", paste(jump.unfixed.nms, collapse=", "), ")")
      new.jump.func <- function(row, index, unfixed.vals, state.vals){transform.jump(yuima@model@jump.coeff, fixed, jump.unfixed.nms, state.nms, row, index, unfixed.vals, state.vals)}
      vector.jump <- c()
      for(row in 1:length(yuima@model@jump.coeff)) {
        vector.jump <- c(vector.jump, paste("new.jump.func(", row, ", ", 1:length(yuima@model@jump.coeff[[row]]), ", ", jump.unfixed.vals, ", ", state.vals, ")", sep=""))
      }
      transformed.jump <- matrix(vector.jump, nrow = length(yuima@model@jump.coeff),byrow=T)

      new.measure = list()
      measure.params <- yuima@model@parameter@measure
      df.measure.params <- measure.params
      if(yuima@model@measure.type=="CP"){
        intensity <- as.character(yuima@model@measure$intensity)
        if(is.element(intensity, names(fixed))) {
          df.measure.params <- df.measure.params[-which(df.measure.params %in% intensity)]
          intensity = as.character(fixed[intensity])
        }
        new.measure[["intensity"]] <- intensity
      }
      df <- yuima@model@measure$df
      expr <- df$expr
      cal <- as.call(expr[[1]])
      params <- as.list(cal[-1])
      for(i in 1:length(params)) {
        if(is.element(c(params[[i]]), names(fixed))) {
          params[[i]] <- fixed[[params[[i]]]]
        }
      }
      cal[-1] <- as.call(params)
      new.measure[["df"]] <- list(as.character(as.expression(cal)))
    } else {
      transformed.jump <- NULL
      new.measure <- list()
    }

    new.ymodel <- setModel(drift = transformed.drift, diffusion = transformed.diffusion, hurst = yuima@model@hurst, 
      jump.coeff = transformed.jump, measure = new.measure, measure.type = yuima@model@measure.type, 
      state.variable = yuima@model@state.variable, jump.variable = yuima@model@jump.variable, 
      time.variable = yuima@model@time.variable, solve.variable = yuima@model@solve.variable, 
      xinit = yuima@model@xinit)
    new.yuima <- setYuima(data = yuima@data, model = new.ymodel, sampling = yuima@sampling, characteristic = yuima@characteristic, functional = yuima@functional)
    return(list(new.yuima=new.yuima, env=env))
}