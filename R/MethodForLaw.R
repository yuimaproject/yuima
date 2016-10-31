aux.funForLaw <- function(object, param, my.env, dummy){
  # param <- unlist(param)
  name.par <- names(param)
  if(length(object@time.var)>=1){
    if(object@time.var%in%name.par){
      assign(object@time.var,param[[object@time.var]], envir = my.env)
    }else{
      yuima.stop("time.var is not assigned")
    }
    param <- param[!(name.par %in% object@time.var)]
    name.par <- names(param)
  }
  if(length(param)>0){
    if(length(param)!=length(object@param.measure)){
      yuima.stop("mismatch arguments")
    }
    for(i in c(1: length(param))){
      cond<-object@param.measure %in% name.par[[i]]
      assign(object@param.measure[cond], param[[i]], envir = my.env)
    }
  }
  res <- eval(parse(text=dummy),envir=my.env)
  return(res)
}


# rng
aux.rand<- function(object, n, param, ...){
  dummy <- deparse(object@rng)[1]
  dummy <- gsub("function ", deparse(substitute(object@rng)), dummy)
  my.env <- new.env()
  assign("n", n, envir = my.env)
  res <- aux.funForLaw(object, param, my.env, dummy)
  return(res)
  # param <- unlist(param)
  # name.par <- names(param)
  # if(length(object@time.var)>=1){
  #   if(object@time.var%in%name.par){
  #     assign(object@time.var,param[object@time.var], envir = my.env)
  #   }else{
  #     yuima.stop("time.var is not assigned")
  #   }
  #   param <- param[!(name.par %in% object@time.var)]
  #   name.par <- names(param)
  # }
  # if(length(param)>0){
  #   if(length(param)!=length(object@param.measure)){
  #     yuima.stop("mismatch arguments")
  #   }
  #   for(i in c(1: length(param))){
  #     cond<-object@param.measure %in% name.par[i]
  #     assign(object@param.measure[cond], param[i], envir = my.env)
  #   }
  # }
  # res <- eval(parse(text=dummy),envir=my.env)
  # return(res)
}


setGeneric("rand",
           function(object, n, param, ...){
             standardGeneric("rand")
           }
)



# dens

setGeneric("dens",
           function(object, x, param, log = FALSE, ...)
             standardGeneric("dens")
)




aux.dens<- function(object, x, param, log, ...){
  dummy <- deparse(object@density)[1]
  dummy <- gsub("function ", deparse(substitute(object@density)), dummy)
  my.env <- new.env()
  assign("x", x, envir = my.env)
  res <- aux.funForLaw(object, param, my.env, dummy)

  if(log){res <- log(res)}
  return(as.numeric(res))
}

# CDF

# cdf

setGeneric("cdf",
           function(object, q, param, ...)
             standardGeneric("cdf")
)




aux.cdf<- function(object, q, param, ...){
  dummy <- deparse(object@cdf)[1]
  dummy <- gsub("function ", deparse(substitute(object@cdf)), dummy)
  my.env <- new.env()
  assign("q", q, envir = my.env)
  res <- aux.funForLaw(object, param, my.env, dummy)
  return(as.numeric(res))
}

# quantile

setGeneric("quant",
           function(object, p, param, ...)
             standardGeneric("quant")
)




aux.quant <- function(object, p, param, ...){
  dummy <- deparse(object@quantile)[1]
  dummy <- gsub("function ", deparse(substitute(object@quantile)), dummy)
  my.env <- new.env()
  assign("p", p, envir = my.env)
  res <- aux.funForLaw(object, param, my.env, dummy)
  return(as.numeric(res))
}

# characteristic

setGeneric("char",
           function(object, u, param, ...)
             standardGeneric("char")
)




aux.char <- function(object, u, param, ...){
  dummy <- deparse(object@characteristic)[1]
  dummy <- gsub("function ", deparse(substitute(object@characteristic)), dummy)
  my.env <- new.env()
  assign("u", u, envir = my.env)
  res <- aux.funForLaw(object, param, my.env, dummy)
  return(as.numeric(res))
}
