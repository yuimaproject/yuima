## Distribution Law

setClass("yuima.law",representation(rng = "function",
                                    density = "function",
                                    cdf = "function",
                                    quantile = "function",
                                    characteristic = "function",
                                    param.measure = "character",
                                    time.var = "character",
                                    dim = "numLike")
           )

setMethod("initialize", "yuima.law",
                function(.Object,
                         rng = function(n,...){},
                         density = function(x,...){},
                         cdf = function(q,...){},
                         quantile = function(p,...){},
                         characteristic = function(u,...){},
                         param.measure = character(),
                         time.var = character(),
                         dim = NA
                ){
                  .Object@rng <- rng
                  .Object@density <- density
                  .Object@cdf <- cdf
                  .Object@quantile <- quantile
                  .Object@characteristic <- characteristic
                  .Object@param.measure <- param.measure
                  .Object@time.var <- time.var
                  .Object@dim <- dim
                  return(.Object)
                }
)

setMethod("rand","yuima.law",
          function(object, n, param, ...){
            res <- aux.rand(object, n, param, ...)
            return(res)
          }
)

setMethod("dens","yuima.law",
          function(object, x, param, log = FALSE, ...){
            res <- aux.dens(object, x, param, log,  ...)
            return(res)
          }
)

setMethod("cdf","yuima.law",
          function(object, q, param,  ...){
            res <- aux.cdf(object, q, param, log,  ...)
            return(res)
          }
)

setMethod("quant","yuima.law",
          function(object, p, param,  ...){
            res <- aux.quant(object, p, param,  ...)
            return(res)
          }
)

setMethod("char","yuima.law",
          function(object, u, param,  ...){
            res <- aux.char(object, u, param,  ...)
            return(res)
          }
)


#  Constructor

setLaw <- function(rng = function(n,...){NULL},
                   density = function(x,...){NULL},
                   cdf = function(q,...){NULL},
                   quant = function(p,...){NULL},
                   characteristic = function(u,...){NULL},
                   time.var="t",
                   dim = NA){
  param <- NULL
  param.rng <- extrapParam(myfun = rng, time.var = time.var, aux.var = "n" )
  CondRng<- FALSE
  if(all(param.rng %in% "...")){
  #  yuima.warn("rng is not defined")
  }else{
    CondRng <- TRUE
    param <- param.rng
  }

  param.dens <- extrapParam(myfun = density, time.var = time.var, aux.var = "x" )
  CondDens<- FALSE
  if(all(param.dens %in% "...")){
   # yuima.warn("density is not defined")
  }else{
    CondDens <- TRUE
    param <- param.dens
  }

  if(CondDens){
    if(CondRng){
      if(!all(param.dens %in%  param.rng)){
        yuima.stop("dens and rng  have different parameters")
      }
    }
  }

  param.cdf <- extrapParam(myfun = cdf, time.var = time.var, aux.var = "q" )
  #Condcdf<- FALSE
  if(all(param.cdf %in% "...")){
    #yuima.warn("cdf is not defined")
  }else{
   # Condcdf <- TRUE
    if(is.null(param)){
      param <- param.cdf
    }else{
      if(!all(param %in%  param.cdf)){
        yuima.stop("cdf has different parameters")
      }
    }
  }

  param.quant <- extrapParam(myfun = quant, time.var = time.var, aux.var = "p" )
#  Condquant<- FALSE
  if(all(param.quant %in% "...")){
    #yuima.warn("cdf is not defined")
  }else{
#    Condquant <- TRUE
    if(is.null(param)){
      param <- param.quant
    }else{
      if(!all(param %in%  param.quant)){
        yuima.stop("quantile has different parameters")
      }
    }
  }

  param.char <- extrapParam(myfun = characteristic,
                            time.var = time.var,
                            aux.var = "u" )

  if(all(param.char %in% "...")){
  #  yuima.warn("char is not defined")
  }else{
    if(is.null(param)){
      param <- param.char
    }else{
      if(!all(param %in%  param.char)){
        yuima.stop("quantile has different parameters")
      }
    }
  }
  if(is.null(param)){
    param<-character()
  }

  res <- new("yuima.law",
             rng = rng,
             density = density,
             cdf = cdf,
             characteristic = characteristic,
             quantile = quant,
             param.measure = param,
             time.var = time.var,
             dim = NA)
  return(res)
}

extrapParam <- function(myfun, time.var, aux.var){
  dummy <- names(as.list(args(myfun)))
  dummy <- dummy[-length(dummy)]
  if(dummy[1] != aux.var){
    yuima.stop("Change rand.var or charac.var ...")
  }
  cond <- dummy %in% time.var
  # if(!any(cond)){
  #   yuima.warn("The yuima.law is the distribution
  #              of the jump size, in a CP process")
  # }

  dummy.par <- dummy[!cond]
  dummy <- dummy.par[!dummy.par%in%aux.var]
  return(dummy)
}
