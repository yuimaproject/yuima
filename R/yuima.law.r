## Distribution Law

setClass("yuima.law",representation(rng = "function",
                                    density = "function",
                                    cdf = "function",
                                    param.measure = "character",
                                    characteristic = "function",
                                    time.var = "character",
                                    rand.var = "character",
                                    charact.var = "character",
                                    dim = "numLike")
           )

setMethod("initialize", "yuima.law",
                function(.Object,
                         rng = function(x,...){},
                         density = function(x,...){},
                         cdf = function(x,...){},
                         param.measure = character(),
                         characteristic = function(u,...){},
                         time.var = character(),
                         rand.var = character(),
                         charact.var = character(),
                         dim = NA
                ){
                  .Object@rng <- rng
                  .Object@density <- density
                  .Object@cdf <- cdf
                  .Object@param.measure <- param.measure
                  .Object@characteristic <- characteristic
                  .Object@time.var <- time.var
                  .Object@rand.var <- rand.var
                  .Object@charact.var <- charact.var
                  .Object@dim <- dim
                  return(.Object)
                }
)

#  Constructor

setLaw <- function(rng = function(x,...){NULL},
                   density = function(x,...){NULL},
                   cdf = function(x,...){NULL},
                   characteristic = function(u,...){NULL},
                   time.var="t",
                   rand.var = "x",
                   character.var = "u",
                   dim = NA){

  param.rng <- extrapParam(myfun = rng, time.var = time.var, aux.var = rand.var )
  CondRng<- FALSE
  if(all(param.rng %in% "...")){
    yuima.warn("rng is not defined")
  }else{
    CondRng <- TRUE
    param <- param.rng
  }

  param.dens <- extrapParam(myfun = density, time.var = time.var, aux.var = rand.var )
  CondDens<- FALSE
  if(all(param.dens %in% "...")){
    yuima.warn("density is not defined")
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

  param.cdf <- extrapParam(myfun = cdf, time.var = time.var, aux.var = rand.var )
  Condcdf<- FALSE
  if(all(param.cdf %in% "...")){
    yuima.warn("cdf is not defined")
  }else{
    Condcdf <- TRUE
  }

  if(Condcdf){
    if(!all(param %in%  param.cdf)){
      yuima.stop("cdf has different parameters")
    }
  }

  param.char <- extrapParam(myfun = characteristic,
                            time.var = time.var,
                            aux.var = character.var )

  Condchar<- FALSE
  if(all(param.char %in% "...")){
    yuima.warn("char is not defined")
  }else{
    Condchar <- TRUE
  }

  if(Condchar){
    if(!all(param %in%  param.char)){
      yuima.stop("char has different parameters")
    }
  }

  res <- new("yuima.law",
             rng = rng,
             density = density,
             cdf = cdf,
             characteristic = characteristic,
             param.measure = param,
             time.var = time.var,
             rand.var = rand.var,
             charact.var = character.var,
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
