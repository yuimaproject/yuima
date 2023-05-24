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

### From a Characteristic Function To yuima.law object

InternalDensity <- function(x, param, mynames, time.names, time.var, 
                            up, low, N_grid, method , myfun, N_Fourier){
  myenv <- new.env()
  for(i in c(1:length(param))){
    assign(mynames[i], param[i] , envir =myenv )
  }
  assign(time.names, time.var, envir =myenv)
  x_old <- x
  x_new <-unique(sort(c(low, x_old, up)))
  x_new <- unique(sort(c(seq(min(x_new)-0.1,max(x_new)+0.1, length.out =N_grid+1),x)))
  alim <- min(x_new)
  blim <- max(x_new)
  
  
  i <- 0:(N_Fourier - 1)
  dx <- (blim - alim)/N_Fourier
  x <- alim + i * dx
  dt <- 2 * pi/(N_Fourier * dx)
  c <- -N_Fourier/2 * dt
  d <- N_Fourier/2 * dt
  u <- c + i * dt
  
  assign("u",u,myenv)
  #dummyphy <- eval(parse(text=myfun))
  #phi <- as.numeric(eval(parse(text=myfun), envir =myenv))
  phi <- eval(parse(text=myfun), envir =myenv)
  #plot(u,phi, type="l")
  
  X <- exp(-(0 + (0+1i)) * i * dt * alim) * phi
  Y <- fft(X)
  density <- dt/(2 * pi) * exp(-(0 + (0+1i)) * c * x) * Y
  invFFT<-data.frame(i = i, u = u, characteristic_function = phi, 
                     x = x, density = Re(density))
  
  #dens <- na.approx(zoo(x=invFFT$density, order.by= invFFT$x), xout=x_old)
  na <- is.na(invFFT$density)
  dens <- approx(x=invFFT$x[!na], y=invFFT$density[!na], xout=x_old)$y
  
  return(dens)
}

InternalCdf <- function(q, param, mynames, time.names,  
                        time.var, up, low, N_grid, method, myfun, N_Fourier){
  x_old <- q
  x_new <-unique(sort(c(low, x_old, up)))
  x_new <- unique(sort(c(seq(min(x_new)-0.1,max(x_new)+0.1, length.out =N_grid+1),x_old)))
  dens <- InternalDensity(x_new, param, mynames, time.names, time.var, up, low, N_grid, method , myfun, N_Fourier)
  cdf <- c(0,cumsum(as.numeric(dens)[-1]*diff(x_new)))
  res <- na.approx(zoo(cdf, order.by=sort(x_new)), xout = q)
  return(res)
}

InternalRnd <- function(n, param, mynames, time.names,  
                        time.var, up, low, N_grid, method, myfun, N_Fourier){
  x_new <- seq(low-0.1,up+0.1, length.out =N_grid+1)
  cdf<-as.numeric(InternalCdf(q=x_new, param, mynames, time.names,time.var, up, low, N_grid, method, myfun, N_Fourier))
  cdf0 <- cdf[cdf>0 & cdf<1]
  x_new0 <- x_new[cdf>0 & cdf<1]
  rndUn0 <- runif(n, min = min(cdf0), max(cdf0))
  res <- approx(y=x_new0, x = cdf0, xout=rndUn0, ties = mean, rule =2)
  if(length(res$y)==n){
    return(as.numeric(res$y))
  }else{
    m <- n-length(res$y)
    rndUn1 <- runif(m, min = min(cdf0), max(cdf0))
    res1 <- approx(y=x_new0, x = cdf0, xout=rndUn1)
    res_new <- c(res$y,res1$y)
    if(length(res_new)==n){
      return(res_new)
    }else{
      m<-n-length(res_new)
      return(c(res_new, sample(res_new, m)))
    }
  }
}

InternalQnt <- function(p, param, mynames, time.names,  
                        time.var, up, low, N_grid, method, myfun, N_Fourier){
  x_new <- seq(low-0.1,up+0.1, length.out =N_grid+1)
  cdf<-as.numeric(InternalCdf(q=x_new, param, mynames, time.names,time.var, up, low, N_grid, method, myfun, N_Fourier))
  cdf0 <- cdf[cdf>0 & cdf<1]
  x_new0 <- x_new[cdf>0 & cdf<1]
  res <- approx(y=x_new0, x = cdf0, xout=p)
  return(res$y)
}

FromCF2yuima_law <- function(myfun, time.names = "t", var_char = "u",
                             up = 45, low = -45, N_grid = 50001, N_Fourier=2^10){
  method <- "FFT"
  if(!myfun %in% names(globalenv())){
    yuima.stop("the characteristi function is not defined in the global enviromnent: check arg myfun")
  }
  dumEval <- parse(text = paste0("dumFun <- ", myfun)) 
  lawpar<-extrapParam(eval(dumEval),time.names, var_char)
  # if(!all(names(true.par) %in% lawpar)){
  #   stop("error massage")
  # }
  true.par<- lawpar
  mynames <- lawpar
  nametime <- time.names
  dumm1 <- paste0("(", paste0(c("u",lawpar,time.names), collapse=", "),")")
  
  mystring <- paste0(paste(paste0(lawpar, collapse=", "), time.names, sep =", " ),"){")
  mystring <- paste(mystring, paste0(" up <- ",up) ,sep="\n")
  mystring <- paste(mystring, paste0(" low <- ", low) ,sep="\n")
  mystring <- paste(mystring, paste0(" N_grid <- ", N_grid) ,sep="\n")
  mystring <- paste(mystring, paste0(" param <- c(", paste0(lawpar, collapse=" ,"), ")") ,sep="\n")
  mystring <- paste(mystring, paste0(" mynames <- c('", paste0(mynames, collapse="' ,'"), "')") ,sep="\n")
  mystring <- paste(mystring, paste0(" time.var <- c(", paste0(nametime, collapse=" ,"), ")") ,sep="\n")
  mystring <- paste(mystring, paste0(" time.names <- c('", paste0(nametime, collapse="' ,"), "')") ,sep="\n")
  mystring <- paste(mystring, paste0(" method <- c('", paste0(method, collapse="' ,"), "')") ,sep="\n")
  mystring <- paste(mystring, paste0(" myfun <- '", paste0(paste0(myfun, collapse="' ,"),dumm1), "'") ,sep="\n")
  mystring <- paste(mystring, paste0(" N_Fourier <- ", N_Fourier) ,sep="\n")
  
  mystring_dens_in <- "dmyLaw <- function(x, "
  mystring_dens_fin <- paste(mystring, paste(" res <- InternalDensity(x", "param", "mynames", "time.names",  
                                             "time.var", "up", "low", "N_grid",  
                                             "method ",
                                             "myfun",
                                             "N_Fourier",
                                             sep = ", " ) ,sep="\n")
  mystring_dens <- paste0(mystring_dens_in,mystring_dens_fin)
  mystring_dens <- paste(mystring_dens,")", sep="") 
  mystring_dens <- paste(mystring_dens, " return(res)","}" ,sep="\n")
  
  mystring_cdf_in <- "pmyLaw <- function(q, "
  mystring_cdf_fin <- paste(mystring, paste(" res <- InternalCdf(q", "param", "mynames", "time.names",  
                                            "time.var", "up", "low", "N_grid",  
                                            "method ",
                                            "myfun",
                                            "N_Fourier",
                                            sep = ", " ) ,sep="\n")
  #cat(mystring)
  mystring_cdf <- paste0(mystring_cdf_in,mystring_cdf_fin)
  mystring_cdf <- paste(mystring_cdf,")", sep="") 
  mystring_cdf <- paste(mystring_cdf, " return(res)","}" ,sep="\n")
  
  mystring_rnd_in <- "rmyLaw <- function(n, "
  mystring_rnd_fin <- paste(mystring, paste(" res <- InternalRnd(n", "param", "mynames", "time.names",  
                                            "time.var", "up", "low", "N_grid",  
                                            "method ",
                                            "myfun",
                                            "N_Fourier",
                                            sep = ", " ) ,sep="\n")
  #cat(mystring)
  mystring_rnd <- paste0(mystring_rnd_in,mystring_rnd_fin)
  mystring_rnd <- paste(mystring_rnd,")", sep="") 
  mystring_rnd <- paste(mystring_rnd, " return(res)","}" ,sep="\n")
  
  mystring_qnt_in <- "qmyLaw <- function(p, "
  mystring_qnt_fin <- paste(mystring, paste(" res <- InternalQnt(p", "param", "mynames", "time.names",  
                                            "time.var", "up", "low", "N_grid",  
                                            "method ",
                                            "myfun",
                                            "N_Fourier",
                                            sep = ", " ) ,sep="\n")
  #cat(mystring)
  mystring_qnt <- paste0(mystring_qnt_in,mystring_qnt_fin)
  mystring_qnt <- paste(mystring_qnt,")", sep="") 
  mystring_qnt <- paste(mystring_qnt, " return(res)","}" ,sep="\n")
  dmyLaw <- NULL
  pmyLaw <- NULL
  rmyLaw <- NULL
  qmyLaw <- NULL
  eval(parse(text=mystring_dens))
  eval(parse(text=mystring_cdf))
  eval(parse(text=mystring_rnd))
  eval(parse(text=mystring_qnt))
  res<-setLaw(density = dmyLaw)
#  res@characteristic<-myChar
  res@cdf <- pmyLaw
  res@rng <- rmyLaw
  res@quantile <- qmyLaw
  return(res)
}




