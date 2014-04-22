yuima.stop <- function(x) 
stop(sprintf("\nYUIMA: %s\n", x))

yuima.warn <- function(x) 
   warning(sprintf("\nYUIMA: %s\n", x))

# 22/11/2013
# We introduce a new utility yuima.simplify that allows us to simplify 
# the expressions in the drift, diffusion and jump terms.

# yuima.Simplify modified from the original code  Simplify.R 
# by Andrew Clausen <clausen@econ.upenn.edu> in 2007. 
# http://economics.sas.upenn.edu/~clausen/computing/Simplify.R

# This isn't a serious attempt at simplification code.  It just does some
# obvious things like 0 + x => x.  It was written to support Deriv.R.

yuima.Simplify <- function(expr, yuima.env){
    
    ###
    
    
    Simplify_ <- function(expr)
    {
        if (is.symbol(expr)) {
            expr
        } else if (is.language(expr) && is.symbol(expr[[1]])) {
            # is there a rule in the table?
            sym.name <- as.character(expr[[1]])
            if (class(try(Simplify.rule <-
            get(sym.name, envir=yuima.env,
            inherits=FALSE), silent=TRUE))
            != "try-error")
            return(Simplify.rule(expr))
        }
        expr
    }
    
    Simplify.function <- function(f, x=names(formals(f)), env=parent.frame())
    {
        stopifnot(is.function(f))
        as.function(c(as.list(formals(f)),
        Simplify_(body(f))),
        envir=env)
    }
    
    `Simplify.+` <- function(expr)
    {
        if (length(expr) == 2)
        {
            if (is.numeric(expr[[2]]))
            return(+expr[[2]])
            return(expr)
        }
        
        a <- Simplify_(expr[[2]])
        b <- Simplify_(expr[[3]])
        
        if (is.numeric(a) && all(a == 0)) {
            b
        } else if (is.numeric(b) && all(b == 0)) {
            a
        } else if (is.numeric(a) && is.numeric(b)) {
            a + b
        } else {
            expr[[2]] <- a
            expr[[3]] <- b
            expr
        }
    }
    
    `Simplify.-` <- function(expr)
    {
        if (length(expr) == 2)
        {
            if (is.numeric(expr[[2]]))
            return(-expr[[2]])
            return(expr)
        }
        
        a <- Simplify_(expr[[2]])
        b <- Simplify_(expr[[3]])
        
        if (is.numeric(a) && all(a == 0)) {
            -b
        } else if (is.numeric(b) && all(b == 0)) {
            a
        } else if (is.numeric(a) && is.numeric(b)) {
            a - b
        } else {
            expr[[2]] <- a
            expr[[3]] <- b
            expr
        }
    }
    
    `Simplify.(` <- function(expr)
    expr[[2]]
    
    `Simplify.*` <- function(expr)
    {
        a <- Simplify_(expr[[2]])
        b <- Simplify_(expr[[3]])
        
        if (is.numeric(a) && all(a == 0)) {
            0
        } else if (is.numeric(b) && all(b == 0)) {
            0
        } else if (is.numeric(a) && all(a == 1)) {
            b
        } else if (is.numeric(b) && all(b == 1)) {
            a
        } else if (is.numeric(a) && is.numeric(b)) {
            a * b
        } else {
            expr[[2]] <- a
            expr[[3]] <- b
            expr
        }
    }
    
    `Simplify.^` <- function(expr)
    {
        a <- Simplify_(expr[[2]])
        b <- Simplify_(expr[[3]])
        
        if (is.numeric(a) && all(a == 0)) {
            0
        } else if (is.numeric(b) && all(b == 0)) {
            1
        } else if (is.numeric(a) && all(a == 1)) {
            1
        } else if (is.numeric(b) && all(b == 1)) {
            a
        } else if (is.numeric(a) && is.numeric(b)) {
            a ^ b
        } else {
            expr[[2]] <- a
            expr[[3]] <- b
            expr
        }
    }
    
    `Simplify.c` <- function(expr)
    {
        args <- expr[-1]
        args.simplified <- lapply(args, Simplify_)
        if (all(lapply(args.simplified, is.numeric))) {
            as.numeric(args.simplified)
        } else {
            for (i in 1:length(args))
            expr[[i + 1]] <- args.simplified[[i]]
            expr
        }
    }
    
    
    assign("+", `Simplify.+`, envir=yuima.env)
    assign("-", `Simplify.-`, envir=yuima.env)
    assign("*", `Simplify.*`, envir=yuima.env)
    assign("(", `Simplify.(`, envir=yuima.env)
    assign("c", `Simplify.c`, envir=yuima.env)
    assign("^", `Simplify.^`, envir=yuima.env)
    

    ###
    
    
    
    
    
    
    
    
  as.expression(Simplify_(expr[[1]]))

}


## Constructor and Initializer of class 'yuima'

# we convert objects to "zoo" internally
# we should change it later to more flexible classes

setMethod("initialize", "yuima",
          function(.Object, data=NULL, model=NULL, sampling=NULL, characteristic=NULL, functional=NULL){  
            eqn <- NULL

            if(!is.null(data)){
              .Object@data <- data
              eqn <- dim(data)
            }
            
            if(!is.null(model)){
              if(!is.null(eqn)){
                if(eqn!=model@equation.number){
                  yuima.warn("Model's equation number missmatch.")
                  return(NULL)
                }
              }else{
                eqn <- model@equation.number
              }
              .Object@model <- model
            }
            
            if(!is.null(sampling)){
              if(!is.null(eqn)){
                if(eqn!=length(sampling@Terminal)){
                  if(length(sampling@Terminal)==1){
                    sampling@Terminal <- rep(sampling@Terminal, eqn)
                    sampling@n <- rep(sampling@n, eqn)
                  }else{
                    yuima.warn("Sampling's equation number missmatch.")
                    return(NULL)
                  }
                }
              }else{
                eqn <- length(sampling@Terminal)
              }
              .Object@sampling <- sampling
            }
            
            if(!is.null(characteristic)){
              if(!is.null(eqn)){
                if(eqn!=characteristic@equation.number){
                  yuima.warn("Characteristic's equation number missmatch.")
                  return(NULL)                  
                }
              }
              .Object@characteristic <- characteristic
            }else if(!is.null(eqn)){
              characteristic <- new("yuima.characteristic", equation.number=eqn, time.scale=1)
              .Object@characteristic <- characteristic
            }
            
			if(!is.null(functional)) .Object@functional <- functional
			
            return(.Object)
          })

# setter
setYuima <-
  function(data=NULL, model=NULL, sampling=NULL, characteristic=NULL, functional=NULL){
    return(new("yuima", data=data, model=model, sampling=sampling, characteristic=characteristic,functional=functional))
  }




setMethod("show", "yuima.functional",
function(object){
    str(object)
} )

setMethod("show", "yuima.sampling",
function(object){
    str(object)
} )


setMethod("show", "yuima.data",
function(object){
    show(setYuima(data=object))
    
} )


setMethod("show", "yuima.model",
function(object){
    show(setYuima(model=object))
    
} )

setMethod("show", "yuima",
function(object){
    
    mod <- object@model
    has.drift <- FALSE
    has.diff <- FALSE
    has.fbm <- FALSE
    has.levy <- FALSE
    is.wienerdiff <- FALSE
    is.fracdiff <- FALSE
    is.jumpdiff <- FALSE
    
    if(length(mod@drift)>0) has.drift <- TRUE
    if(length(mod@diffusion)>0) has.diff <- TRUE
    if(length(mod@jump.coeff)>0) has.levy <- TRUE
    if(!is.null(mod@hurst)) has.fbm <- TRUE
    
    if( has.drift | has.diff ) is.wienerdiff <- TRUE
    if( has.fbm  ) is.fracdiff <- TRUE
    if( has.levy ) is.jumpdiff <- TRUE
    
    if( is.wienerdiff | is.fracdiff | is.jumpdiff  ){
        if( is.wienerdiff )
        cat("\nDiffusion process")
        if( is.fracdiff & mod@hurst!=0.5)
        cat(sprintf(" with Hurst index:%.2f", mod@hurst))
        if( is.jumpdiff ){
            if( is.wienerdiff ){
                cat(" with Levy jumps")
            } else {
                cat("Levy jump process")
            }
        }
        
        cat(sprintf("\nNumber of equations: %d", mod@equation.number))
        cat(sprintf("\nNumber of Wiener noises: %d", length(mod@diffusion)))
    }
    
    if(length(object@data@original.data)>0){
        n.series <- 1
        if(!is.null(dim(object@data@original.data))){
            n.series <- dim(object@data@original.data)[2]
            n.length <- dim(object@data@original.data)[1]
        } else {
            n.length <- length(object@data@original.data)
        }
        
        cat(sprintf("\n\nNumber of original time series: %d\nlength = %d, time range [%s ; %s]", n.series, n.length, min(time(object@data@original.data)), max(time(object@data@original.data))))
    }
    if(length(object@data@zoo.data)>0){
        n.series <- length(object@data@zoo.data)
        n.length <- unlist(lapply(object@data@zoo.data, length))
        t.min <- unlist(lapply(object@data@zoo.data, function(u) as.character(round(time(u)[which.min(time(u))],3))))
        t.max <- unlist(lapply(object@data@zoo.data, function(u) as.character(round(time(u)[which.max(time(u))],3))))
        
        delta <- NULL
        for(i in 1:n.series){
            tmp <- deltat(object@data@zoo.data[[i]])
            if(is.null(tmp)){
                delta <- c(delta, NA)
            } else {
                delta <- c(delta, tmp)
            }
        }
        
        
        cat(sprintf("\n\nNumber of zoo time series: %d\n", n.series))
        tmp <- data.frame(length=n.length, time.min = t.min, time.max =t.max, delta=delta)
        rownames(tmp) <- sprintf("Series %d",1:n.series)
        print(tmp)
    }
    
})

