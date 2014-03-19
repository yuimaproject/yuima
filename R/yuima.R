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
