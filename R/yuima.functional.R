# constructor and initializer of class 'yuima.functional'

setMethod("initialize", "yuima.functional",
           function(.Object, F, f, xinit, e){
             .Object@F <- F
             .Object@f <- f
             .Object@xinit <- xinit
             #.Object@Terminal <- Terminal
             #.Object@division <- division
             .Object@e <- e
             return(.Object)
           })


# setter

setGeneric("setFunctional",
           function(model, F, f, xinit, e)
           standardGeneric("setFunctional")
           )
setMethod("setFunctional", "yuima",
  function(model, F, f, xinit, e){
    model@functional <- setFunctional(model@model,F,f,xinit,e)@functional
    return(model)
  })

setMethod("setFunctional", "yuima.model",
  function(model, F, f, xinit, e){
    # error check
    if( missing(model)){
      cat("\nyuima.model is missing.\n")
      return(NULL)
    }
    if( missing(xinit)){
      cat("\nInitial value of state variable is missing\n")
      return(NULL)
    }
    if( missing(f) || missing(F) || missing(e)){
      cat("\nFunctional specification is incomplete.\n")
      return(NULL)
    }

    r.size <- model@noise.number
    d.size <- model@equation.number
    k.size <- length(F)

    if( length(f) != (r.size + 1)){
      cat("\nFunctional needs r+1 f_alphas\n")
      return(NULL)
    }
    if( length(f[[1]]) != k.size){
      cat("\nMissmatch in dimension of functional\n")
      return(NULL)
    }
    if( length(xinit) != d.size){
      cat("\nMissmatch in dimension of functional and state variables\n")
      return(NULL)
    }
    
	# instanciate
    return(setYuima(model = model,functional = new("yuima.functional", F=F, f=f, xinit=xinit, e=e )))
  })

# getter of each variables
setGeneric("getF",
           function(x)
           standardGeneric("getF")
           )
setMethod("getF", "yuima.functional",
          function(x){
            return(x@F)
          })
setGeneric("getf",
           function(x)
           standardGeneric("getf")
           )
setMethod("getf", "yuima.functional",
          function(x){
            return(x@f)
          })
setGeneric("getxinit",
           function(x)
           standardGeneric("getxinit")
           )
setMethod("getxinit", "yuima.functional",
          function(x){
            return(x@xinit)
          })
setGeneric("gete",
           function(x)
           standardGeneric("gete")
           )
setMethod("gete", "yuima.functional",
          function(x){
            return(x@e)
          })
