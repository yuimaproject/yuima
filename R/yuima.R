yuima.warn <- function(x) 
   cat(sprintf("\nYUIMA: %s\n", x))

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
