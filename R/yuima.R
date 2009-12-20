## Constructor and Initializer of class 'yuima'

# we convert objects to "zoo" internally
# we should change it later to more flexible classes

setMethod("initialize", "yuima",
          function(.Object, data=NULL, model=NULL, sampling=NULL, characteristic=NULL){  
            eqn <- NULL

            if(!is.null(data)){
              .Object@data <- data
              eqn <- dim(data)
            }
            
            if(!is.null(model)){
              if(!is.null(eqn)){
                if(eqn!=model@equation.number){
                  cat("\nModel's equation number missmatch.\n")
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
                    cat("\nSampling's equation number missmatch.\n")
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
                  cat("\nCharacteristic's equation number missmatch.\n")
                  return(NULL)                  
                }
              }
              .Object@characteristic <- characteristic
            }else if(!is.null(eqn)){
              characteristic <- new("yuima.characteristic", equation.number=eqn, time.scale=1)
              .Object@characteristic <- characteristic
            }
            
            return(.Object)
          })

# setter
setYuima <-
  function(data=NULL, model=NULL, sampling=NULL, characteristic=NULL){
    return(new("yuima", data=data, model=model, sampling=sampling, characteristic=characteristic))
  }
