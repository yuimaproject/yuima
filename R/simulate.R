##:: function simulate
##:: solves SDE and returns result
setGeneric("simulate",
			function(yuima, xinit, true.parameter, space.discretized=FALSE, increment.W=NULL, increment.L=NULL)
			standardGeneric("simulate")
           )

setMethod("simulate", "yuima",
			function(yuima, xinit, true.parameter, space.discretized=FALSE, increment.W=NULL, increment.L=NULL){


##:: errors checks

##:1: error on yuima model
  if(missing(yuima)){
    cat("\nyuima object is missing.\n")
    return(NULL)
  }
  
				sdeModel <- yuima@model
				Terminal <- yuima@sampling@Terminal[1]
				division <- yuima@sampling@division[1]
				r.size <- sdeModel@noise.number
				d.size <- sdeModel@equation.number				
				
##:2: error on xinit 
  if(missing(xinit)){
	xinit <- sdeModel@xinit
  }else if(length(xinit) != d.size){
					if(length(xinit)==1){
						xinit <- rep(xinit, d.size)
					}else{
						cat("\nDimension of xinit variables missmuch.\n")
						return(NULL)
					}
				}
				

  
  ##   str(sdeModel)
  ##   print(Terminal)
  ##   print(yuima@sampling@division)
  ##   readline()
  
  if(missing(true.parameter)){
    true.parameter <- numeric(length(sdeModel@parameter@all))
  }
  
 
  
  
  
  if(space.discretized){
    if(r.size>1){
      warning("Space-discretized EM cannot be used for multi-dimentional models. Use standard method.")
      space.discretized <- FALSE
    }
    if(length(sdeModel@jump.coeff)){
      warning("Space-discretized EM is for only Wiener Proc. Use standard method.")
      space.discretized <- FALSE
    }
  }

  ##:: Error check for increment specified version.
  if(!missing(increment.W)){
    if(space.discretized == TRUE){
      cat("\nParameter increment must be invalid if space.discretized=TRUE.\n")
      return(NULL)
    }else if(dim(increment.W)[1] != r.size){
      cat("\nLength of increment's row must be same as yuima@model@noise.number.\n")
      return(NULL)
    }else if(dim(increment.W)[2] != division){
      cat("\nLength of increment's column must be same as yuima@sampling@division[1].\n")
      return(NULL)
    }
  }

    ##:: Error check for increment specified version.
  if(!missing(increment.L)){
    if(space.discretized == TRUE){
      cat("\nParameter increment must be invalid if space.discretized=TRUE.\n")
      return(NULL)
    }else if(dim(increment.L)[1] != r.size){
      cat("\nLength of increment's row must be same as yuima@model@noise.number.\n")
      return(NULL)
    }else if(dim(increment.L)[2] != division){
      cat("\nLength of increment's column must be same as yuima@sampling@division[1].\n")
      return(NULL)
    }
  }
  
    
  par.len <- length(sdeModel@parameter@all)
  if(par.len>0){
    for(i in 1:par.len){
      pars <- sdeModel@parameter@all[i]
      assign(pars, true.parameter[i])
    }
  }

				
  
  if(space.discretized){   	  
	  ##:: using Space-discretized Euler-Maruyama method	  
	  yuima@data <- space.discretized(xinit, yuima)
	  return(yuima)  
  }
	
	
    
  ##:: using Euler-Maruyama method
  
  ##:: Diffusion terms
  if( missing(increment.W)){
  	
    dW <- rnorm(division * r.size, 0, sqrt(delta))
    dW <- t(matrix(dW, nrow=division, ncol=r.size))
  }else{
    dW <- increment.W
  }
    
  yuima@data <- euler(xinit,yuima,dW)
  return(yuima)
})
