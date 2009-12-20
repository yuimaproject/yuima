## We have splitted the simulate function into blocks to allow for future 
## methods to be added. S.M.I. & A.B.
## Interface to simulate() changed to match the S3 generic funciton in the 
## package stats
## added an environment to let R find the proper values defined in the main
## body of the function, which in turn calls different simulation methods
## All new simulation methods should look into the yuimaEnv for local variables
## when they need to "eval" R expressions

##:: function simulate
##:: solves SDE and returns result
setGeneric("simulate",
			function(object, nsim, seed, xinit, true.parameter, space.discretized=FALSE, increment.W=NULL, increment.L=NULL, methodfGn="Cholesky")
			standardGeneric("simulate")
           )

setMethod("simulate", "yuima",
			function(object, nsim=1, seed=NULL, xinit, true.parameter, space.discretized=FALSE, increment.W=NULL, increment.L=NULL,methodfGn="Cholesky"){


##:: errors checks

##:1: error on yuima model
  yuima <- object
				
  if(missing(yuima)){
    cat("\nyuima object is missing.\n")
    return(NULL)
  }
  
				sdeModel <- yuima@model
				Terminal <- yuima@sampling@Terminal[1]
				n <- yuima@sampling@n[1]
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
				


				
  par.len <- length(sdeModel@parameter@all)
				
  if(missing(true.parameter) & par.len>0){
	true.parameter <- vector(par.len, mode="list")
	for(i in 1:par.len)
     true.parameter[[i]] <- 0
	names(true.parameter) <-   sdeModel@parameter@all
  }
    
  yuimaEnv <- new.env()
				
  if(par.len>0){
	for(i in 1:par.len){
	 pars <- sdeModel@parameter@all[i]
	 assign(sdeModel@parameter@all[i], true.parameter[[i]], yuimaEnv)
	}
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
    }else if(dim(increment.W)[2] != n){
      cat("\nLength of increment's column must be same as yuima@sampling@n[1].\n")
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
    }else if(dim(increment.L)[2] != n){
      cat("\nLength of increment's column must be same as yuima@sampling@n[1].\n")
      return(NULL)
    }
  }
  
    

				
  
  if(space.discretized){   	  
	  ##:: using Space-discretized Euler-Maruyama method	  
	  yuima@data <- space.discretized(xinit, yuima, yuimaEnv)
	  return(yuima)  
  }
	
	
##:: using Euler-Maruyama method
  delta <- Terminal/n 
				

				
				if(missing(increment.W)){
					if( sdeModel@hurst!=0.5 ){ 
						grid<-sampling2grid(yuima@sampling)	
						
						if(methodfGn=="Cholesky"){
						dW<-CholeskyfGn(grid, sdeModel@hurst,r.size)
						}else{
						cat("\nNot done presently\n")
						return(NULL)	
						}
						
					} else {
						delta<-Terminal/n
						dW <- rnorm(n * r.size, 0, sqrt(delta))
						dW <- matrix(dW, ncol=n, nrow=r.size,byrow=TRUE)  
					}
				} else {
					dW <- increment.W
				}
				
  yuima@data <- euler(xinit, yuima, dW, yuimaEnv)
  return(yuima)
})
