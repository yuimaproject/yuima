## We have splitted the simulate function into blocks to allow for future 
## methods to be added. S.M.I. & A.B.
## Interface to simulate() changed to match the S3 generic function in the 
## package stats
## added an environment to let R find the proper values defined in the main
## body of the function, which in turn calls different simulation methods
## All new simulation methods should look into the yuimaEnv for local variables
## when they need to "eval" R expressions

##:: function simulate
##:: solves SDE and returns result
subsampling <- function(x,y) return(x)

setGeneric("simulate",
           function(object, nsim=1, seed=NULL, xinit, true.parameter, space.discretized=FALSE, 
                    increment.W=NULL, increment.L=NULL, hurst, methodfGn="WoodChan", 
                    sampling=sampling, subsampling=subsampling, ...
                    #		Initial = 0, Terminal = 1, n = 100, delta, 
                    #		grid = as.numeric(NULL), random = FALSE, sdelta=as.numeric(NULL), 
                    #		sgrid=as.numeric(NULL), interpolation="none"
           )
             standardGeneric("simulate")
)


setMethod("simulate", "yuima.model",
          function(object, nsim=1, seed=NULL, xinit, true.parameter, 
                   space.discretized=FALSE, increment.W=NULL, increment.L=NULL,
                   hurst, methodfGn="WoodChan",
                   sampling, subsampling,
                   #Initial = 0, Terminal = 1, n = 100, delta, 
                   #	grid, random = FALSE, sdelta=as.numeric(NULL), 
                   #	sgrid=as.numeric(NULL), interpolation="none"
                   ...){
            
            tmpsamp <- NULL
            if(missing(sampling)){
              tmpsamp <- setSampling(...)
              #		 tmpsamp <- setSampling(Initial = Initial, Terminal = Terminal, n = n, 
              #				delta = delta, grid = grid, random = random, sdelta=sdelta, 
              #				sgrid=sgrid, interpolation=interpolation)
            } else {
              tmpsamp <- sampling
            }
            tmpyuima <- setYuima(model=object, sampling=tmpsamp) 	
            
            out <- simulate(tmpyuima, nsim=nsim, seed=seed, xinit=xinit, 
                            true.parameter=true.parameter, 
                            space.discretized=space.discretized, 
                            increment.W=increment.W, increment.L=increment.L,
                            hurst=hurst,methodfGn=methodfGn, subsampling=subsampling)
            return(out)	
          })


# We rewrite the code of simulate method. We build a new internal function aux.simulate containing 
# the old code. This function starts directly if the model is an object of yuima.model-class 
# or yuima.carma-class. If the model is an object of class yuima.cogarch, simulate method runs  
# the internal function aux.simulateCogarch that generates first a path of the underlying levy and then
# uses this path to construct the trajectories of the Cogarch model by calling the aux.simulate function.

setMethod("simulate", "yuima",
          function(object, nsim=1, seed=NULL, xinit, true.parameter, 
                   space.discretized=FALSE, increment.W=NULL, increment.L=NULL,
                   hurst,methodfGn="WoodChan",
                   sampling, subsampling,
                   Initial = 0, Terminal = 1, n = 100, delta, 
                   grid = as.numeric(NULL), random = FALSE, sdelta=as.numeric(NULL), 
                   sgrid=as.numeric(NULL), interpolation="none"){
            
            if(is(object@model,"yuima.cogarch")){
              res<-aux.simulateCogarch(object, nsim, seed, xinit, true.parameter, 
                           space.discretized, increment.W, increment.L,
                           hurst,methodfGn,
                           sampling, subsampling,
                           Initial, Terminal, n, delta, 
                           grid, random, sdelta, 
                           sgrid, interpolation)
            }else{
              res<-aux.simulate(object, nsim, seed, xinit, true.parameter, 
                                   space.discretized, increment.W, increment.L,
                                   hurst,methodfGn,
                                   sampling, subsampling,
                                   Initial, Terminal, n, delta, 
                                   grid, random, sdelta, 
                                   sgrid, interpolation)
            }
            return(res)
            
#             ##:: errors checks
#             
#             ##:1: error on yuima model
#             yuima <- object
#             
#             if(missing(yuima)){
#               yuima.warn("yuima object is missing.")
#               return(NULL)
#             }
#             
#             tmphurst<-yuima@model@hurst      
#             
#             if(!missing(hurst)){
#               yuima@model@hurst=hurst
#             }      
#             
#             if (is.na(yuima@model@hurst)){
#               yuima.warn("Specify the hurst parameter.")
#               return(NULL)
#             }      
#             
#             tmpsamp <- NULL
#             if(is.null(yuima@sampling)){
#               if(missing(sampling)){
#                 tmpsamp <- setSampling(Initial = Initial, 
#                                        Terminal = Terminal, n = n, delta = delta, 
#                                        grid = grid, random = random, sdelta=sdelta, 
#                                        sgrid=sgrid, interpolation=interpolation)
#               } else {
#                 tmpsamp <- sampling
#               }
#             } else {
#               tmpsamp <- yuima@sampling
#             }
#             
#             yuima@sampling <- tmpsamp		
#             
#             sdeModel <- yuima@model
#             Terminal <- yuima@sampling@Terminal[1]
#             n <- yuima@sampling@n[1]
#             r.size <- sdeModel@noise.number
#             d.size <- sdeModel@equation.number				
#             
#             ##:2: error on xinit 
#             if(missing(xinit)){
#               xinit <- sdeModel@xinit
#             } else {
#               if(length(xinit) != d.size){
#                if(length(xinit)==1){
#                  xinit <- rep(xinit, d.size)
#                } else {
#                 yuima.warn("Dimension of xinit variables missmatch.")
#                 return(NULL)
#                }
#               }
#             }
#             
#             xinit <- as.expression(xinit)  # force xinit to be an expression
#             
# 
#             par.len <- length(sdeModel@parameter@all)
#             
#             if(missing(true.parameter) & par.len>0){
#               true.parameter <- vector(par.len, mode="list")
#               for(i in 1:par.len)
#                 true.parameter[[i]] <- 0
#               names(true.parameter) <-   sdeModel@parameter@all
#             }
#             
#             yuimaEnv <- new.env()
#             
#             if(par.len>0){
#               for(i in 1:par.len){
#                 pars <- sdeModel@parameter@all[i]
#                 
#                 for(j in 1:length(true.parameter)){
#                   if( is.na(match(pars, names(true.parameter)[j]))!=TRUE){
#                     assign(sdeModel@parameter@all[i], true.parameter[[j]],yuimaEnv)
#                   }
#                 }
#                 #assign(sdeModel@parameter@all[i], true.parameter[[i]], yuimaEnv)
#               }
#             }
#             
#             
#             if(space.discretized){
#               if(r.size>1){
#                 warning("Space-discretized EM cannot be used for multi-dimentional models. Use standard method.")
#                 space.discretized <- FALSE
#               }
#               if(length(sdeModel@jump.coeff)){
#                 warning("Space-discretized EM is for only Wiener Proc. Use standard method.")
#                 space.discretized <- FALSE
#               }
#             }
#             
#             ##:: Error check for increment specified version.
#             if(!missing(increment.W) & !is.null(increment.W)){
#               if(space.discretized == TRUE){
#                 yuima.warn("Parameter increment must be invalid if space.discretized=TRUE.")
#                 return(NULL)
#               }else if(dim(increment.W)[1] != r.size){
#                 yuima.warn("Length of increment's row must be same as yuima@model@noise.number.")
#                 return(NULL)
#               }else if(dim(increment.W)[2] != n){
#                 yuima.warn("Length of increment's column must be same as sampling@n[1].")
#                 return(NULL)
#               }
#             }
#             
#             ##:: Error check for increment specified version.
#             if(!missing(increment.L) & !is.null(increment.L)){
#               if(space.discretized == TRUE){
#                 yuima.warn("Parameter increment must be invalid if space.discretized=TRUE.")
#                 return(NULL)
#               }else if(dim(increment.L)[1] != length(yuima@model@jump.coeff[[1]]) ){ #r.size){
#                 yuima.warn("Length of increment's row must be same as yuima@model@noise.number.")
#                 return(NULL)
#               }else if(dim(increment.L)[2] != n){
#                 yuima.warn("Length of increment's column must be same as sampling@n[1].")
#                 return(NULL)
#               }
#             }
#             
#             
#             yuimaEnv$dL <- increment.L
#             
#             
#             if(space.discretized){   	  
#               ##:: using Space-discretized Euler-Maruyama method	  
#               yuima@data <- space.discretized(xinit, yuima, yuimaEnv)
#               
#               yuima@model@hurst<-tmphurst
#               return(yuima)  
#             }
#             
#             
#             ##:: using Euler-Maruyama method
#             delta <- Terminal/n 
#             
#             if(missing(increment.W) | is.null(increment.W)){
#               
#               if( sdeModel@hurst!=0.5 ){ 
#                 
#                 grid<-sampling2grid(yuima@sampling)	
#                 isregular<-yuima@sampling@regular 
#                 
#                 if((!isregular) || (methodfGn=="Cholesky")){
#                   dW<-CholeskyfGn(grid, sdeModel@hurst,r.size)
#                   yuima.warn("Cholesky method for simulating fGn has been used.")
#                 } else {
#                   dW<-WoodChanfGn(grid, sdeModel@hurst,r.size)
#                 }
#                 
#               } else {
#                 
#                 delta<-Terminal/n
#                 if(!is.Poisson(sdeModel)){ # if pure CP no need to setup dW
#                  dW <- rnorm(n * r.size, 0, sqrt(delta))
#                  dW <- matrix(dW, ncol=n, nrow=r.size,byrow=TRUE)
#                 } else {
#                     dW <- matrix(0,ncol=n,nrow=1)  # maybe to be fixed
#                 }
#               }
#               
#             } else {
#               dW <- increment.W
#             }
#             
#             if(is.Poisson(sdeModel)){
#                 yuima@data <- simCP(xinit, yuima, yuimaEnv)
#             } else {
#                 yuima@data <- euler(xinit, yuima, dW, yuimaEnv)
#             }
#             
#             for(i in 1:length(yuima@data@zoo.data)) 
#               index(yuima@data@zoo.data[[i]]) <- yuima@sampling@grid[[1]]  ## to be fixed
#             yuima@model@xinit <- xinit
#             yuima@model@hurst <-tmphurst
#             
#             if(missing(subsampling))
#               return(yuima)
#             subsampling(yuima, subsampling)
#             
           }
          )

aux.simulate<-function(object, nsim, seed, xinit, true.parameter, 
         space.discretized, increment.W, increment.L,
         hurst,methodfGn,
         sampling, subsampling,
         Initial, Terminal, n, delta, 
         grid, random, sdelta, 
         sgrid, interpolation){
  
  ##:: errors checks
  
  ##:1: error on yuima model
  yuima <- object
  
  if(missing(yuima)){
    yuima.warn("yuima object is missing.")
    return(NULL)
  }
  
  tmphurst<-yuima@model@hurst      
  
  if(!missing(hurst)){
    yuima@model@hurst=hurst
  }      
  
  if (is.na(yuima@model@hurst)){
    yuima.warn("Specify the hurst parameter.")
    return(NULL)
  }      
  
  tmpsamp <- NULL
  if(is.null(yuima@sampling)){
    if(missing(sampling)){
      tmpsamp <- setSampling(Initial = Initial, 
                             Terminal = Terminal, n = n, delta = delta, 
                             grid = grid, random = random, sdelta=sdelta, 
                             sgrid=sgrid, interpolation=interpolation)
    } else {
      tmpsamp <- sampling
    }
  } else {
    tmpsamp <- yuima@sampling
  }
  
  yuima@sampling <- tmpsamp  	
  
  sdeModel <- yuima@model
  Terminal <- yuima@sampling@Terminal[1]
  n <- yuima@sampling@n[1]
  r.size <- sdeModel@noise.number
  d.size <- sdeModel@equation.number				
  
  ##:2: error on xinit 
  if(missing(xinit)){
    xinit <- sdeModel@xinit
  } else {
    if(length(xinit) != d.size){
      if(length(xinit)==1){
        xinit <- rep(xinit, d.size)
      } else {
        yuima.warn("Dimension of xinit variables missmatch.")
        return(NULL)
      }
    }
  }
  
  xinit <- as.expression(xinit)  # force xinit to be an expression
  
  
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
      
      for(j in 1:length(true.parameter)){
        if( is.na(match(pars, names(true.parameter)[j]))!=TRUE){
          assign(sdeModel@parameter@all[i], true.parameter[[j]],yuimaEnv)
        }
      }
      #assign(sdeModel@parameter@all[i], true.parameter[[i]], yuimaEnv)
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
  if(!missing(increment.W) & !is.null(increment.W)){
    if(space.discretized == TRUE){
      yuima.warn("Parameter increment must be invalid if space.discretized=TRUE.")
      return(NULL)
    }else if(dim(increment.W)[1] != r.size){
      yuima.warn("Length of increment's row must be same as yuima@model@noise.number.")
      return(NULL)
    }else if(dim(increment.W)[2] != n){
      yuima.warn("Length of increment's column must be same as sampling@n[1].")
      return(NULL)
    }
  }
  
  ##:: Error check for increment specified version.
  if(!missing(increment.L) & !is.null(increment.L)){
    if(space.discretized == TRUE){
      yuima.warn("Parameter increment must be invalid if space.discretized=TRUE.")
      return(NULL)
    }else if(dim(increment.L)[1] != length(yuima@model@jump.coeff[[1]]) ){ #r.size){
      yuima.warn("Length of increment's row must be same as yuima@model@noise.number.")
      return(NULL)
    }else if(dim(increment.L)[2] != n){
      yuima.warn("Length of increment's column must be same as sampling@n[1].")
      return(NULL)
    }
  }
  
  
  yuimaEnv$dL <- increment.L
  
  
  if(space.discretized){   	  
    ##:: using Space-discretized Euler-Maruyama method	  
    yuima@data <- space.discretized(xinit, yuima, yuimaEnv)
    
    yuima@model@hurst<-tmphurst
    return(yuima)  
  }
  
  
  ##:: using Euler-Maruyama method
  delta <- Terminal/n 
  
  if(missing(increment.W) | is.null(increment.W)){
    
    if( sdeModel@hurst!=0.5 ){ 
      
      grid<-sampling2grid(yuima@sampling)	
      isregular<-yuima@sampling@regular 
      
      if((!isregular) || (methodfGn=="Cholesky")){
        dW<-CholeskyfGn(grid, sdeModel@hurst,r.size)
        yuima.warn("Cholesky method for simulating fGn has been used.")
      } else {
        dW<-WoodChanfGn(grid, sdeModel@hurst,r.size)
      }
      
    } else {
      
      delta<-Terminal/n
      if(!is.Poisson(sdeModel)){ # if pure CP no need to setup dW
        dW <- rnorm(n * r.size, 0, sqrt(delta))
        dW <- matrix(dW, ncol=n, nrow=r.size,byrow=TRUE)
      } else {
        dW <- matrix(0,ncol=n,nrow=1)  # maybe to be fixed
      }
    }
    
  } else {
    dW <- increment.W
  }
  
  if(is.Poisson(sdeModel)){
    yuima@data <- simCP(xinit, yuima, yuimaEnv)
  } else {
    yuima@data <- euler(xinit, yuima, dW, yuimaEnv)
  }
  
  for(i in 1:length(yuima@data@zoo.data)) 
    index(yuima@data@zoo.data[[i]]) <- yuima@sampling@grid[[1]]  ## to be fixed
  yuima@model@xinit <- xinit
  yuima@model@hurst <-tmphurst
  
  if(missing(subsampling))
    return(yuima)
  subsampling(yuima, subsampling)
  
}

aux.simulateCogarch<-function(object, nsim, seed, xinit, true.parameter, 
                              space.discretized, increment.W, increment.L,
                              hurst,methodfGn,
                              sampling, subsampling,
                              Initial, Terminal, n, delta, 
                              grid, random, sdelta, 
                              sgrid, interpolation){
  yuimaCogarch<-object
  model<-yuimaCogarch@model
  info<-model@info
  samp <- yuimaCogarch@sampling
    aux.Noise<-setModel(drift="1",
                        diffusion="0",
                        jump.coeff="1",
                        measure=info@measure,
                        measure.type=info@measure.type)
  
#   aux.samp<-setSampling(Initial = samp@Initial, Terminal = samp@Terminal[1], n = samp@n[1], delta = samp@delta, 
#                         grid=samp@grid, random = samp@random, sdelta=samp@sdelta, 
#                         sgrid=samp@sgrid, interpolation=samp@interpolation )

    aux.samp<-setSampling(Initial = samp@Initial, Terminal = samp@Terminal[1], n = samp@n[1])
    auxModel<-setYuima(model=aux.Noise, sampling= aux.samp)

  if(length(model@parameter@measure)==0){
    aux.incr2<-aux.simulate(object=auxModel, nsim=nsim, seed=seed, 
                           space.discretized=space.discretized, increment.W=increment.W, 
                           increment.L=increment.L,
                           hurst=0.5,methodfGn=methodfGn)
  }else{
    aux.incr2<-aux.simulate(object=auxModel, nsim=nsim, seed=seed, 
                            true.parameter = true.parameter[model@parameter@measure], 
                            space.discretized=space.discretized, increment.W=increment.W, 
                            increment.L=increment.L,
                            hurst=0.5,methodfGn=methodfGn)
  }    
  increment<-diff(as.numeric(get.zoo.data(aux.incr2)[[1]]))
  # Using the simulated increment for generating the quadratic variation
  # As first step we compute it in a crude way. A more fine approach is based on 
  # the mpv function.
  quadratVariation <- increment^2
  incr.L <- t(matrix(c(increment,quadratVariation),ncol=2))
  incr.W <- matrix(0, nrow=1,ncol=length(increment))
  # Simulate trajectories Cogarch
  result <- aux.simulate(object=yuimaCogarch, nsim=nsim, seed=seed, xinit=xinit,
                         true.parameter = true.parameter, 
                         space.discretized = space.discretized,increment.W =incr.W,
                         increment.L=incr.L,
                    hurst=hurst,methodfGn=methodfGn,
                    sampling=sampling, subsampling=subsampling,
                    Initial=Initial, Terminal=Terminal, n=n, delta=delta, 
                    grid=grid, random=random, sdelta=sdelta, 
                    sgrid=sgrid, interpolation=interpolation)
  return(result)
}

