setMethod("simulate", "yuima.Hawkes",
          function(object, nsim=1, seed=NULL, xinit, true.parameter,
                   space.discretized=FALSE, increment.W=NULL, increment.L=NULL, method="euler",
                   hurst, methodfGn="WoodChan",
                   sampling, subsampling,
                   #Initial = 0, Terminal = 1, n = 100, delta,
                   #	grid, random = FALSE, sdelta=as.numeric(NULL),
                   #	sgrid=as.numeric(NULL), interpolation="none"
                   ...){
            res <- aux.simulatHawkes(object, nsim = nsim, seed = seed,
                                     xinit = xinit, true.parameter = true.parameter,
                                     space.discretized = space.discretized, increment.W = increment.W,
                                     increment.L = increment.L, method = method, hurst = hurst,
                                     methodfGn = methodfGn, sampling = sampling, subsampling = subsampling)

            return(res)
          }
)

aux.simulatHawkes<- function(object, nsim, seed,
  xinit, true.parameter, space.discretized, increment.W,
  increment.L, method, hurst, methodfGn, sampling, subsampling){
  # Here we can construct specific algorithm for the standard Hawkes process
  res <- aux.simulatPPR(object, nsim = nsim, seed = seed,
                        xinit = xinit, true.parameter = true.parameter,
                        space.discretized = space.discretized, increment.W = increment.W,
                        increment.L = increment.L, method = method, hurst = hurst,
                        methodfGn = methodfGn, sampling = sampling, subsampling = subsampling)

  return(res)

#   object@Kernel@param.Integral@allparam
#   simOzaki.aux(gFun=object@gFun@formula,a,cCoeff, Time, numJump)
}

setMethod("simulate", "yuima.PPR",
          function(object, nsim=1, seed=NULL, xinit, true.parameter,
                   space.discretized=FALSE, increment.W=NULL, increment.L=NULL, method="euler",
                   hurst, methodfGn="WoodChan",
                   sampling, subsampling,
                   #Initial = 0, Terminal = 1, n = 100, delta,
                   #	grid, random = FALSE, sdelta=as.numeric(NULL),
                   #	sgrid=as.numeric(NULL), interpolation="none"
                   ...){
            res <- aux.simulatPPR(object, nsim = nsim, seed = seed,
                                     xinit = xinit, true.parameter = true.parameter,
                                     space.discretized = space.discretized, increment.W = increment.W,
                                     increment.L = increment.L, method = method, hurst = hurst,
                                     methodfGn = methodfGn, sampling = sampling, subsampling = subsampling)

            return(res)
          }
)

constHazIntPr <- function(g.Fun , Kern.Fun, covariates, counting.var,statevar=NULL){
  numb.Int <- length(g.Fun)
  Int.Intens <- list()
  dum.g <- character(length=numb.Int)
  for(i in c(1:numb.Int)){
    dum.g0 <- as.character(g.Fun[i])
    dum.g0 <- gsub("(", "( ", fixed=TRUE,x = dum.g0)
    dum.g0 <- gsub(")", " )", fixed=TRUE,x = dum.g0)
    
    if(length(counting.var)>0){  
      for(j in c(1:length(counting.var))){
        my.countOld <- paste0(counting.var[j] ," ")
        #my.countNew <- paste0("as.numeric(", counting.var[i] ,")")
        my.countNew <- paste0("(", counting.var[j] ,")")
        dum.g0 <- gsub(my.countOld, my.countNew, x = dum.g0, fixed=TRUE)
        my.countOld <- paste0(counting.var[j] ,"[",Kern.Fun@variable.Integral@upper.var,"]")
        
        my.countNew <- paste0("(", counting.var[j] ,")")
        dum.g0 <- gsub(my.countOld, my.countNew, x = dum.g0, fixed=TRUE)
      }
    }
    if(length(covariates)>0){
      for(j in c(1:length(covariates))){
        my.covarOld <- paste0(covariates[j] ," ")
        my.covarNew <-  covariates[j]
        dum.g0 <- gsub(my.covarOld, my.covarNew, x = dum.g0, fixed=TRUE)
        my.covarOld <- paste0(covariates[j] ,"[",Kern.Fun@variable.Integral@upper.var,"]")
        my.covarNew <- covariates[j] 
        dum.g0 <- gsub(my.covarOld, my.covarNew, x = dum.g0, fixed=TRUE)
      }
    }
    
    dum.g[i] <- paste("tail(",dum.g0,", n=1L)")
    
  }
  dum.Ker <- as.character(unlist(Kern.Fun@Integrand@IntegrandList))
  dum.Ker <- gsub("(", "( ", fixed=TRUE,x = dum.Ker)
  dum.Ker <- gsub(")", " )", fixed=TRUE,x = dum.Ker)
  
  dimKernel<-length(dum.Ker)
  condIntInKer <- FALSE
  for(i in c(1:dimKernel)){
    if(!condIntInKer)
      condIntInKer <- statevar%in%all.vars(parse(text=dum.Ker[i]))
  }
  if(condIntInKer){
    for(i in c(1:length(statevar))){
      my.countOld <- paste0(statevar[i] ," ")
      my.countNew <- paste0( statevar[i] ,
                             "ForKernel[CondJumpGrid]")
      dum.Ker <- gsub(my.countOld, my.countNew, x = dum.Ker, fixed=TRUE)
      my.countOld <- paste0(statevar[i] ,"[",Kern.Fun@variable.Integral@upper.var,"]")
      # my.countNew <- paste0( counting.var[i] ,
      #                       "[ as.character( ",Kern.Fun@variable.Integral@upper.var ," ) ]")
      my.countNew <- paste0( "tail(",statevar[i] ,"ForKernel ,n=1L) ")
      dum.Ker <- gsub(my.countOld, my.countNew, x = dum.Ker, fixed=TRUE)
      
      my.countOld <- paste0(statevar[i] ,"[",Kern.Fun@variable.Integral@var.time,"]")
      my.countNew <- paste0(statevar[i] ,
                            "ForKernel[CondJumpGrid]")
      dum.Ker <- gsub(my.countOld, my.countNew, x = dum.Ker, fixed=TRUE)
      
    }
  }
  
    
  if(length(counting.var)>0){  
    for(i in c(1:length(counting.var))){
      my.countOld <- paste0(counting.var[i] ," ")
      my.countNew <- paste0( counting.var[i] ,
                              "[CondJumpGrid]")
      dum.Ker <- gsub(my.countOld, my.countNew, x = dum.Ker, fixed=TRUE)
      my.countOld <- paste0(counting.var[i] ,"[",Kern.Fun@variable.Integral@upper.var,"]")
        # my.countNew <- paste0( counting.var[i] ,
        #                       "[ as.character( ",Kern.Fun@variable.Integral@upper.var ," ) ]")
      my.countNew <- paste0( "tail(",counting.var[i] ,",n=1L) ")
      dum.Ker <- gsub(my.countOld, my.countNew, x = dum.Ker, fixed=TRUE)
        
      my.countOld <- paste0(counting.var[i] ,"[",Kern.Fun@variable.Integral@var.time,"]")
      my.countNew <- paste0(counting.var[i] ,
                              "[CondJumpGrid]")
      dum.Ker <- gsub(my.countOld, my.countNew, x = dum.Ker, fixed=TRUE)
        
    }
  }
  if(length(covariates)>0){
    for(i in c(1:length(covariates))){
        
      my.countOld <- paste0(covariates[i] ," ")
      my.countNew <- paste0( covariates[i] ,
                            "[CondJumpGrid]")
      dum.Ker <- gsub(my.countOld, my.countNew, x = dum.Ker, fixed=TRUE)
        
      my.countOld <- paste0(covariates[i] ,"[",Kern.Fun@variable.Integral@upper.var,"]")
      my.countNew <- paste0("tail(", covariates[i] , ", n=1L ) ")
      dum.Ker <- gsub(my.countOld, my.countNew, x = dum.Ker, fixed=TRUE)
        
      my.countOld <- paste0(covariates[i] ,"[",Kern.Fun@variable.Integral@var.time,"]")
      my.countNew <- paste0( covariates[i] ,
                              "[CondJumpGrid]")
      dum.Ker <- gsub(my.countOld, my.countNew, x = dum.Ker, fixed=TRUE)
    }
  }
  dif.dx <- paste("d",Kern.Fun@variable.Integral@var.dx, sep="")
  if(Kern.Fun@Integrand@dimIntegrand[1]==1){
    dum.Ker <- paste(dum.Ker,dif.dx, sep = "*")
  }else{
      # dum.Ker <- paste(dum.Ker,rep(dif.dx, Kern.Fun@Integrand@dimIntegrand[1]), sep = "*")
    dum.Ker <- matrix(dum.Ker,Kern.Fun@Integrand@dimIntegrand[1],
    Kern.Fun@Integrand@dimIntegrand[2], byrow=T)
    dum.Ker <- paste(dum.Ker,dif.dx, sep = "*")
    dum.Ker <- matrix(dum.Ker,Kern.Fun@Integrand@dimIntegrand[1],
                      Kern.Fun@Integrand@dimIntegrand[2], byrow=T)
  }
  cond.Sup <- paste(Kern.Fun@variable.Integral@var.time, "<", Kern.Fun@variable.Integral@upper.var)
  dum.Ker <- paste("(",dum.Ker, ") * (", cond.Sup, ")")
  dum.Ker <- paste0("sum(",dum.Ker,")")
  if(Kern.Fun@Integrand@dimIntegrand[2]>1 & Kern.Fun@Integrand@dimIntegrand[1]==1){
    dum.Ker <- paste(dum.Ker,collapse = " + ")
  }
  if(Kern.Fun@Integrand@dimIntegrand[1]>1){
    mydum <- matrix(dum.Ker,Kern.Fun@Integrand@dimIntegrand[1],
                      Kern.Fun@Integrand@dimIntegrand[2])
    dum.Ker <- character(length = Kern.Fun@Integrand@dimIntegrand[1])
    for (i in c(1:Kern.Fun@Integrand@dimIntegrand[1])){
      dum.Ker[i] <- paste(mydum[i,],collapse = " + ")
    }
      #yuima.stop("Check")
      
  }
    # dum.Ker <- paste("(",dum.Ker,") * (")
    # cond.Sup <- paste(Kern.Fun@variable.Integral@var.time, "<", Kern.Fun@variable.Integral@upper.var)
    # dum.Ker <- paste(dum.Ker, cond.Sup, ")")
  for(i in c(1:numb.Int)){   
    Int.Intens[[i]] <- parse(text = paste(dum.g[i], dum.Ker[i], sep = " + "))
  }
  res <- list(Intens = Int.Intens)
}




aux.simulatPPR<- function(object, nsim = nsim, seed = seed,
               xinit = xinit, true.parameter = true.parameter,
               space.discretized = space.discretized, increment.W = increment.W,
               increment.L = increment.L, method = method, hurst = hurst,
               methodfGn = methodfGn, sampling = sampling,
               subsampling = subsampling){
  ROLDVER<-!(is(object@model@measure$df,"yuima.law"))
  if(ROLDVER){
    object <- aux.simulatPPRROldVersion(object, nsim = nsim, seed = seed,
                                          xinit = xinit, true.parameter = true.parameter,
                                          space.discretized = space.discretized, increment.W = increment.W,
                                          increment.L = increment.L, method = method, hurst = hurst,
                                          methodfGn = methodfGn, sampling = sampling,
                                          subsampling = subsampling)
  }else{
    if(object@PPR@RegressWithCount){
      object <-aux.simulatPPRWithCount(object, nsim = nsim, seed = seed,
                                       xinit = xinit, true.parameter = true.parameter,
                                       space.discretized = space.discretized, increment.W = increment.W,
                                       increment.L = increment.L, method = method, hurst = hurst,
                                       methodfGn = methodfGn, sampling = sampling,
                                       subsampling = subsampling)
      
    }else{
        posLambda <- object@model@state.variable %in% object@PPR@additional.info
        condInteFeedbackCov <- any(posLambda) 
      if(condInteFeedbackCov){
        object <- aux.simulatPPRWithIntesFeedBack(object, nsim = nsim, seed = seed,
                                                  xinit = xinit, true.parameter = true.parameter,
                                                  space.discretized = space.discretized, increment.W = increment.W,
                                                  increment.L = increment.L, method = method, hurst = hurst,
                                                  methodfGn = methodfGn, sampling = sampling,
                                                  subsampling = subsampling,
                                                  posLambda=posLambda)
        
      }else{
       object <- aux.simulatPPRROldNew(object, nsim = nsim, seed = seed,
                                        xinit = xinit, true.parameter = true.parameter,
                                        space.discretized = space.discretized, increment.W = increment.W,
                                        increment.L = increment.L, method = method, hurst = hurst,
                                        methodfGn = methodfGn, sampling = sampling,
                                        subsampling = subsampling)
      }
    }
  }
  return(object)
}

eulerPPR<-function(xinit,yuima,Initial,Terminal, dW, n, env){
  
  sdeModel<-yuima@model
  
  modelstate <- sdeModel@solve.variable
  modeltime <- sdeModel@time.variable
  V0 <- sdeModel@drift
  V <- sdeModel@diffusion
  r.size <- sdeModel@noise.number
  d.size <- sdeModel@equation.number
  # Terminal <- yuima@sampling@Terminal[1]
  # Initial <- yuima@sampling@Initial[1]
  
  #n <- ceiling((Terminal-Initial)/yuima@sampling@delta)
  dL <- env$dL
  
  if(length(unique(as.character(xinit)))==1 &&
     is.numeric(tryCatch(eval(xinit[1],env),error=function(...) FALSE))){
    dX_dummy<-xinit[1]
    dummy.val<-eval(dX_dummy, env)
    if(length(dummy.val)==1){dummy.val<-rep(dummy.val,length(xinit))}
    for(i in 1:length(modelstate)){
      assign(modelstate[i],dummy.val[i] ,env)
    }
    dX<-vector(mode="numeric",length(dX_dummy))
    
    for(i in 1:length(xinit)){
      dX[i] <- dummy.val[i]
    }
  }else{
    dX_dummy <- xinit
    if(length(modelstate)==length(dX_dummy)){
      for(i in 1:length(modelstate)) {
        if(is.numeric(tryCatch(eval(dX_dummy[i],env),error=function(...) FALSE))){
          assign(modelstate[i], eval(dX_dummy[i], env),env)
        }else{
          assign(modelstate[i], 0, env)
        }
      }
    }else{
      yuima.warn("the number of model states do not match the number of initial conditions")
      return(NULL)
    }
    
    dX<-vector(mode="numeric",length(dX_dummy))
    
    for(i in 1:length(dX_dummy)){
      dX[i] <- eval(dX_dummy[i], env)
    }
  }
  
  delta <- yuima@sampling@delta
  
  if(!length(sdeModel@measure.type)){
    
    b <- parse(text=paste("c(",paste(as.character(V0),collapse=","),")"))
    vecV <- parse(text=paste("c(",paste(as.character(unlist(V)),collapse=","),")"))
    
    X_mat <- .Call("euler", dX, Initial, as.integer(r.size), 
                   rep(1, n) * delta, dW, modeltime, modelstate, quote(eval(b, env)), 
                   quote(eval(vecV, env)), env, new.env())
    tsX <- ts(data=t(X_mat), deltat=delta , start = Initial) #LM
  }else{
    has.drift <- sum(as.character(sdeModel@drift) != "(0)")
    var.in.diff <- is.logical(any(match(unlist(lapply(sdeModel@diffusion, all.vars)), sdeModel@state.variable)))
    
    p.b <- function(t, X=numeric(d.size)){
      for(i in 1:length(modelstate)){
        assign(modelstate[i], X[i], env)
      }
      assign(modeltime, t, env)
      if(has.drift){
        tmp <- matrix(0, d.size, r.size+1)
        for(i in 1:d.size){
          tmp[i,1] <- eval(V0[i], env)
          for(j in 1:r.size){
            tmp[i,j+1] <- eval(V[[i]][j],env)
          }
        }
      } else {  ##:: no drift term (faster)
        tmp <- matrix(0, d.size, r.size)
        if(!is.Poisson(sdeModel)){ # we do not need to evaluate diffusion
          for(i in 1:d.size){
            for(j in 1:r.size){
              tmp[i,j] <- eval(V[[i]][j],env)
            } # for j
          } # foh i
        } # !is.Poisson
      } # else
      return(tmp)
    }
    
    X_mat <- matrix(0, d.size, (n+1))
    X_mat[,1] <- dX
    
    if(has.drift){  
      dW <- rbind( rep(1, n)*delta , dW)
    }
    
    JP <- sdeModel@jump.coeff
    mu.size <- length(JP)
    
    p.b.j <- function(t, X=numeric(d.size)){
      for(i in 1:length(modelstate)){
        assign(modelstate[i], X[i], env)
      }
      assign(modeltime, t, env)
      
      j.size <- length(JP[[1]])
      tmp <- matrix(0, mu.size, j.size)
      
      for(i in 1:mu.size){
        for(j in 1:j.size){
          tmp[i,j] <- eval(JP[[i]][j],env)
        }
      }
      return(tmp)
    }
    
    dZ <- dL
    
    if(is.null(dim(dZ)))
      dZ <- matrix(dZ,nrow=1)
    for(i in 1:n){
      # if(i==720 & n==720){
      #   aa<-NULL
      # }
      assign(sdeModel@jump.variable, dZ[,i], env)
      
      if(sdeModel@J.flag){
        dZ[,i] <- 1
      }
      tmp.j <- p.b.j(t=Initial+(i - 1)*delta, X=dX) # YK
      if(sum(dim(tmp.j))==2)
        tmp.j <- as.numeric(tmp.j)
      dX <- dX + p.b(t=Initial+(i - 1)*delta, X=dX) %*% dW[, i] +tmp.j %*% dZ[,i] # YK
      X_mat[, i+1] <- dX
    }
    #tsX <- ts(data=t(X_mat), deltat=delta, start=yuima@sampling@Initial)
    
  }
  return(X_mat)
}

eulerPPRwithInt<-function(xinit,yuima,Initial,Terminal, dW, dL, n, env){
  
  sdeModel<-yuima@model
  
  modelstate <- sdeModel@solve.variable
  modeltime <- sdeModel@time.variable
  V0 <- sdeModel@drift
  V <- sdeModel@diffusion
  r.size <- sdeModel@noise.number
  d.size <- sdeModel@equation.number
  # Terminal <- yuima@sampling@Terminal[1]
  # Initial <- yuima@sampling@Initial[1]
  
  #n <- ceiling((Terminal-Initial)/yuima@sampling@delta)
  #dL <- env$dL
  
  if(length(unique(as.character(xinit)))==1 &&
     is.numeric(tryCatch(eval(xinit[1],env),error=function(...) FALSE))){
    dX_dummy<-xinit[1]
    dummy.val<-eval(dX_dummy, env)
    if(length(dummy.val)==1){dummy.val<-rep(dummy.val,length(xinit))}
    for(i in 1:length(modelstate)){
      assign(modelstate[i],dummy.val[i] ,env)
    }
    dX<-vector(mode="numeric",length(dX_dummy))
    
    for(i in 1:length(xinit)){
      dX[i] <- dummy.val[i]
    }
  }else{
    dX_dummy <- xinit
    if(length(modelstate)==length(dX_dummy)){
      for(i in 1:length(modelstate)) {
        if(is.numeric(tryCatch(eval(dX_dummy[i],env),error=function(...) FALSE))){
          assign(modelstate[i], eval(dX_dummy[i], env),env)
        }else{
          assign(modelstate[i], 0, env)
        }
      }
    }else{
      yuima.warn("the number of model states do not match the number of initial conditions")
      # return(NULL)
      dX_dummy <- c(dX_dummy,env[[yuima@PPR@additional.info]])
      for(i in 1:length(modelstate)){
        if(is.numeric(tryCatch(eval(dX_dummy[i],env),error=function(...) FALSE))){
          assign(modelstate[i], eval(dX_dummy[i], env),env)
        }else{
          assign(modelstate[i], 0, env)
        }        
      }
    }
    
    dX<-vector(mode="numeric",length(dX_dummy))
    
    for(i in 1:length(dX_dummy)){
      dX[i] <- eval(dX_dummy[i], env)
    }
  }
  
  delta <- yuima@sampling@delta
  
  if(!length(sdeModel@measure.type)){
    
    b <- parse(text=paste("c(",paste(as.character(V0),collapse=","),")"))
    vecV <- parse(text=paste("c(",paste(as.character(unlist(V)),collapse=","),")"))
    
    X_mat <- .Call("euler", dX, Initial, as.integer(r.size), 
                   rep(1, n) * delta, dW, modeltime, modelstate, quote(eval(b, env)), 
                   quote(eval(vecV, env)), env, new.env())
    tsX <- ts(data=t(X_mat), deltat=delta , start = Initial) #LM
  }else{
    has.drift <- sum(as.character(sdeModel@drift) != "(0)")
    var.in.diff <- is.logical(any(match(unlist(lapply(sdeModel@diffusion, all.vars)), sdeModel@state.variable)))
    
    p.b <- function(t, X=numeric(d.size)){
      for(i in 1:length(modelstate)){
        assign(modelstate[i], X[i], env)
      }
      assign(modeltime, t, env)
      if(has.drift){
        tmp <- matrix(0, d.size, r.size+1)
        for(i in 1:d.size){
          tmp[i,1] <- eval(V0[i], env)
          for(j in 1:r.size){
            tmp[i,j+1] <- eval(V[[i]][j],env)
          }
        }
      } else {  ##:: no drift term (faster)
        tmp <- matrix(0, d.size, r.size)
        if(!is.Poisson(sdeModel)){ # we do not need to evaluate diffusion
          for(i in 1:d.size){
            for(j in 1:r.size){
              tmp[i,j] <- eval(V[[i]][j],env)
            } # for j
          } # foh i
        } # !is.Poisson
      } # else
      return(tmp)
    }
    
    X_mat <- matrix(0, d.size, (n+1))
    X_mat[,1] <- dX
    
    if(has.drift){  
      dW <- rbind( rep(1, n)*delta , dW)
    }
    
    JP <- sdeModel@jump.coeff
    mu.size <- length(JP)
    
    p.b.j <- function(t, X=numeric(d.size)){
      for(i in 1:length(modelstate)){
        assign(modelstate[i], X[i], env)
      }
      assign(modeltime, t, env)
      
      j.size <- length(JP[[1]])
      tmp <- matrix(0, mu.size, j.size)
      
      for(i in 1:mu.size){
        for(j in 1:j.size){
          tmp[i,j] <- eval(JP[[i]][j],env)
        }
      }
      return(tmp)
    }
    
    dZ <- dL
    
    if(is.null(dim(dZ)))
      dZ <- matrix(dZ,nrow=1)
    for(i in 1:n){
      # if(i==720 & n==720){
      #   aa<-NULL
      # }
      assign(sdeModel@jump.variable, dZ[,i], env)
      
      if(sdeModel@J.flag){
        dZ[,i] <- 1
      }
      tmp.j <- p.b.j(t=Initial+(i - 1)*delta, X=dX) # YK
      if(sum(dim(tmp.j))==2)
        tmp.j <- as.numeric(tmp.j)
      dX <- dX + p.b(t=Initial+(i - 1)*delta, X=dX) %*% dW[, i] +tmp.j %*% dZ[,i] # YK
      X_mat[, i+1] <- dX
    }
    #tsX <- ts(data=t(X_mat), deltat=delta, start=yuima@sampling@Initial)
    
  }
  return(X_mat)
}


aux.simulatPPRWithCount<-function(object, nsim = nsim, seed = seed,
                                xinit = xinit, true.parameter = true.parameter,
                                space.discretized = space.discretized, increment.W = increment.W,
                                increment.L = increment.L, method = method, hurst = 0.5,
                                methodfGn = methodfGn, sampling = sampling,
                                subsampling = subsampling){
  
  samp <- sampling
  Model <- object@model
  gFun <- object@gFun
  Kern <- object@Kernel
  object@sampling <- samp
  randomGenerator<-object@model@measure$df
  
  if(missing(increment.W) | is.null(increment.W)){
    
    if( Model@hurst!=0.5 ){
      
      grid<-sampling2grid(object@sampling)
      isregular<-object@sampling@regular
      
      if((!isregular) || (methodfGn=="Cholesky")){
        dW<-CholeskyfGn(grid, Model@hurst,Model@noise.number)
        yuima.warn("Cholesky method for simulating fGn has been used.")
      } else {
        dW<-WoodChanfGn(grid, Model@hurst,Model@noise.number)
      }
      
    } else {
      
      # delta<-(Terminal-Initial)/n
      delta <- samp@delta
      if(!is.Poisson(Model)){ # if pure CP no need to setup dW
        dW <- rnorm(samp@n * Model@noise.number, 0, sqrt(delta))
        dW <- matrix(dW, ncol=samp@n, nrow=Model@noise.number,byrow=TRUE)
      } else {
        dW <- matrix(0,ncol=samp@n,nrow=1)  # maybe to be fixed
      }
    }
    
  } else {
    dW <- increment.W
  }
  if(missing(xinit)){
    if(length(object@model@xinit)!=0){
      xinit<-numeric(length=length(object@model@xinit))
      for(i in c(1:object@model@equation.number))
        xinit[i] <- eval(object@model@xinit[i])
    }else{
      xinit <- rep(0,object@model@equation.number)
      object@model@xinit<-xinit
    }
  }
  if(missing(hurst)){
    hurst<-0.5
  }
  if(samp@regular){
    tForMeas<-samp@delta
    NumbIncr<-samp@n
    if(missing(true.parameter)){
      eval(parse(text= paste0("measureparam$",
                              object@model@time.variable," <- tForMeas",collapse="")))
    }else{
      measureparam<-true.parameter[object@model@parameter@measure]
      eval(parse(text= paste0("measureparam$",
                              object@model@time.variable," <- tForMeas",collapse="")))
      
    }
    Noise.L <- t(rand(object = randomGenerator, n=NumbIncr, param=measureparam))
    rownames(Noise.L)<-Model@solve.variable
    #dIncr <- apply(cbind(0,Noise.L),1,diff)
    Noise.count <- Noise.L[object@PPR@counting.var,]
    Noise.Laux <- Noise.L
    for(i in c(1:length(object@PPR@counting.var))){
      Noise.Laux[object@PPR@counting.var[i],]<-0
    }  
  }
    myenv<-new.env()
    par.len <- length(object@PPR@allparam)
    if(par.len>0){
      for(i in 1:par.len){
        pars <- object@PPR@allparam[i]
        
        for(j in 1:length(true.parameter)){
          if( is.na(match(pars, names(true.parameter)[j]))!=TRUE){
            assign(object@PPR@allparam[i], true.parameter[[j]],myenv)
          }
        }
      }
    }
    assign("dL",Noise.Laux,myenv)
    
    
    condMatrdW <- is.matrix(dW)
    if(condMatrdW){
      dimdW <- dim(dW)[2]
    }else{
      dimdW <- length(dW)
    }
    
    CovariateSim<- eulerPPR(xinit=xinit,yuima=object,dW=dW, 
             Initial=samp@Initial,Terminal=samp@Terminal,n=samp@n,
             env=myenv)
    rownames(CovariateSim)<- Model@solve.variable
    assign("info.PPR", object@PPR, myenv)
    
    dimCov <- length(object@PPR@covariates)
    if (dimCov>0){
      for(j in c(1:dimCov)){
        assign(object@PPR@covariates[j],
               as.numeric(CovariateSim[object@PPR@covariates[j],1]),
               envir = myenv)
      }
    }
    
    dimNoise<-dim(Noise.Laux)
    dimCovariateSim <- dim(CovariateSim)
    
    
    ExprHaz <- constHazIntPr(g.Fun = object@gFun@formula,
       Kern.Fun = object@Kernel, covariates = object@PPR@covariates,
       counting.var = object@PPR@counting.var,
       statevar = object@model@state.variable)$Intens
    
    
    # Start Simulation PPR
    compErrHazR4 <- function(samp, Kern,
                             capitalTime, Model, 
                             my.env, ExprHaz, Time, dN, 
                             Index, pos){
      assign(Kern@variable.Integral@var.time, Time, envir = my.env)
      
      assign(Model@time.variable, capitalTime, envir = my.env)
      l <- 1
      for(i in c(1:length(Kern@variable.Integral@var.dx)) ){
        if(any(Kern@variable.Integral@var.dx[i]==my.env$info.PPR@counting.var)){
          assign(paste0("d",Kern@variable.Integral@var.dx[i]), dN[l,], envir =my.env)
          l <- l + 1
        }
        if(Kern@variable.Integral@var.dx[i]%in%my.env$info.PPR@covariates){  
          assign(paste0("d",Kern@variable.Integral@var.dx[i]),
                 diff(c(0,my.env[[Kern@variable.Integral@var.dx[i]]])) , 
                 envir =my.env)
        }
      }
      condPointIngrid <- samp@grid[[1]]<=my.env$t
      PointIngridInt <- samp@grid[[1]][condPointIngrid]
      CondJumpGrid <- PointIngridInt %in% my.env$s
      assign("CondJumpGrid", CondJumpGrid, envir = my.env)
      Lambda <- NULL
      #  for(h in c(1:Index)){
      #      Lambdadum <- eval(ExprHaz[[h]], envir = my.env)
      #      Lambda <- rbind(Lambda,Lambdadum)
      #  
      Lambda <- eval(ExprHaz[[pos]], envir = my.env)
      #    rownames(Lambda) <- my.env$info.PPR@counting.var
      return(Lambda)
      
    }
    dN <- matrix(0,object@gFun@dimension[1],object@gFun@dimension[2])
    grid <- samp@grid[[1]]
    const <- -log(runif(gFun@dimension[1]))
    condMyTR <- const<delta
    while(any(condMyTR)){
      if(sum(condMyTR)==0){
        const <- -log(runif(length(condMyTR)))
        condMyTR <- const<delta
      }else{
        const[condMyTR] <- -log(runif(sum(condMyTR)))
        condMyTR <- const<delta
      }
    }
    jumpT<-NULL
    i <- 1
    dimGrid <-length(grid)
    cond <- const
    Index <- gFun@dimension[1]
    inter_i <- rep(i,Index)
    noExit<-rep(T,Index)
    while(any(noExit)){
      
      for(j in c(1:Index)){
        HazardRate<-0
        while(cond[j]>0 && noExit[j]){
          lambda<-compErrHazR4(samp, Kern, capitalTime=samp@grid[[1]][inter_i[j]], 
                               Model, myenv, ExprHaz, Time=jumpT, dN, Index, j)
    
          incrlambda <- lambda*delta
          HazardRate <- HazardRate+incrlambda
          cond[j] <- const[j]-HazardRate
          inter_i[j]<-inter_i[j]+1
          if(inter_i[j]>=(dimGrid-1)){
            noExit[j] <- FALSE
          }
          if(inter_i[j]<dim(CovariateSim)[2]){  
            dimCov <- length(object@PPR@covariates)
            
            if (dimCov>0){
              for(j in c(1:dimCov)){
                assign(object@PPR@covariates[j],
                       as.numeric(CovariateSim[object@PPR@covariates[j],1:inter_i[j]]),
                       envir = myenv)
              }
            }  
          }
        }
      }
      
      i <- min(inter_i)
      
      if(any(noExit)){
        if(i<dim(CovariateSim)[2]){ 
          jumpT<-c(jumpT,grid[i])
          
          if(dim(dN)[2]==1 & all(dN[,1]==0)){
            dN[i==inter_i,1] <- Noise.count[i-1]
            Noise.Laux[object@PPR@counting.var,i-1]<-Noise.count[i-1]
            dumdN <- dN
          }else{
            dumdN <- rep(0,Index)
            dumdN[i==inter_i] <- Noise.count[i-1] 
            Noise.Laux[object@PPR@counting.var,i-1] <- dumdN[i==inter_i] 
            dN <- cbind(dN,dumdN)
          }
          # cat("\n ", i, grid[i])
          
          # assign("dL",Noise.Laux,myenv)
          # 
          # CovariateSim<- eulerPPR(xinit=xinit,yuima=object,dW=dW,
          #                         Initial=samp@Initial,Terminal=samp@Terminal,
          #                         env=myenv)
          
          assign("dL",Noise.Laux[,c((i-1):dimNoise[2])],myenv)
          xinit <- CovariateSim[,i-1]


          if(condMatrdW){
            CovariateSim[,(i-1):dimCovariateSim[2]] <- eulerPPR(xinit=xinit,
              yuima=object,dW=dW[,(i-1):dimdW],
              Initial=samp@grid[[1]][i-1],Terminal=samp@Terminal,n=(samp@n-(i-1)+1),
              env=myenv)
          }else{
            CovariateSim[,(i-1):dimCovariateSim[2]] <- eulerPPR(xinit=xinit,
              yuima=object, dW=dW[(i-1):dimdW],
             Initial=samp@grid[[1]][i-1],Terminal=samp@Terminal,n=(samp@n-(i-1)+1),
             env=myenv)
          }
          
          rownames(CovariateSim)<- Model@solve.variable
          
          
          const <- -log(runif(object@gFun@dimension[1]))
          condMyTR <- const<delta
          while(any(condMyTR)){
            if(sum(condMyTR)==0){
              const <- -log(runif(length(condMyTR)))
              condMyTR <- const<delta
            }else{
              const[condMyTR] <- -log(runif(sum(condMyTR)))
              condMyTR <- const<delta
            }
          }
          
          cond <- const
          if(all(noExit)){
            inter_i <- rep(i, Index)
          }else{
            if(any(noExit)){
              inter_i[noExit] <- i
              inter_i[!noExit] <- samp@n+1 
            }
          }
          
        }
      }  
    }
    tsX <- ts(data=t(CovariateSim), deltat=delta, start=object@sampling@Initial)
    object@data <- setData(original.data=tsX)
    for(i in 1:length(object@data@zoo.data))
      index(object@data@zoo.data[[i]]) <- object@sampling@grid[[1]]  ## to be fixed
     
    #object@model@hurst <-tmphurst
    
    if(missing(subsampling))
      return(object)
    subsampling(object, subsampling)
    
}

aux.simulatPPRWithIntesFeedBack<-function(object, nsim = nsim, seed = seed,
                                xinit = xinit, true.parameter = true.parameter,
                                space.discretized = space.discretized,
                                increment.W = increment.W,
                                increment.L = increment.L, method = method, hurst = hurst,
                                methodfGn = methodfGn, sampling = sampling,
                                subsampling = subsampling,
                                posLambda=posLambda){
  samp <- sampling
  Model <- object@model
  gFun <- object@gFun
  Kern <- object@Kernel
  object@sampling <- samp
  randomGenerator<-object@model@measure$df
  nameIntensityProc <- object@PPR@additional.info
  if(missing(increment.W) | is.null(increment.W)){
    
    if( Model@hurst!=0.5 ){
      
      grid<-sampling2grid(object@sampling)
      isregular<-object@sampling@regular
      
      if((!isregular) || (methodfGn=="Cholesky")){
        dW<-CholeskyfGn(grid, Model@hurst,Model@noise.number)
        yuima.warn("Cholesky method for simulating fGn has been used.")
      } else {
        dW<-WoodChanfGn(grid, Model@hurst,Model@noise.number)
      }
      
    } else {
      
      # delta<-(Terminal-Initial)/n
      delta <- samp@delta
      if(!is.Poisson(Model)){ # if pure CP no need to setup dW
        dW <- rnorm(samp@n * Model@noise.number, 0, sqrt(delta))
        dW <- matrix(dW, ncol=samp@n, nrow=Model@noise.number,byrow=TRUE)
      } else {
        dW <- matrix(0,ncol=samp@n,nrow=1)  # maybe to be fixed
      }
    }
    
  } else {
    dW <- increment.W
  }
  if(missing(xinit)){
    if(length(object@model@xinit)!=0){
      xinit<-numeric(length=length(object@model@xinit))
      for(i in c(1:object@model@equation.number))
        xinit[i] <- eval(object@model@xinit[i])
    }else{
      xinit <- rep(0,object@model@equation.number)
      object@model@xinit<-xinit
    }
  }
  if(missing(hurst)){
    hurst<-0.5
  }
  if(samp@regular){
    tForMeas<-samp@delta
    NumbIncr<-samp@n
    if(missing(true.parameter)){
      eval(parse(text= paste0("measureparam$",
                              object@model@time.variable," <- tForMeas",collapse="")))
    }else{
      measureparam<-true.parameter[object@model@parameter@measure]
      eval(parse(text= paste0("measureparam$",
                              object@model@time.variable," <- tForMeas",collapse="")))
      
    }
    Noise.L <- t(rand(object = randomGenerator, n=NumbIncr, param=measureparam))
    rownames(Noise.L)<-Model@solve.variable
    #dIncr <- apply(cbind(0,Noise.L),1,diff)
    Noise.count <- Noise.L[object@PPR@counting.var,]
    Noise.Laux <- Noise.L
    for(i in c(1:length(object@PPR@counting.var))){
      Noise.Laux[object@PPR@counting.var[i],]<-0
    }  
  }
  myenv<-new.env()
  par.len <- length(object@PPR@allparam)
  if(par.len>0){
    for(i in 1:par.len){
      pars <- object@PPR@allparam[i]
      
      for(j in 1:length(true.parameter)){
        if( is.na(match(pars, names(true.parameter)[j]))!=TRUE){
          assign(object@PPR@allparam[i], true.parameter[[j]],myenv)
        }
      }
    }
  }
  assign("dL",Noise.Laux,myenv)
  
  
  condMatrdW <- is.matrix(dW)
  if(condMatrdW){
    dimdW <- dim(dW)[2]
  }else{
    dimdW <- length(dW)
  }
  assign("info.PPR", object@PPR, myenv)
  
  dimCov <- length(object@PPR@covariates)
  dimNoise<-dim(Noise.Laux)
  
  
  
  
  # Start Simulation PPR
  compErrHazR4 <- function(samp, Kern,
                           capitalTime, Model, 
                           my.env, ExprHaz, Time, dN, 
                           Index, pos){
    assign(Kern@variable.Integral@var.time, Time, envir = my.env)
    
    assign(Model@time.variable, capitalTime, envir = my.env)
    l <- 1
    for(i in c(1:length(Kern@variable.Integral@var.dx)) ){
      if(any(Kern@variable.Integral@var.dx[i]==my.env$info.PPR@counting.var)){
        assign(paste0("d",Kern@variable.Integral@var.dx[i]), dN[l,], envir =my.env)
        l <- l + 1
      }
      if(Kern@variable.Integral@var.dx[i]%in%my.env$info.PPR@covariates){  
        assign(paste0("d",Kern@variable.Integral@var.dx[i]),
               diff(c(0,my.env[[Kern@variable.Integral@var.dx[i]]])) , 
               envir =my.env)
      }
    }
    condPointIngrid <- samp@grid[[1]]<=my.env$t
    PointIngridInt <- samp@grid[[1]][condPointIngrid]
    CondJumpGrid <- PointIngridInt %in% my.env$s
    assign("CondJumpGrid", CondJumpGrid, envir = my.env)
    Lambda <- NULL
    #  for(h in c(1:Index)){
    #      Lambdadum <- eval(ExprHaz[[h]], envir = my.env)
    #      Lambda <- rbind(Lambda,Lambdadum)
    #  
    Lambda <- eval(ExprHaz[[pos]], envir = my.env)
    #    rownames(Lambda) <- my.env$info.PPR@counting.var
    return(Lambda)
    
  }
  if (dimCov>0){
    for(j in c(1:dimCov)){
      assign(object@PPR@covariates[j],
             as.numeric(xinit[j]),
             envir = myenv)
    }
  }
  CovariateSim <-matrix(0,Model@equation.number,(samp@n+1))
  #IntensityProcInter <- matrix(0,length(nameIntensityProc),(samp@n+1))
  ExprHaz <- constHazIntPr(g.Fun = object@gFun@formula,
                           Kern.Fun = object@Kernel, covariates = object@PPR@covariates,
                           counting.var = object@PPR@counting.var, statevar=nameIntensityProc)$Intens
  IntensityProcInter <- as.matrix(tryCatch(eval(object@gFun@formula,envir=myenv),error =function(){1}))
  dN <- matrix(0,object@gFun@dimension[1],object@gFun@dimension[2])
  rownames(CovariateSim)<- Model@solve.variable
  assign(object@PPR@counting.var,CovariateSim[object@PPR@counting.var,1],envir=myenv)
  grid <- samp@grid[[1]]
  const <- -log(runif(gFun@dimension[1]))
  condMyTR <- const<delta
  AllnameIntensityProc <- paste0(nameIntensityProc,"ForKernel")
  assign(AllnameIntensityProc,IntensityProcInter,envir=myenv)
  while(any(condMyTR)){
    if(sum(condMyTR)==0){
      const <- -log(runif(length(condMyTR)))
      condMyTR <- const<delta
    }else{
      const[condMyTR] <- -log(runif(sum(condMyTR)))
      condMyTR <- const<delta
    }
  }
  jumpT<-NULL
  i <- 1
  Initial_i <- i-1
  dimGrid <-length(grid)
  cond <- const
  Index <- gFun@dimension[1]
  inter_i <- rep(i,Index)
  noExit<-rep(T,Index)
  while(any(noExit)){
    
    for(j in c(1:Index)){
      HazardRate<-0
      while(cond[j]>0 && noExit[j]){
        lambda<-compErrHazR4(samp, Kern, capitalTime=samp@grid[[1]][inter_i[j]], 
                             Model, myenv, ExprHaz, Time=jumpT, dN, Index, j)
        if(is.matrix(posLambda)){}else{
          #assign(object@model@state.variable[posLambda],lambda, envir = myenv)
          assign(nameIntensityProc,lambda[j], envir = myenv)
          # myenv[[AllnameIntensityProc]][j,]<-cbind(myenv[[AllnameIntensityProc]][j,],
          #                                                   lambda[j])
          assign(AllnameIntensityProc,
                 cbind(t(myenv[[AllnameIntensityProc]][j,]),
                       lambda[j]),
                 envir=myenv)
        }
        
        
        
        
        incrlambda <- lambda*delta
        HazardRate <- HazardRate+incrlambda
        cond[j] <- const[j]-HazardRate
        # if(cond[j]>0){
        #   dN<-cbind(dN,rep(0,Index))
        # }
        inter_i[j]<-inter_i[j]+1
        if(inter_i[j]-1==1){
          CovariateSim[,c((inter_i[j]-1):inter_i[j])]<- eulerPPRwithInt(xinit=xinit,yuima=object,dW=dW[,(inter_i[j]-1)], 
                                       dL=as.matrix(myenv$dL[,c(i-Initial_i)]),Initial=samp@Initial,Terminal=samp@grid[[1]][inter_i[j]],n=1,
                                       env=myenv)
          rownames(CovariateSim)<- Model@solve.variable
        }else{
          CovariateSim[,inter_i[j]]<- eulerPPRwithInt(xinit=CovariateSim[,(inter_i[j]-1)],
                                                      yuima=object,dW=dW[,(inter_i[j]-1)],
                                                      dL=as.matrix(myenv$dL[,c(inter_i[j]-1-Initial_i)]), 
                                                      Initial=samp@grid[[1]][(inter_i[j]-1)],
                                                      Terminal=samp@grid[[1]][inter_i[j]],n=1,
                                                      env=myenv)[,-1]
        }
        # if(inter_i[j]==66){
        #   aaaaa<-1
        # }
        
        if(inter_i[j]>=(dimGrid)){
          noExit[j] <- FALSE
        }
        if(inter_i[j]<=dimGrid){  
          assign(object@PPR@counting.var,CovariateSim[object@PPR@counting.var[j],1:inter_i[j]],envir=myenv)
          dimCov <- length(object@PPR@covariates)
          if (dimCov>0){
            for(jj in c(1:dimCov)){
              assign(object@PPR@covariates[jj],
                     as.numeric(CovariateSim[object@PPR@covariates[jj],1:inter_i[j]]),
                     envir = myenv)
            }
          }  
        }
      }
    }
    
    i <- min(inter_i)
    Initial_i <- i-1
    if(any(noExit)){
      if(i<dim(CovariateSim)[2]){ 
        jumpT<-c(jumpT,grid[i])
        
        if(dim(dN)[2]==1 & all(dN[,1]==0)){
          dN[i==inter_i,1] <- Noise.count[i-1]
          Noise.Laux[object@PPR@counting.var,i-1]<-Noise.count[i-1]
          dumdN <- dN
        }else{
          dumdN <- rep(0,Index)
          dumdN[i==inter_i] <- Noise.count[i-1] 
          Noise.Laux[object@PPR@counting.var,i-1] <- dumdN[i==inter_i] 
          dN <- cbind(dN,dumdN)
        }
        # cat("\n ", i, grid[i])
        
        # assign("dL",Noise.Laux,myenv)
        # 
        # CovariateSim<- eulerPPR(xinit=xinit,yuima=object,dW=dW,
        #                         Initial=samp@Initial,Terminal=samp@Terminal,
        #                         env=myenv)
        
        assign("dL",Noise.Laux[,c((i-1):dimNoise[2])],myenv)
        xinit <- CovariateSim[,i-1]
        
        
        # if(condMatrdW){
        #   CovariateSim[,(i-1):dimCovariateSim[2]] <- eulerPPRwithInt(xinit=xinit,
        #                                                       yuima=object,dW=dW[,(i-1):dimdW],
        #                                                       Initial=samp@grid[[1]][i-1],Terminal=samp@Terminal,n=(samp@n-(i-1)+1),
        #                                                       env=myenv)
        # }else{
        #   CovariateSim[,(i-1):dimCovariateSim[2]] <- eulerPPRwithInt(xinit=xinit,
        #                                                       yuima=object, dW=dW[(i-1):dimdW],
        #                                                       Initial=samp@grid[[1]][i-1],Terminal=samp@Terminal,n=(samp@n-(i-1)+1),
        #                                                       env=myenv)
        # }
        # 
        # rownames(CovariateSim)<- Model@solve.variable
        
        
        const <- -log(runif(object@gFun@dimension[1]))
        condMyTR <- const<delta
        while(any(condMyTR)){
          if(sum(condMyTR)==0){
            const <- -log(runif(length(condMyTR)))
            condMyTR <- const<delta
          }else{
            const[condMyTR] <- -log(runif(sum(condMyTR)))
            condMyTR <- const<delta
          }
        }
        
        cond <- const
        if(all(noExit)){
          inter_i <- rep(i, Index)
        }else{
          if(any(noExit)){
            inter_i[noExit] <- i
            inter_i[!noExit] <- samp@n+1 
          }
        }
        
      }
    }  
  }
  tsX <- ts(data=t(CovariateSim), deltat=delta, start=object@sampling@Initial)
  object@data <- setData(original.data=tsX)
  for(i in 1:length(object@data@zoo.data))
    index(object@data@zoo.data[[i]]) <- object@sampling@grid[[1]]  ## to be fixed
  
  #object@model@hurst <-tmphurst
  
  if(missing(subsampling))
    return(object)
  subsampling(object, subsampling)
}

aux.simulatPPRROldNew<-function(object, nsim = nsim, seed = seed,
                          xinit = xinit, true.parameter = true.parameter,
                          space.discretized = space.discretized, increment.W = increment.W,
                          increment.L = increment.L, method = method, hurst = 0.5,
                          methodfGn = methodfGn, sampling = sampling,
                          subsampling = subsampling){
  myhawkesP <- function(simMod, Kern,
                        samp, Model, my.env, ExprHaz,
                        Time, dN){
    noExit<-TRUE
    const <- -log(runif(1))
    delta <- samp@delta
    grid <- samp@grid[[1]]
    while(const<delta){
      const <- -log(runif(1))
    }
    jumpT<-NULL
    i <- 1
    dimGrid <-length(grid)
    cond <- const
    allconst <- NULL
    allcond <- NULL
    allhaz <- NULL
    while(noExit){
      HazardRate<-0
      while(cond>0 && noExit){
        #lastJump <- tail(jumpT,n=1L)
        lambda<-compErrHazR2(simMod, Kern, capitalTime=samp@grid[[1]][i], Model, my.env, ExprHaz,
                             Time=jumpT, dN)
        # lambda<-hawkesInt(mu=mu, alpha=alpha, beta=beta,
        #                   timet=grid[i], JumpT=jumpT)
        incrlambda <- lambda*delta
        HazardRate <- HazardRate+incrlambda
        cond <- const-HazardRate
        i<-i+1
        if(i>=(dimGrid-1)){
          noExit <- FALSE
        }
        if(i<dim(simMod@data@original.data)[1]){  
          dimCov <- length(object@PPR@covariates)
          
          if (dimCov>0){
            for(j in c(1:dimCov)){
              # my.covdata <- simMod@data@original.data[1:i,object@PPR@covariates[j]]
              # names(my.covdata) <-simMod@sampling@grid[[1]][1:i]
              # 
              # assign(object@PPR@covariates[j],
              #        my.covdata,
              #        envir = my.env)
              
              assign(object@PPR@covariates[j],
                     as.numeric(simMod@data@original.data[1:i,object@PPR@covariates[j]]),
                     envir = my.env)
            }
          }  
          
          
          # Line 354 necessary for the development of the code.
          # cat("\n ", i, grid[i])
        }
      }
      if(i<dim(simMod@data@original.data)[1]){ 
        jumpT<-c(jumpT,grid[i])
        # if(i==7001){
        #   cat("\n",noExit)
        # }
        if(dN[1]==0){
          #dN <- 1
          dN <- simMod@data@original.data[i,object@PPR@counting.var]-simMod@data@original.data[i-1,object@PPR@counting.var]
        }else{
          dN <- c(dN,
                  simMod@data@original.data[i,object@PPR@counting.var]-simMod@data@original.data[i-1,object@PPR@counting.var])
        }
        #names(dN)<-jumpT
        allhaz <- c(allhaz,HazardRate)
        allcond <- c(allcond,cond)
        cond <- const
        allconst <- c(allconst, const)
        const <- -log(runif(1))
        while(const<delta){
          const <- -log(runif(1))
        }
      }
    }
    return(list(jumpT=jumpT,allcond=allcond,allconst=allconst, allhaz=allhaz))
  }
  
  #myhawkesPMulti
  
  myhawkesPMulti <- function(simMod, Kern,
                        samp, Model, my.env, ExprHaz,
                        Time, dN, Index){
    #noExit<-TRUE
    
    delta <- samp@delta
    grid <- samp@grid[[1]]
    const <- -log(runif(object@gFun@dimension[1]))
    condMyTR <- const<delta
    while(any(condMyTR)){
      if(sum(condMyTR)==0){
        const <- -log(runif(length(condMyTR)))
        condMyTR <- const<delta
      }else{
        const[condMyTR] <- -log(runif(sum(condMyTR)))
        condMyTR <- const<delta
      }
    }
    # while(const<delta){
    #   const <- -log(runif(1))
    # }
    jumpT<-NULL
    i <- 1
    dimGrid <-length(grid)
    cond <- const
    inter_i <- rep(i,Index)
    noExit<-rep(T,Index)
    while(any(noExit)){
        
        for(j in c(1:Index)){
          HazardRate<-0
          while(cond[j]>0 && noExit[j]){
            lambda<-compErrHazR3(simMod, Kern, capitalTime=samp@grid[[1]][inter_i[j]], 
              Model, my.env, ExprHaz, Time=jumpT, dN, Index, j)
        # lambda<-hawkesInt(mu=mu, alpha=alpha, beta=beta,
        #                   timet=grid[i], JumpT=jumpT)
            incrlambda <- lambda*delta
            HazardRate <- HazardRate+incrlambda
            cond[j] <- const[j]-HazardRate
            inter_i[j]<-inter_i[j]+1
        if(inter_i[j]>=(dimGrid-1)){
          noExit[j] <- FALSE
        }
        if(inter_i[j]<dim(simMod@data@original.data)[1]){  
          dimCov <- length(object@PPR@covariates)
          
          if (dimCov>0){
            for(jj in c(1:dimCov)){
              # my.covdata <- simMod@data@original.data[1:i,object@PPR@covariates[j]]
              # names(my.covdata) <-simMod@sampling@grid[[1]][1:i]
              # 
              # assign(object@PPR@covariates[j],
              #        my.covdata,
              #        envir = my.env)
              
              assign(object@PPR@covariates[jj],
                     as.numeric(simMod@data@original.data[1:inter_i[j],object@PPR@covariates[jj]]),
                     envir = my.env)
            }
          }  
          
          
          # Line 354 necessary for the development of the code.
          # cat("\n ", i, grid[i])
        }
      }
        }
      
        i <- min(inter_i)

      if(any(noExit)){
        if(i<dim(simMod@data@original.data)[1]){ 
        jumpT<-c(jumpT,grid[i])

        if(dim(dN)[2]==1 & all(dN[,1]==0)){
          dN[i==inter_i,1] <- 1
          dumdN <- dN
        }else{
          dumdN <- rep(0,Index)
          dumdN[i==inter_i] <- 1 
          dN <- cbind(dN,dumdN)
        }
        #names(dN)<-jumpT
        # const <- -log(runif(1))
        # while(const<delta){
        #   const <- -log(runif(1))
        # }
        const <- -log(runif(object@gFun@dimension[1]))
        condMyTR <- const<delta
        while(any(condMyTR)){
          if(sum(condMyTR)==0){
            const <- -log(runif(length(condMyTR)))
            condMyTR <- const<delta
          }else{
            const[condMyTR] <- -log(runif(sum(condMyTR)))
            condMyTR <- const<delta
          }
        }
        
        cond <- const
        if(all(noExit)){
          inter_i <- rep(i, Index)
        }else{
          if(any(noExit)){
            inter_i[noExit] <- i
            inter_i[!noExit] <- samp@n+1 
          }
        }
          
        }
      }  
    }
    return(list(jumpT=jumpT,dN = dN))
  }
  
  
  # compErrHazR2 <- function(simMod, Kern,
  #                          capitalTime, Model, my.env, ExprHaz,
  #                          Time, dN){
  #   #  dummyLambda <- numeric(length=(samp@n+1))
  #   if(length(Kern@variable.Integral@var.dx)==1){
  #     #   MyPos <- sum(samp@grid[[1]]<=tail(Time,n=1L))
  #     assign(Kern@variable.Integral@var.time, Time, envir = my.env)
  #     #  cond <- -log(cost)-sum(dummyLambda)*samp@delta
  #     
  #     assign(Model@time.variable, capitalTime, envir = my.env)
  #     assign(paste0("d",Kern@variable.Integral@var.dx), dN, envir =my.env)
  #     
  #     condPointIngrid <- simMod@sampling@grid[[1]]<=my.env$t
  #     PointIngridInt <- simMod@sampling@grid[[1]][condPointIngrid]
  #     CondJumpGrid <- PointIngridInt %in% my.env$s
  #     assign("CondJumpGrid", CondJumpGrid, envir = my.env)
  #     
  #     Lambda <- eval(ExprHaz[[1]], envir=my.env)
  #     return(Lambda)
  #   }else{
  #     if(Kern@Integrand@dimIntegrand[1]==1){
  #       assign(Kern@variable.Integral@var.time, Time, envir = my.env)
  #     #  cond <- -log(cost)-sum(dummyLambda)*samp@delta
  #     
  #       assign(Model@time.variable, capitalTime, envir = my.env)
  #       for(i in c(1:length(Kern@variable.Integral@var.dx)) ){
  #         if(Kern@variable.Integral@var.dx[i]==my.env$info.PPR@counting.var){
  #           assign(paste0("d",Kern@variable.Integral@var.dx[i]), dN, envir =my.env)
  #         }
  #         if(Kern@variable.Integral@var.dx[i]%in%my.env$info.PPR@covariates){  
  #           assign(paste0("d",Kern@variable.Integral@var.dx[i]),
  #                  diff(c(0,my.env[[Kern@variable.Integral@var.dx[i]]])) , 
  #                  envir =my.env)
  #         }
  #         if(Kern@variable.Integral@var.dx[i]%in%my.env$info.PPR@var.dt){
  #           assign(paste0("d",Kern@variable.Integral@var.dx[i]),
  #                  diff(c(0,my.env[[Kern@variable.Integral@var.dx[i]]])) , 
  #                  envir =my.env)
  #         }
  #       } 
  #       condPointIngrid <- simMod@sampling@grid[[1]]<=my.env$t
  #       PointIngridInt <- simMod@sampling@grid[[1]][condPointIngrid]
  #       CondJumpGrid <- PointIngridInt %in% my.env$s
  #       assign("CondJumpGrid", CondJumpGrid, envir = my.env)
  #     
  #       Lambda <- eval(ExprHaz[[1]], envir=my.env)
  #     return(Lambda)
  #     }
  #   }
  # }

  compErrHazR2 <- function(simMod, Kern,
                           capitalTime, Model, my.env, ExprHaz,
                           Time, dN){
    #  dummyLambda <- numeric(length=(samp@n+1))
    if(length(Kern@variable.Integral@var.dx)==1){
      #   MyPos <- sum(samp@grid[[1]]<=tail(Time,n=1L))
      assign(Kern@variable.Integral@var.time, Time, envir = my.env)
      #  cond <- -log(cost)-sum(dummyLambda)*samp@delta
      
      assign(Model@time.variable, capitalTime, envir = my.env)
      assign(paste0("d",Kern@variable.Integral@var.dx), dN, envir =my.env)
      
      condPointIngrid <- simMod@sampling@grid[[1]]<=my.env$t
      PointIngridInt <- simMod@sampling@grid[[1]][condPointIngrid]
      CondJumpGrid <- PointIngridInt %in% my.env$s
      assign("CondJumpGrid", CondJumpGrid, envir = my.env)
      
      Lambda <- eval(ExprHaz[[1]], envir=my.env)
      return(Lambda)
    }else{
      if(Kern@Integrand@dimIntegrand[1]==1){
        assign(Kern@variable.Integral@var.time, Time, envir = my.env)
        #  cond <- -log(cost)-sum(dummyLambda)*samp@delta
        
        assign(Model@time.variable, capitalTime, envir = my.env)
        for(i in c(1:length(Kern@variable.Integral@var.dx)) ){
          if(Kern@variable.Integral@var.dx[i]==my.env$info.PPR@counting.var){
            assign(paste0("d",Kern@variable.Integral@var.dx[i]), dN, envir =my.env)
          }
          if(Kern@variable.Integral@var.dx[i]%in%my.env$info.PPR@covariates){  
            assign(paste0("d",Kern@variable.Integral@var.dx[i]),
                   diff(c(0,my.env[[Kern@variable.Integral@var.dx[i]]])) , 
                   envir =my.env)
          }
          if(Kern@variable.Integral@var.dx[i]%in%my.env$info.PPR@var.dt){
            assign(paste0("d",Kern@variable.Integral@var.dx[i]),
                   diff(c(0,my.env[[Kern@variable.Integral@var.dx[i]]])) , 
                   envir =my.env)
          }
        } 
        condPointIngrid <- simMod@sampling@grid[[1]]<=my.env$t
        PointIngridInt <- simMod@sampling@grid[[1]][condPointIngrid]
        CondJumpGrid <- PointIngridInt %in% my.env$s
        assign("CondJumpGrid", CondJumpGrid, envir = my.env)
        
        Lambda <- eval(ExprHaz[[1]], envir=my.env)
        return(Lambda)
      }
    }
  }  
  
  compErrHazR3 <- function(simMod, Kern,
    capitalTime, Model, my.env, ExprHaz, Time, dN, Index, pos){
    assign(Kern@variable.Integral@var.time, Time, envir = my.env)
    
    assign(Model@time.variable, capitalTime, envir = my.env)
    l <- 1
    for(i in c(1:length(Kern@variable.Integral@var.dx)) ){
      if(any(Kern@variable.Integral@var.dx[i]==my.env$info.PPR@counting.var)){
        assign(paste0("d",Kern@variable.Integral@var.dx[i]), dN[l,], envir =my.env)
        l <- l + 1
      }
      if(Kern@variable.Integral@var.dx[i]%in%my.env$info.PPR@covariates){  
        assign(paste0("d",Kern@variable.Integral@var.dx[i]),
               diff(c(0,my.env[[Kern@variable.Integral@var.dx[i]]])) , 
               envir =my.env)
      }
    }
    condPointIngrid <- simMod@sampling@grid[[1]]<=my.env$t
    PointIngridInt <- simMod@sampling@grid[[1]][condPointIngrid]
    CondJumpGrid <- PointIngridInt %in% my.env$s
    assign("CondJumpGrid", CondJumpGrid, envir = my.env)
    Lambda <- NULL
#  for(h in c(1:Index)){
#      Lambdadum <- eval(ExprHaz[[h]], envir = my.env)
#      Lambda <- rbind(Lambda,Lambdadum)
#  
    Lambda <- eval(ExprHaz[[pos]], envir = my.env)
#    rownames(Lambda) <- my.env$info.PPR@counting.var
    return(Lambda)
    
  }
  if(missing(hurst)){
    hurst<-0.5
  }
  samp <- sampling
  Model <- object@model
  gFun <- object@gFun
  Kern <- object@Kernel

  if(missing(xinit)){
    if(object@PPR@RegressWithCount){

      yuima.warn("Counting Variables are also covariates.
                 In this case, the algorthim will be implemented
                 as soon as possible.")
      return(NULL)
    }
  }else{
    if(object@PPR@RegressWithCount){
      yuima.warn("Counting Variables are also covariates.
                 In this case, the algorthim will be implemented
                 as soon as possible.")
      return(NULL)
    }
  }
  if(!object@PPR@RegressWithCount && !object@PPR@IntensWithCount){
    auxg <- setMap(func = gFun@formula, yuima =Model)
    dummyKernIntgrand <- Kern@Integrand@IntegrandList
    dummyUpperTime<- paste0(Kern@variable.Integral@upper.var,
                       Kern@variable.Integral@upper.var,
                       collapse = "")
    dummyTime <-Model@time.variable
    for(i in c(1:length(dummyKernIntgrand))){
      if(Kern@variable.Integral@upper.var %in% all.vars(dummyKernIntgrand[[i]])){
        dumExpr <- paste0("substitute(expression(",
               dummyKernIntgrand[[i]],"), list(",
               Kern@variable.Integral@upper.var,
               " =  as.symbol(dummyUpperTime), ",
               Kern@variable.Integral@var.time,
               " =  as.symbol(Model@time.variable)))")
        dummyKernIntgrand[[i]] <- eval(parse(text=dumExpr))
      }
    }

    auxIntMy <- unlist(lapply(dummyKernIntgrand, FUN = function(X){as.character(X)[2]}))
    auxIntMy <- matrix(auxIntMy, Kern@Integrand@dimIntegrand[1],
      Kern@Integrand@dimIntegrand[2], byrow=T)

    if(object@Kernel@variable.Integral@var.dx==object@Kernel@variable.Integral@var.time){
      auxInt <- setIntegral(yuima = Model,
        integrand = auxIntMy,
        var.dx = Model@time.variable,
        upper.var = dummyUpperTime,
        lower.var = Kern@variable.Integral@lower.var)
    }else{
      auxInt <- setIntegral(yuima = Model,
                            integrand = auxIntMy,
                            var.dx =object@Kernel@variable.Integral@var.dx ,
                            upper.var = dummyUpperTime,
                            lower.var = Kern@variable.Integral@lower.var)
    }
    randomGenerator<-object@model@measure$df
    if(samp@regular){
      tForMeas<-samp@delta
      NumbIncr<-samp@n
      if(missing(true.parameter)){
        eval(parse(text= paste0("measureparam$",
                                object@model@time.variable," <- tForMeas",collapse="")))
      }else{
        measureparam<-true.parameter[object@model@parameter@measure]
        eval(parse(text= paste0("measureparam$",
                                object@model@time.variable," <- tForMeas",collapse="")))

      }
      Noise.L <- t(rand(object = randomGenerator, n=NumbIncr, param=measureparam))
      Noise.W <- t(rnorm(NumbIncr, 0,tForMeas))
      if(length(object@model@diffusion[[1]])>1){
        for(i in c(2:length(object@model@diffusion[[1]]))){
          Noise.W <- rbind(Noise.W, rnorm(NumbIncr, 0,tForMeas))
        }
      }
    if(missing(xinit)){
      simg <- simulate(object = auxg, true.parameter = true.parameter[auxg@Output@param@allparam],
                     sampling = samp, hurst = hurst,
                     increment.W = Noise.W, increment.L = Noise.L)
      simK <- simulate(object = auxInt, true.parameter = true.parameter[auxInt@Integral@param.Integral@allparam],
                       sampling = samp, hurst = hurst,
                       increment.W = Noise.W,
                       increment.L = Noise.L)
      Lambda.data <- simg@data@original.data+simK@data@original.data
      Pos<-0
      globPos<-Pos
      condWhile <- TRUE
      while(condWhile){
        Hazard<--cumsum(as.numeric(Lambda.data)[c(Pos:(samp@n+1))])*samp@delta
        U<-runif(1)
        CondPos <- log(U)<=Hazard
        Pos <- Pos+sum(CondPos)
        if(Pos > (samp@n+1)){
          condWhile <- FALSE
        }else{
          globPos <- c(globPos,Pos)
        }
      }
      globPos <- unique(globPos)
      globPos <- globPos[(globPos<=samp@n)]
      NewNoise.L <- Noise.L
      cod <-Model@solve.variable%in%object@PPR@counting.var
      NeWNoise.W<-Noise.W
      NeWNoise.W[cod,] <- 0
      NewNoise.L[cod,] <- 0
      NewNoise.L[cod,globPos[-1]] <- Noise.L[cod,globPos[-1]]
      simM <- simulate(object = Model, true.parameter = true.parameter[Model@parameter@all],
                       sampling = samp, hurst = hurst,
                       increment.W = NeWNoise.W,
                       increment.L = NewNoise.L)
      object@data <- simM@data
      object@sampling <- samp
      return(object)
      #Lambda.data <- simg@data@original.data+simK@data@original.data
     }else{
       simg <- simulate(object = auxg, xinit=xinit,
                       sampling = samp)
      }
    }
  }else{
    if(!object@PPR@RegressWithCount && object@PPR@IntensWithCount){
      ## Here we consider the case where we have a counting variable in the intensity but
      ## we haven't it in the coefficients of the covariates.

      # Simulation of the noise
      DummyT <- c(true.parameter[Model@parameter@measure], samp@delta)
      names(DummyT) <- c(names(true.parameter[Model@parameter@measure]),
                         Model@time.variable)
      increment.L <- rand(object = Model@measure$df,
               n = samp@n ,
               param = DummyT)
      if(!is.matrix(increment.L)){
        increment.L <- matrix(increment.L,ncol = 1)
      }
      if(missing(xinit)){
        simMod <- simulate(object = Model, hurst = hurst,
          sampling = samp,
          true.parameter = true.parameter[Model@parameter@all],
          increment.L = t(increment.L))
      }else{
        simMod <- simulate(object = Model, hurst = hurst,
          sampling = samp, xinit =xinit,
          true.parameter = true.parameter[Model@parameter@all],
          increment.L = t(increment.L))
      }

      colnames(simMod@data@original.data) <- Model@solve.variable

      Data.tot <- as.matrix(simMod@data@original.data)

      ExprHaz <- constHazIntPr(g.Fun = object@gFun@formula,
        Kern.Fun = object@Kernel, covariates = object@PPR@covariates,
        counting.var = object@PPR@counting.var,
        statevar = object@model@state.variable)$Intens
     # if(FALSE){
     if(length(ExprHaz)>=1){

        Time <- samp@Initial

        my.env <- new.env()
        
        assign("info.PPR", object@PPR, my.env)
        
        for(i in c(1:length(object@PPR@allparam))){
          assign(object@PPR@allparam[i],
            as.numeric(true.parameter[object@PPR@allparam[i]]),
            envir = my.env)
        }

        dimCov <- length(object@PPR@covariates)
        if (dimCov>0){
          for(j in c(1:dimCov)){
              assign(object@PPR@covariates[j],
                   as.numeric(simMod@data@original.data[1,object@PPR@covariates[j]]),
                   envir = my.env)
          }
        }
        
        

        if(object@gFun@dimension[1]==1){
        #if(FALSE){
          IntensityProc <- 0
        #set.seed(1)
          dN <- 0
          prova1 <- myhawkesP(simMod, Kern,
            samp, Model, my.env, ExprHaz, Time, dN)
        }else{
          CPP<-FALSE
          
          IntensityProc <- matrix(0,object@gFun@dimension[1],object@gFun@dimension[2])
          #set.seed(1)
          dN <- matrix(0,object@gFun@dimension[1],object@gFun@dimension[2])
          prova1 <- myhawkesPMulti(simMod, Kern,
            samp, Model, my.env, ExprHaz, Time, dN, 
            Index = object@gFun@dimension[1])
          
          Time<-unique(prova1$jumpT)
          dN <- prova1$dN[,1:length(Time)]
          cond <- samp@grid[[1]][-1] %in% Time
          countVar <- Model@solve.variable %in%  object@PPR@counting.var
          increment.L[!cond, countVar]<-0
          increment.L[cond, countVar]<-t(dN)
          if(missing(xinit)){
            simModNew <- simulate(object = Model, hurst = hurst,
                                  sampling = samp,
                                  true.parameter = true.parameter[Model@parameter@all],
                                  increment.L = t(increment.L))
          }else{
            simModNew <- simulate(object = Model, hurst = hurst,
                                  sampling = samp, xinit =xinit,
                                  true.parameter = true.parameter[Model@parameter@all],
                                  increment.L = t(increment.L))
          }
          object@data<-simModNew@data
          object@sampling<-simModNew@sampling
          
          return(object)
        }
       
        #cond <- samp@grid[[1]][-1] %in% Time[-1]
        Time<-prova1$jumpT  
        cond <- samp@grid[[1]][-1] %in% Time
        countVar <- Model@solve.variable %in%  object@PPR@counting.var
        increment.L[!cond, countVar]<-0
        if(missing(xinit)){
          simModNew <- simulate(object = Model, hurst = hurst,
                             sampling = samp,
                             true.parameter = true.parameter[Model@parameter@all],
                             increment.L = t(increment.L))
        }else{
          simModNew <- simulate(object = Model, hurst = hurst,
                             sampling = samp, xinit =xinit,
                             true.parameter = true.parameter[Model@parameter@all],
                             increment.L = t(increment.L))
        }
        object@data<-simModNew@data
        object@sampling<-simModNew@sampling

        return(object)

      }else{
        my.env <- new.env()
        assign("info.PPR", object@PPR, my.env)
        # ExprHaz
        
        for(i in c(1:length(object@PPR@allparam))){
          assign(object@PPR@allparam[i],
                 as.numeric(true.parameter[object@PPR@allparam[i]]),
                 envir = my.env)
        }
        
        dimCov <- length(object@PPR@covariates)
        if (dimCov>0){
          for(j in c(1:dimCov)){
            assign(object@PPR@covariates[j],
                   as.numeric(simMod@data@original.data[1,object@PPR@covariates[j]]),
                   envir = my.env)
          }
        }
        
        CPP<-FALSE
        
        IntensityProc <- matrix(0,object@gFun@dimension[1],object@gFun@dimension[2])
        #set.seed(1)
        #dN<-matrix(0,object@gFun@dimension[1],object@gFun@dimension[2])
        dN <- NULL
        CondGlobal <- TRUE
        CondWhile <- TRUE
        # JumpTime <- samp@grid[[1]]
        u_bar <- samp@Initial
        u_old <- samp@Initial
        jumpT<-NULL
        posInitial <- 1
        while(CondGlobal){
          CondWhile <- TRUE
          const <- log(runif(object@gFun@dimension[1]))
          delta <- samp@delta
          grid <- samp@grid[[1]]
          condMyTR <- -const<delta
          while(any(condMyTR)){
            if(sum(condMyTR)==0){
              const <- log(runif(length(condMyTR)))
              condMyTR <- -const<delta
            }else{
              const[condMyTR] <- log(runif(sum(condMyTR)))
              condMyTR <- -const<delta
            }
          }
          
          posfin <- sum(samp@grid[[1]]<=u_bar)
          dimGrid <-length(grid)
          cond <- const
          allconst <- NULL
          allcond <- NULL
          allhaz <- NULL
          # if(u_bar>=47.83){
          #   aaa<-1
          # }
          checkFunDum_old <- const
          checkFunDum <- const
          count <- 0
          while(CondWhile & count<20){
            HowManylambda <- posfin-posInitial+1
            lambda <- matrix(NA,object@gFun@dimension[1],HowManylambda)
            for(hh in c(1:HowManylambda)){
              lambda[,hh] <- compErrHazR3(simMod, Kern, 
                                          capitalTime=samp@grid[[1]][hh+(posInitial-1)], 
                                          Model, my.env, ExprHaz,
                                          Time=jumpT, dN, object@gFun@dimension[1])
            }
            # Find the optimal u_i next as minimum
            u_next_comp <- numeric(length = object@gFun@dimension[1])
            FunDum <- numeric(length=object@gFun@dimension[1])
            for(l in c(1:object@gFun@dimension[1])){
              FunDum[l] <- const[l] + sum(lambda[l, ]*delta)
              denomi <- lambda[l,HowManylambda]
              u_next_comp[l] <- u_bar-FunDum[l]/denomi
            }
            u_next <- min(u_next_comp)
            if(abs(tail(grid[grid<=u_next],1L) - tail(grid[grid<=u_bar],1L))<delta/2){
              CondWhile<-FALSE
            }
            condpos <- u_next_comp %in% u_next
            checkFunDumAll <- FunDum[condpos]
            checkFunDum <- checkFunDumAll
            
            if(u_next > u_old){
              if(checkFunDum_old<=0){
                if(checkFunDum <= 0){
                  u_old<-u_bar
                  checkFunDum_old <- checkFunDum
                }else{
                  checkFunDum_old <- checkFunDum
                }
              }
                u_bar <- u_next
            }else{
              if(CondWhile){
                u_bar  <- (u_bar + u_old)/2
              }else{
                u_bar <- u_next
              }  
            }
            
            posfin <- sum(samp@grid[[1]]<=u_bar)
            count <- count+1
          #end while  
          }  
          next_jump <- tail(grid[grid<=u_next],1L)
          dummydN <- rep(0,object@gFun@dimension[1])
          for(hhh in c(1:object@gFun@dimension[1])){
            #condJumpComp <- tail(grid[grid<=u_next_comp[hhh]],1L)==u_bar
            condJumpComp <- u_next == u_next_comp[hhh]
            if(condJumpComp)
              dummydN[hhh] <- 1
          }
          dN <- cbind(dN,dummydN)
          # if(length(jumpT)>0){
          #   if(abs(tail(jumpT,1L)-u_bar)<(delta-10^(-12))){
          #     u_bar <- u_bar + delta
          #   }
          # }
          if(length(jumpT)>0){
            if(tail(jumpT, 1L)+delta >= next_jump){
              next_jump <- next_jump+delta 
            }
          }else{
            if(next_jump < delta){
              next_jump <- next_jump+delta
            }
          }
        
          jumpT<-c(jumpT,next_jump)
         # cat("\n ", c(next_jump, checkFunDum,count))
          
          u_bar <- tail(jumpT,1L)
          posInitial<- sum(grid<=next_jump)
          posfin <- posInitial
          u_old <- next_jump
          if((next_jump+delta)>=samp@Terminal-delta){
            CondGlobal <- FALSE
          }
          
        # end First while
        }
        #return(list(dN=dN,jumpT=jumpT))
        Time<-jumpT  
        cond <- samp@grid[[1]][-1] %in% Time
        countVar <- Model@solve.variable %in%  object@PPR@counting.var
        increment.L[!cond, countVar]<-0
        if(missing(xinit)){
          simModNew <- simulate(object = Model, hurst = hurst,
                                sampling = samp,
                                true.parameter = true.parameter[Model@parameter@all],
                                increment.L = t(increment.L))
        }else{
          simModNew <- simulate(object = Model, hurst = hurst,
                                sampling = samp, xinit =xinit,
                                true.parameter = true.parameter[Model@parameter@all],
                                increment.L = t(increment.L))
        }
        object@data<-simModNew@data
        object@sampling<-simModNew@sampling
        
        return(object)
        
      }
    }
  }
  return(NULL)
}

# SolvePPR <- function(posMid, posLeft, posRight, solveLeft = NULL, solveRight = NULL,
#                      cost, Kern, simMod, samp, Model, ExprHaz,
#                       my.env, Time, IntensityProc){
#
#   if((posMid+1)>=(samp@n+1)){
#     mylist <- list(VeryExit = TRUE)
#     return(mylist)
#   }
#   if((posMid+1)>=(samp@n+1)){
#     mylist <- list(VeryExit = TRUE)
#     return(mylist)
#   }
#
#
#    solveMid<- compErrHazR(posMid, simMod, Kern, samp, Model, my.env, ExprHaz, cost, Time)
#    if(solveMid$solveLambda <= 0){
#      # first check
#      if(solveMid$solveLambda<0 ){
#        if(posLeft == (posMid-1)){
#          if(solveLeft*solveMid$solveLambda<0){
#             mylist <- list()
#             mylist$exit <- TRUE
#             mylist$left <- TRUE
#             mylist$posLeft <- posMid
#             mylist$posRight <- samp@n+1
#             mylist$solveLeft <- solveMid$solveLambda
#             mylist$solveRight <- NULL
#             mylist$Time <- c(Time,samp@grid[[1]][-1][posMid])
#             mylist$IntensityProc <- c(IntensityProc, solveMid$dummyLambda)
#
#            return(mylist)
#          }
#        }
#      solveMidLeft <- compErrHazR(posMid-1, simMod, Kern, samp, Model, my.env, ExprHaz, cost, Time)
#       if(solveMidLeft$solveLambda >=0){
#         mylist <- list()
#         mylist$exit <- TRUE
#         mylist$left <- TRUE
#         mylist$posLeft <- posMid-1
#         mylist$posRight <- samp@n+1
#         mylist$solveLeft <- solveMidLeft$solveLambda
#         mylist$solveRight <- NULL
#         mylist$Time <- c(Time,samp@grid[[1]][-1][posMid-1])
#         mylist$IntensityProc <- c(IntensityProc, solveMidLeft$dummyLambda)
#         return(mylist)
#       }else{
#         mylist <- list()
#         mylist$exit <- FALSE
#         mylist$left <- TRUE
#         mylist$posLeft <- posLeft
#         mylist$posRight <- posMid
#         mylist$solveLeft <- solveLeft
#         mylist$solveRight <-solveMidLeft$solveLambda
#         mylist$Time <- Time
#         mylist$IntensityProc <- c(IntensityProc)
#         return(mylist)
#       }
#      }
#    }
#      if(solveMid$solveLambda==0){
#        mylist <- list()
#        mylist$exit <- TRUE
#        mylist$left <- FALSE
#        mylist$posLeft <-posMid
#        mylist$posRight <- samp@n+1
#        mylist$solveLeft <- solveMid$solveLambda
#        mylist$solveRight <- solveRight
#        mylist$Time <- c(Time,samp@grid[[1]][-1][posMid-1])
#        mylist$IntensityProc <- c(IntensityProc, solveMid$dummyLambda)
#        return(mylist)
#      }
#      if(solveMid$solveLambda > 0 && (posMid+1) <(samp@n+1)){
#        solveMidRight <- compErrHazR(posMid+1, simMod, Kern, samp, Model, my.env, ExprHaz, cost, Time)
#        if(solveMidRight$solveLambda <=0){
#          mylist <- list()
#          mylist$exit <- TRUE
#          mylist$left <- FALSE
#          mylist$posLeft <- posMid+1
#          mylist$posRight <- samp@n+1
#          mylist$solveLeft <-  solveMidRight$solveLambda
#          mylist$solveRight <- solveRight
#          mylist$Time <- c(Time,samp@grid[[1]][-1][posMid+1])
#          mylist$IntensityProc <- c(IntensityProc, solveMidRight$dummyLambda)
#          return(mylist)
#        }else{
#          mylist <- list()
#          mylist$exit <- FALSE
#          mylist$left <- FALSE
#          mylist$posLeft <- posMid+1
#          mylist$posRight <- posRight
#          mylist$solveLeft <- solveMidRight$solveLambda
#          mylist$solveRight <-solveRight
#          mylist$Time <- Time
#          mylist$IntensityProc <- c(IntensityProc)
#          return(mylist)
#        }
#       }
# }


# SolvePPR <- function(TopposInGridIn, OldTimePoint, solveLambdaInOld,
#                      cost, Kern, simMod, samp, Model, ExprHaz, dN,
#                      LastTime, my.env, Time, IntensityProc, checkside = FALSE,
#                      solveLeft=NULL, solveRight=NULL){
#
#   if(is.null(solveLambdaInOld)){
#     solveLambdaOld <- -log(cost)
#     solveLeft <- solveLambdaOld
#     solveRight <- NULL
#     dummyLambda <- numeric(length=(TopposInGridIn+1))
#     if(length(Kern@variable.Integral@var.dx)==1){
#       dN <- rep(0, (TopposInGridIn+1))
#       #if(length(Time)==1){
#         con <- (samp@grid[[1]] %in% Time)
#         dN[c(FALSE, con)[c(1:length(dN))]] <- as.numeric(simMod@data@original.data[c(FALSE, con[-length(con)]),Kern@variable.Integral@var.dx]
#                                          -simMod@data@original.data[con,Kern@variable.Integral@var.dx])
#       #}
#     }else{
#
#     }
#     for(i in c(2:(TopposInGridIn+1))){
#       posInGrid <- i
#       LastTime <- samp@grid[[1]][(posInGrid)]
#       LastStime <- samp@grid[[1]][c(1:(posInGrid-1))]
#       assign(Model@time.variable, LastTime, envir = my.env)
#       assign(Kern@variable.Integral@var.time, LastStime, envir = my.env)
#       assign(paste0("d",Kern@variable.Integral@var.dx), dN[c(2:posInGrid)], envir =my.env)
#       dummyLambda[i] <- eval(ExprHaz[[1]], envir=my.env)
#     }
#     solveLambdaOld00 <- -log(cost)-sum(dummyLambda[c(sum(samp@grid[[1]]<=tail(Time,n=1L)):(TopposInGridIn+1))])
#     if(solveLambdaOld*solveLambdaOld00<0){
#       TotposInGrid<-samp@n
#       mylist <- list(InfTopposInGridInOld = min(TopposInGridIn,TotposInGrid),
#                      supTopposInGridInOld = max(TopposInGridIn,TotposInGrid))
#       if(mylist$InfTopposInGridInOld==TopposInGridIn){
#         mylist$left <- TRUE
#       }else{
#         mylist$left <- FALSE
#       }
#       mylist$TotposInGrid <- TopposInGridIn+1
#       mylist$OldSolveLambda <- solveLambdaOld00
#       mylist$solveLeft <- solveLambdaOld00
#       mylist$solveRight <- solveRight
#       mylist$exit <- TRUE
#       mylist$Time <- c(Time,samp@grid[[1]][-1][mylist$TotposInGrid])
#       mylist$IntensityProc <- c(IntensityProc, tail(dummyLambda,n=1L))
#       return(mylist)
#     }
#
#   }else{
#     if(TopposInGridIn>1){
#       solveLambdaOld <- solveLambdaInOld
#     }else{
#
#       dummyLambda <- numeric(length=(TopposInGridIn-1))
#       if(length(Kern@variable.Integral@var.dx)==1){
#         dN <- rep(0, (TopposInGridIn))
#         dN[(TopposInGridIn)] <- as.numeric(simMod@data@original.data[TopposInGridIn,Kern@variable.Integral@var.dx]
#                                            -simMod@data@original.data[TopposInGridIn-1,Kern@variable.Integral@var.dx])
#       }else{
#
#       }
#       for(i in c(2:(TopposInGridIn))){
#         posInGrid <- i
#         LastTime <- samp@grid[[1]][(posInGrid)]
#         LastStime <- samp@grid[[1]][c(1:(posInGrid-1))]
#         assign(Model@time.variable, LastTime, envir = my.env)
#         assign(Kern@variable.Integral@var.time, LastStime, envir = my.env)
#         assign(paste0("d",Kern@variable.Integral@var.dx), dN[c(1:posInGrid)], envir =my.env)
#         dummyLambda[i] <- eval(ExprHaz[[1]], envir=my.env)
#       }
#
#
#       solveLambdaOld <- -log(cost)-sum(dummyLambda)
#
#     }
#   }
#
#   TotposInGrid <- floor(abs((OldTimePoint)-TopposInGridIn)/2)+min(TopposInGridIn,(OldTimePoint))
#
#   cat(sprintf("\n%.5f ", TotposInGrid))
#
#
#   dummyLambda <- numeric(length=(TotposInGrid-1))
#   if(length(Kern@variable.Integral@var.dx)==1){
#     dN <- rep(0, (TotposInGrid))
#     con <- (samp@grid[[1]] %in% Time)
#     con[TotposInGrid-1] <- TRUE
#     dN[c(FALSE, con)[c(1:length(dN))]] <- as.numeric(simMod@data@original.data[c(FALSE, con[-length(con)]),Kern@variable.Integral@var.dx]
#                                                      -simMod@data@original.data[con,Kern@variable.Integral@var.dx])
#   }else{
#
#   }
#   for(i in c(2:(TotposInGrid))){
#     posInGrid <- i
#     LastTime <- samp@grid[[1]][(posInGrid)]
#     LastStime <- samp@grid[[1]][c(1:(posInGrid-1))]
#     assign(Model@time.variable, LastTime, envir = my.env)
#     assign(Kern@variable.Integral@var.time, LastStime, envir = my.env)
#     assign(paste0("d",Kern@variable.Integral@var.dx), dN[c(2:posInGrid)], envir =my.env)
#     dummyLambda[i] <- eval(ExprHaz[[1]], envir=my.env)
#   }
#
#
#   Solvelambda1 <- -log(cost)-sum(dummyLambda[c(sum(samp@grid[[1]]<=tail(Time,n=1L)):(TotposInGrid))])
#   TotposInGridFin <- TotposInGrid
#
#   if(Solvelambda1*solveLambdaOld < 0 | Solvelambda1*solveLambdaOld > 0){
#
#       if(solveLeft*Solvelambda1>0){
#         #solveLeft<-Solvelambda1
#         TotposInGridFin <- TotposInGridFin+1
#         dummyLambda <- numeric(length=(TotposInGridFin))
#       }else{
#         #solveRight <- Solvelambda1
#         TotposInGridFin <- TotposInGridFin-1
#         dummyLambda <- numeric(length=(TotposInGridFin-1))
#       }
#
#     if(length(Kern@variable.Integral@var.dx)==1){
#       dN <- rep(0, (TotposInGridFin))
#       # dN[(TotposInGridFin)] <- as.numeric(simMod@data@original.data[TotposInGridFin,Kern@variable.Integral@var.dx]
#       #                                     -simMod@data@original.data[TotposInGridFin-1,Kern@variable.Integral@var.dx])
#       con <- (samp@grid[[1]] %in% Time)
#       con[TotposInGridFin-1] <- TRUE
#       dN[c(FALSE, con)[c(1:length(dN))]] <- as.numeric(simMod@data@original.data[c(FALSE, con[-length(con)]),Kern@variable.Integral@var.dx]
#                                                        -simMod@data@original.data[con,Kern@variable.Integral@var.dx])
#     }else{
#
#     }
#     for(i in c(2:(TotposInGridFin))){
#       posInGrid <- i
#       LastTime <- samp@grid[[1]][(posInGrid)]
#       LastStime <- samp@grid[[1]][c(1:(posInGrid-1))]
#       assign(Model@time.variable, LastTime, envir = my.env)
#       assign(Kern@variable.Integral@var.time, LastStime, envir = my.env)
#       assign(paste0("d",Kern@variable.Integral@var.dx), dN[c(2:posInGrid)], envir =my.env)
#       dummyLambda[i] <- eval(ExprHaz[[1]], envir=my.env)
#     }
#
#
#     Solvelambda2 <- -log(cost)-sum(dummyLambda[c(sum(samp@grid[[1]]<=tail(Time,n=1L)):(TotposInGridFin))])
#     if(Solvelambda2*Solvelambda1<0){
#       mylist <- list(InfTopposInGridInOld = min(TopposInGridIn,TotposInGridFin),
#                      supTopposInGridInOld = max(TopposInGridIn,TotposInGridFin))
#       if(mylist$InfTopposInGridInOld==TopposInGridIn){
#         mylist$left <- TRUE
#         mylist$solveLeft <- solveLeft
#         mylist$solveRight <- Solvelambda2
#       }else{
#         mylist$left <- FALSE
#         mylist$solveRight <- solveRight
#         mylist$solveLeft <- Solvelambda2
#       }
#       #  TotposInGrid <- floor(abs(TotposInGrid-TopposInGridIn)/2)+min(TotposInGrid,TopposInGridIn)
#
#       mylist$TotposInGrid <- TotposInGridFin
#       mylist$OldSolveLambda <- Solvelambda2
#
#       mylist$exit <- TRUE
#      # mylist$Time <- c(Time,my.env$t)
#       mylist$Time <- c(Time,samp@grid[[1]][-1][mylist$TotposInGrid])
#       mylist$IntensityProc<- c(IntensityProc,tail(dummyLambda,n=1L))
#       return(mylist)
#     }else{
#       mylist <- list(InfTopposInGridInOld = min(TopposInGridIn,TotposInGridFin),
#                      supTopposInGridInOld = max(TopposInGridIn,TotposInGridFin))
#
#       if(solveLambdaOld>Solvelambda2){
#           mylist$left <- TRUE
#           mylist$solveLeft <- solveLeft
#           mylist$solveRight <- Solvelambda2
#         }else{
#           mylist$left <- FALSE
#           mylist$solveRight <- solveRight
#           mylist$solveLeft <- Solvelambda2
#         }
#       }
#       # if(solveLambdaOld>0){
#       #   if(solveLambdaOld>Solvelambda2){
#       #     mylist$left <- TRUE
#       #   }else{
#       #     TotposInGridFin <- TotposInGridFin-1
#       #     dummyLambda <- numeric(length=(TotposInGridFin-1))
#       #   }
#       # }else{
#       #   if(solveLambdaOld>Solvelambda1){
#       #     TotposInGridFin <- TotposInGridFin+1
#       #     dummyLambda <- numeric(length=(TotposInGridFin))
#       #   }else{
#       #     TotposInGridFin <- TotposInGridFin-1
#       #     dummyLambda <- numeric(length=(TotposInGridFin-1))
#       #   }
#       # }
#
#     #  TotposInGrid <- floor(abs(TotposInGrid-TopposInGridIn)/2)+min(TotposInGrid,TopposInGridIn)
#
#       mylist$TotposInGrid <- TotposInGridFin
#       mylist$OldSolveLambda <- Solvelambda2
#       mylist$exit <- FALSE
#       mylist$Time <- Time
#       mylist$IntensityProc <- IntensityProc
#       return(mylist)
#       #repeat
#     }
#
#   if(Solvelambda1 == 0){
#     mylist <- list(InfTopposInGridInOld = min(TopposInGridIn,TotposInGridFin),
#                    supTopposInGridInOld = max(TopposInGridIn,TotposInGridFin))
#     if(solveLambdaOld>=Solvelambda1){
#       mylist$left <- TRUE
#       mylist$solveLeft <- solveLeft
#       mylist$solveRight <- Solvelambda1
#     }else{
#       mylist$left <- FALSE
#       mylist$solveRight <- solveRight
#       mylist$solveLeft <- Solvelambda1
#     }
#     #  TotposInGrid <- floor(abs(TotposInGrid-TopposInGridIn)/2)+min(TotposInGrid,TopposInGridIn)
#
#     mylist$TotposInGrid <- TotposInGridFin
#     mylist$OldSolveLambda <- Solvelambda2
#     mylist$exit <- TRUE
#   #  mylist$Time <- c(Time,my.env$t)
#     mylist$Time <- c(Time,samp@grid[[1]][-1][mylist$TotposInGrid])
#     mylist$IntensityProc<- c(IntensityProc,tail(dummyLambda,n=1L))
#     return(mylist)
#   }
# }

# compErrHazR <- function(TopposInGrid, simMod, Kern,
#                         samp, Model, my.env, ExprHaz,
#                         cost, Time){
#   dummyLambda <- numeric(length=(TopposInGrid))
#   if(length(Kern@variable.Integral@var.dx)==1){
#     dN <- rep(0, TopposInGrid)
#
#     con <- (samp@grid[[1]] %in% c(Time[-1],samp@grid[[1]][TopposInGrid]))
#     dN[con[c(1:length(dN))]] <- as.numeric(simMod@data@original.data[c(FALSE, con[-length(con)]),Kern@variable.Integral@var.dx]
#                                                      -simMod@data@original.data[con,Kern@variable.Integral@var.dx])
#   }else{}
#   #for(i in c(1:TopposInGrid)){
#   #MyPos
#   MyPos <- sum(samp@grid[[1]]<=tail(Time,n=1L))
#   #dummyLambda <- numeric(length=TopposInGrid)
#   assign(Kern@variable.Integral@var.time, Time, envir = my.env)
#   for(i in c(MyPos:TopposInGrid)){
#     posInGrid <- i
#     LastTime <- samp@grid[[1]][-1][(posInGrid)]
#     #LastStime <- samp@grid[[1]][c(1:posInGrid)]
#     assign(Model@time.variable, LastTime, envir = my.env)
#     #assign(Kern@variable.Integral@var.time, LastStime, envir = my.env)
#     #assign(paste0("d",Kern@variable.Integral@var.dx), dN[c(1:posInGrid)], envir =my.env)
#     assign(paste0("d",Kern@variable.Integral@var.dx), 1, envir =my.env)
#     dummyLambda[i] <- eval(ExprHaz[[1]], envir=my.env)
#   }
#  # solveLambda <- -log(cost)-sum(dummyLambda[c(sum(samp@grid[[1]]<=tail(Time,n=1L)):(TopposInGrid))])*samp@delta
#   solveLambda <- -log(cost)-sum(dummyLambda[c(MyPos:(TopposInGrid))])*samp@delta
#   res <- list(solveLambda = solveLambda, dummyLambda = tail(dummyLambda,n=1L))
#   return(res)
# }




aux.simulatPPRROldVersion <- function(object, nsim = nsim, seed = seed,
                                      xinit = xinit, true.parameter = true.parameter,
                                      space.discretized = space.discretized, increment.W = increment.W,
                                      increment.L = increment.L, method = method, hurst = hurst,
                                      methodfGn = methodfGn, sampling = sampling,
                                      subsampling = subsampling){
  Time <- sampling@Terminal
  numbVardx <- length(object@PPR@var.dx)
  numbCountVar <- length(object@PPR@counting.var)
  U <- runif(numbCountVar)

  my.env<- new.env()

  true.parameter <- unlist(true.parameter)

  if(!all(names(true.parameter)==object@PPR@allparam)){
    yuima.stop("true.parameters mismatch the model parameters")
  }
  for(i in c(1:length(object@PPR@allparam))){
    assign(object@PPR@allparam[i],true.parameter[object@PPR@allparam[i]], envir = my.env)
  }

  assign("t",object@gFun@param@time.var, envir = my.env)


  nameu <- object@gFun@param@time.var
  assign("dt",sampling@delta, envir = my.env)

  if(is.null(increment.W)){
    dimW <- length(object@model@diffusion[[1]])
    W <- matrix(rnorm(dimW*sampling@n,mean=0,sd= sqrt(sampling@delta)),nrow=dimW,ncol=sampling@n)
  }
  Condcovariate <- TRUE
  if(is.null(increment.L)){
    dimL <- length(object@model@jump.coeff[[1]])
    L <- matrix(0,nrow=dimL,ncol=sampling@n)
    Condcovariate <- FALSE
    # if(length(object@PPR@covariates)!=0)
    #    Condcovariate <- TRUE
    cond <- !(object@model@solve.variable %in% object@PPR@counting.var)
    if(any(cond)){
      Condcovariate <- TRUE
    }
    dimMd <- length(object@model@solve.variable)
    dumMod <- setModel(drift = rep("0",dimMd),
                       diffusion = matrix("0",dimMd,1),
                       jump.coeff = diag("1",dimMd,dimMd),
                       measure = object@PPR@Info.measure$measure,
                       measure.type = object@PPR@Info.measure$type,
                       solve.variable = object@model@solve.variable)
    if(length(object@model@parameter@measure)!=0){
      simMod <- simulate(object = dumMod,
                         true.parameter = true.parameter[object@model@parameter@measure],
                         sampling = sampling)
    }else{
      simMod <- simulate(object = dumMod,
                         sampling = sampling)
    }

    L <- t(diff(simMod@data@original.data))
  }

  assign("Condcovariate",Condcovariate, envir = my.env)
  assign("W", W, envir = my.env)


  rownames(L)<- object@model@solve.variable

  assign("L", L, envir = my.env)

  assign("All.labKern",object@Kernel@variable.Integral,envir = my.env)
  assign("All.labgFun",object@gFun@param,envir = my.env)



  Fun1 <- function(u,env){
    part <- seq(0,u,by=env$dt)
    env$t<-part[-length(part)]
    if(Condcovariate){
      yuima<- object@model
      for(i in c(1:length(object@PPR@covariates))){
        assign(object@PPR@covariates[i],
               eval(yuima@xinit[yuima@solve.variable==object@PPR@covariates[i]],
                    envir = env), envir = env)
      }
      if(u!=0){
        # Mat<-matrix(0,length(yuima@solve.variable),length(env$t)+1)
        # for(i in c(1:length(yuima@solve.variable))){
        #   Mat[i,1] = eval(yuima@xinit[i],envir = env)
        # }
        Linc <- env$L[,c(1:(length(part)-1))]
        # Linc[yuima@solve.variable!=object@PPR@covariates,]<-matrix(0,
        #   sum(yuima@solve.variable!=object@PPR@covariates), dim(Linc)[2])
        Linc[yuima@solve.variable!=object@PPR@covariates,] <- 0
        DumUnderlMod <- simulate(yuima, true.parameter = true.parameter,
                                 increment.L = env$L[,c(1:(length(part)-1))],
                                 sampling = setSampling(Terminal = u, n= (length(part)-1)))


        for(i in c(1:length(object@PPR@covariates))){
          VariableDum <- DumUnderlMod@data@original.data[,yuima@solve.variable==object@PPR@covariates[i]]
          assign(object@PPR@covariates[i], as.numeric(VariableDum), envir = env)
        }
      }
    }
    (log(env$U)+sum(eval(env$gFun,envir = env)*env$dt))^2
  }

  Fun2 <- function(u,env){
    u <- max(env$old_u,u)
    dumpart <- seq(0,env$old_u, by=env$dt)
    part <- seq(env$old_u,u,by=env$dt)
    t_k <- env$t
    env$t<-part[-length(part)]
    if(u>=sampling@Terminal){
      # Think a better solution
      my.env$utrue<-u
      return(0)
    }
    if(Condcovariate){
      LevIncr <- env$L[, length(dumpart)+c(1:(length(env$t)))]
      LevIncr[object@PPR@counting.var,]<-0
      yuima<- object@model
      xinit<- numeric(length(object@PPR@covariates))
      names(xinit)<- object@PPR@covariates
      for(i in c(1:length(object@PPR@covariates))){
        xinit[i] <- env[[object@PPR@covariates[i]]]
      }

      xinitCount <- numeric(length(object@PPR@counting.var))
      names(xinitCount) <- object@PPR@counting.var
      for(i in c(1:length(xinitCount))){
        xinitCount[i] <- tail(env[[object@PPR@counting.var[i]]],n = 1)
      }
      xinit <- c(xinit,xinitCount)
      if(part[length(part)]-part[1]!=0){
        DumVarCov  <- simulate(yuima,
                               true.parameter = true.parameter,
                               increment.L = LevIncr,
                               sampling =  setSampling(Terminal = (part[length(part)]-part[1]),
                                                       n = dim(LevIncr)[2]),
                               xinit=xinit[yuima@solve.variable])
        for(i in c(1:length(object@PPR@covariates))){
          VariableDum <- DumVarCov@data@original.data[,yuima@solve.variable==object@PPR@covariates[i]]
          assign(object@PPR@covariates[i], as.numeric(VariableDum), envir = env)
        }
      }else{
        for(i in c(1:length(object@PPR@covariates))){
          VariableDum <- xinit[yuima@solve.variable==object@PPR@covariates[i]]
          assign(object@PPR@covariates[i], as.numeric(VariableDum), envir = env)
        }
      }
      #Insert Here simulation Covariate
    }
    integG <-sum(eval(env$gFun,envir = env)*env$dt)
    env$s <- unique(c(env$s,t_k))[-length(env$s)]
    dumt <- env$t
    num <- length(env$Kern)
    integKer <- 0
    for(j in c(1:length(dumt))){
      env$t <- dumt[j]
      dumKernInt <- 0
      for(i in c(1:num)){
        lab.dx <- env$All.labKern@var.dx[i]
        dumKernInt <- dumKernInt+sum(eval(env$Kern,envir=env)*diff(eval(env[[lab.dx]])))
      }
      integKer <- integKer + dumKernInt
    }
    NewTerm <- 0
    if(env$Condcovariate){
      ## Insert Her
    }
    my.env$utrue<-u
    (log(env$U)+ integG + integKer+NewTerm)^2
  }


  u <- numeric(length = numbCountVar)
  names(u) <- object@PPR@counting.var
  for(i in c(1:numbCountVar)){
    assign("gFun", object@gFun@formula[[i]], envir=my.env)
    assign("U",runif(1),envir = my.env)
    u[i]<- as.numeric(optim(0,Fun1,env=my.env)$par)
  }

  t_1 <- min(u)



  if(t_1>Time){
    yuima.stop("No jump occurs in the considered time interval.
               Increasing Terminal in setSampling is suggested")
  }

  condt1<- u%in%t_1

  namesContVarJump <- names(u[condt1])

  JUMP <- matrix(0,nrow=numbCountVar,ncol=sampling@n)

  rownames(JUMP)<- object@PPR@counting.var
  pos<-sum(sampling@grid[[1]][-1]<=t_1)
  t_1 <- sampling@grid[[1]][-1][pos]
  recordTime<-c(0,t_1)
  pos0<-0

  JUMP[namesContVarJump, pos] <- L[namesContVarJump, pos]
  ntot <- sampling@n
  dL <- L
  dL[object@PPR@counting.var,c((pos0+1):pos)]<-JUMP[object@PPR@counting.var,c((pos0+1):pos)]

  X_mat <- matrix(0, length(object@model@solve.variable),
                  ntot)
  rownames(X_mat) <- object@model@solve.variable

  dummyX <- simulate(object@model, true.parameter = true.parameter,
                     increment.W = if(is.matrix(W[,1:pos])){W[,1:pos]}else{t(as.matrix(W[,1:pos]))},
                     increment.L = if(is.matrix(dL[,1:pos])){dL[,1:pos]}else{t(as.matrix(dL[,1:pos]))},
                     sampling = setSampling(Terminal = t_1,
                                            n = t_1/sampling@delta))
  X_mat[,1:pos] <- t(dummyX@data@original.data)[,-1]

  t_jump <- t_1
  if(length(object@Kernel@variable.Integral@var.dx)==1){
    Comulat.dx <- apply(t(X_mat[object@Kernel@variable.Integral@var.dx,
                                c((pos0+1):pos)]), 1, diff)
  }else{
    Comulat.dx <- apply(t(X_mat[object@Kernel@variable.Integral@var.dx,
                                c((pos0+1):pos)]), 2, diff)
  }



  Index <- matrix(c(1:prod(object@Kernel@Integrand@dimIntegrand)),
                  nrow = object@Kernel@Integrand@dimIntegrand[1],
                  ncol = object@Kernel@Integrand@dimIntegrand[2])

  assign(object@Kernel@variable.Integral@var.time,
         sampling@grid[[1]][c((pos0+1):(pos))],
         envir = my.env)

  assign(object@gFun@param@time.var, t_1, envir = my.env)
  for(i in c(1:object@Kernel@Integrand@dimIntegrand[2])){
    assign(object@Kernel@variable.Integral@var.dx[i],
           as.numeric(Comulat.dx[,i]),
           envir = my.env)
  }
  KernDum <- list()
  for(i in c(1:object@Kernel@Integrand@dimIntegrand[1])){
    dumKern <- expression()
    for(j in c(1:object@Kernel@Integrand@dimIntegrand[2])){
      id <- as.numeric(Index[i,j])
      dumKern <- c(dumKern,object@Kernel@Integrand@IntegrandList[[id]])

    }
    KernDum[[i]] <- dumKern
  }


  udumm <- numeric(length = numbCountVar)
  names(udumm) <- object@Kernel@variable.Integral@var.dx

  assign("L",dL,envir = my.env)
  pos0 <- pos
  assign("pos0", pos, envir = my.env)
  assign("old_u",t_1, envir = my.env)

  while(t_jump<Time){


    oldt_1<-t_1
    for(i in c(1:numbCountVar)){
      assign("gFun", object@gFun@formula[[i]], envir=my.env)
      assign("Kern", KernDum[[i]], envir=my.env)
      my.env$utrue<-0
      while(my.env$utrue<oldt_1){
        assign("U",runif(1),envir = my.env)
        optim((t_1+2*my.env$dt),Fun2,method = "Nelder-Mead",
              env=my.env)$par
        u[i] <- as.numeric(my.env$utrue)
      }
    }

    t_1 <- min(u)

    condt1<- u%in%t_1
    namesContVarJump <- names(u[condt1])

    mypos<-sum(sampling@grid[[1]][-1]<=t_1)
    if((pos0+1)<mypos){
      pos<-sum(sampling@grid[[1]][-1]<=t_1)
      t_jump<- t_1
      t_1 <- sampling@grid[[1]][-1][pos]
      recordTime<-c(recordTime,t_1)


      #if(t_1!=sampling@Terminal){

      pos <- min(pos,dim(L)[2])
      JUMP[namesContVarJump, pos] <- L[namesContVarJump, pos]
      dL[object@PPR@counting.var,c((pos0+1):pos)]<-JUMP[object@PPR@counting.var,c((pos0+1):pos)]
      aa<-setSampling(Terminal = (t_1-my.env$old_u),
                      n = length((pos0+1):pos))
      dummyX <- simulate(object@model, true.parameter = true.parameter,
                         increment.W = if(is.matrix(W[,(pos0+1):pos])){W[,(pos0+1):pos]}else{t(as.matrix(W[,(pos0+1):pos]))},
                         increment.L = if(is.matrix(dL[,(pos0+1):pos])){dL[,(pos0+1):pos]}else{t(as.matrix(dL[,(pos0+1):pos]))},
                         sampling = aa,
                         xinit=X_mat[,(pos0)])
      X_mat[,(pos0+1):pos] <- t(dummyX@data@original.data)[,-1]
      if(length(object@Kernel@variable.Integral@var.dx)==1){
        Comulat.dx <- apply(t(X_mat[object@Kernel@variable.Integral@var.dx,
                                    c((pos0+1):pos)]), 1, diff)
      }else{
        Comulat.dx <- apply(t(X_mat[object@Kernel@variable.Integral@var.dx,
                                    c((pos0+1):pos)]), 2, diff)
      }
      if(!is.matrix(Comulat.dx)){
        Comulat.dx <-t(as.matrix(Comulat.dx))
      }

      Index <- matrix(c(1:prod(object@Kernel@Integrand@dimIntegrand)),
                      nrow = object@Kernel@Integrand@dimIntegrand[1],
                      ncol = object@Kernel@Integrand@dimIntegrand[2])
      assign(object@Kernel@variable.Integral@var.time,
             sampling@grid[[1]][c((pos0+1):(pos))],
             envir = my.env)
      assign(object@gFun@param@time.var, t_1, envir = my.env)
      for(i in c(1:object@Kernel@Integrand@dimIntegrand[2])){

        assign(object@Kernel@variable.Integral@var.dx[i],
               as.numeric(Comulat.dx[,i]),
               envir = my.env)
      }
      pos0<-pos
      assign("pos0", pos, envir = my.env)
      assign("old_u",t_1, envir = my.env)

      #}
    }
    assign("L",dL,envir = my.env)
  }
  X_mat[namesContVarJump,pos]<-X_mat[namesContVarJump,pos]
  res.dum <- list(X_mat=X_mat,timeJump = recordTime, grid=sampling)

  solve.variable <-unique(c(object@model@solve.variable))
  N.VarPPR<-length(solve.variable)

  dummy.mod <- setModel(drift=rep("0",N.VarPPR),
                        diffusion = NULL, jump.coeff = diag(rep("1",N.VarPPR)),
                        measure = object@PPR@Info.measure$measure,
                        measure.type = object@PPR@Info.measure$type,
                        solve.variable = solve.variable, xinit=c(object@model@xinit))

  mynewincr <- if(is.matrix(res.dum$X_mat)){t(as.matrix(apply(cbind(0,res.dum$X_mat),1,diff)))}else{apply(cbind(0,res.dum$X_mat),1,diff)}

  interResMod <- simulate(object = dummy.mod,
                          true.parameter = true.parameter,
                          sampling = sampling,
                          increment.L = mynewincr)

  resGfun<-new("yuima.Map",
               Output = object@gFun,
               yuima=setYuima(model=dummy.mod,sampling = sampling))

  interResGfun <- simulate(object = resGfun,
                           true.parameter = true.parameter,
                           sampling = sampling,
                           increment.L = mynewincr)
  dummyObject <- object@Kernel
  dummyObject@variable.Integral@out.var <-object@PPR@additional.info
  resInt <- new("yuima.Integral",
                Integral = dummyObject,
                yuima = setYuima(model=dummy.mod,sampling = sampling))

  interResInt <- simulate(object = resInt,
                          true.parameter = true.parameter,
                          sampling = sampling,
                          increment.L = mynewincr)
  DataIntensity <- interResGfun@data@original.data + interResInt@data@original.data
  InterMDia<-zoo(interResMod@data@original.data, order.by = index(DataIntensity))
  Alldata <-merge(InterMDia,DataIntensity)
  colnames(Alldata)<-c(solve.variable,object@PPR@additional.info)
  # for(i in c(1:N.VarPPR)){
  #   assign(solve.variable[i],interRes@data@original.data[,i],envir=my.env)
  # }
  # dummy<-NULL
  # for(t in c(1:length(object@PPR@additional.info))){
  #   dummy <-eval(object@gFun)
  #   assign(object@PPR@additional.info[[]])
  # }
  object@data<-setData(Alldata)
  return(object)
}


# simOzaki.aux<-function(gFun,a,cCoeff, Time, numJump){
#   t_k<-0
#   N<-0
#   S<-1
#
#   T_k<-c(t_k)
#
#   N_k<-c(N)
#   U<-runif(1)
#   t_k <- -log(U)/gFun
#   if(t_k<Time){
#     T_k<-c(T_k,t_k)
#     N<-N+numJump
#     N_k<-c(N_k, N)
#   }
#   while(t_k<=Time){
#     U<-runif(1)
#     optim.env<-new.env()
#     assign("U",U,envir=optim.env)
#     assign("t_k",t_k,envir=optim.env)
#     assign("c",cCoeff,envir=optim.env)
#     assign("a",a,envir=optim.env)
#     assign("S",S,envir=optim.env)
#     assign("gFun",gFun,envir=optim.env)
#
#     min<-function(u,env){
#       U<-env$U
#       t_k<-env$t_k
#       c<-env$c
#       a<-env$a
#       S<-env$S
#       gFun<-env$gFun
#       y<-(log(U)+gFun*(u-t_k)+c/a*S*(1-exp(-a*(u-t_k))))^2
#     }
#     y<-optim(par=t_k,min, env=optim.env )$par
#     S<- exp(-a*(y-t_k))*S+1
#     t_k<-y
#     T_k<-c(T_k,t_k)
#     N<-N+numJump
#     N_k<-c(N_k, N)
#   }
#   return(list(T_k=T_k,N_k=N_k))
# }
