## Here we write all auxiliar functions for the Point Process
## Regression Model
is.Ppr <- function(yuimaPpr){is(yuimaPpr,"yuima.Ppr")}

quasiLogLik.Ppr <- function(yuimaPpr, parLambda=list(), method=method, fixed = list(),
                lower, upper, call, ...){
  PprData<-yuimaPpr@data
  Time <- index(yuimaPpr@data@zoo.data[[1]])
  envPpr <- list()
  dY <- paste0("d",yuimaPpr@Ppr@var.dx)
  for(i in (c(1:(length(Time)-1)))){
    envPpr[[i]]<-new.env()
    assign(yuimaPpr@gFun@param@time.var,rep(Time[i+1],i),envir=envPpr[[i]])
    assign(yuimaPpr@Ppr@var.dt,Time[1:i],envir=envPpr[[i]])
    if(length(yuimaPpr@Ppr@covariates)>0){
      for(j in c(1:length(yuimaPpr@Ppr@covariates))){
        cond<-colnames(yuimaPpr@data@original.data)%in%yuimaPpr@Ppr@covariates[[j]]
        assign(yuimaPpr@Ppr@covariates[[j]],
               as.numeric(yuimaPpr@data@original.data[1:(i+1),cond]),
               envir=envPpr[[i]])
      }
    }
    for(j in c(1:length(yuimaPpr@Ppr@counting.var))){
      cond<-colnames(yuimaPpr@data@original.data)%in%yuimaPpr@Ppr@counting.var[[j]]
      assign(yuimaPpr@Ppr@counting.var[[j]],
             as.numeric(yuimaPpr@data@original.data[1:(i+1),cond]),
             envir=envPpr[[i]])
    }

    for(j in c(1:length(yuimaPpr@Ppr@var.dx))){
      cond<-c(colnames(yuimaPpr@data@original.data),yuimaPpr@Ppr@var.dt)%in%c(yuimaPpr@Ppr@var.dx,yuimaPpr@Ppr@var.dt)[[j]]
      if(any(cond[-length(cond)])){
        assign(paste0("d",yuimaPpr@Ppr@var.dx[[j]]),
               diff(as.numeric(yuimaPpr@data@original.data[1:(i+1),cond[-length(cond)]])),
               envir=envPpr[[i]])
      }
      if(tail(cond,n=1L)){
        assign(paste0("d",yuimaPpr@Ppr@var.dx[[j]]),
               as.numeric(diff(Time[1:(i+1),cond[-length(cond)]])),
               envir=envPpr[[i]])
      }
    }
  }
  IntKernExpr<- function(Kern,dy){
    dum<-paste(Kern,dy,sep="*")
    dum<-paste0(dum, collapse = " + ")
    dum <- paste0("sum( ", dum, " )")
    return(parse(text = dum))
  }

  IntegKern <- lapply(yuimaPpr@Kernel@Integrand@IntegrandList,IntKernExpr,dY)
  Integrator <- t(as.matrix(eval(parse(text=dY[1]),envir=envPpr[[length(envPpr)]])))
  if(length(dY)>1){
    for(i in c(2:length(dY))){
      Integrator <- rbind(Integrator,
                          t(as.matrix(eval(parse(text=dY[1]),envir=envPpr[[length(envPpr)]]))))
    }
  }
  assign("Integrator",Integrator,envir=envPpr[[length(envPpr)]])
  assign("Nlamb",length(yuimaPpr@Ppr@counting.var),envir=envPpr[[length(envPpr)]])

  out<-NULL
  param1<-unlist(parLambda)
  my.env <- envPpr
  if(length(lower)==0 && length(upper)>0 && length(fixed)==0){
    out <- optim(par=param1, fn=aux.lambdaFromData,
      method = method, upper=upper, envPpr = my.env, gFun=yuimaPpr@gFun,
      Kern =IntegKern, intensityParm = yuimaPpr@Ppr@allparamPpr,
      logLikelihood=TRUE, ...)

  }

  if(length(lower)==0 && length(upper)==0 && length(fixed)>0){
    out <- optim(par = param1, fn = aux.lambdaFromData,
                 method = method, fixed = fixed, envPpr = my.env,
                 gFun = yuimaPpr@gFun, Kern = IntegKern,
                 intensityParm = yuimaPpr@Ppr@allparamPpr,
                 logLikelihood = TRUE, ...)

  }


  if(length(lower)>0 && length(upper)==0 && length(fixed)==0){
    out <- optim(par = param1, fn = aux.lambdaFromData,
                 method = method, lower=lower, envPpr = my.env,
                 gFun = yuimaPpr@gFun, Kern = IntegKern,
                 intensityParm = yuimaPpr@Ppr@allparamPpr,
                 logLikelihood = TRUE, ...)
  }

  if(length(lower)>0 && length(upper)>0 && length(fixed)==0){
    out <- optim(par = param1, fn = aux.lambdaFromData,
                 method = method, upper = upper,
                 lower=lower, envPpr = my.env,
                 gFun = yuimaPpr@gFun, Kern = IntegKern,
                 intensityParm = yuimaPpr@Ppr@allparamPpr,
                 logLikelihood = TRUE, ...)
  }


  if(length(lower)==0 && length(upper)>0 && length(fixed)>0){
    out <- optim(par = param1, fn = aux.lambdaFromData,
                 method = method, upper = upper,
                 fixed = fixed, envPpr = my.env,
                 gFun = yuimaPpr@gFun, Kern = IntegKern,
                 intensityParm = yuimaPpr@Ppr@allparamPpr,
                 logLikelihood = TRUE, ...)
  }

  if(length(lower)>0 && length(upper)==0 && length(fixed)>0){
    out <- optim(par=param1, fn=aux.lambdaFromData,
                 method = method, lower = lower,
                 fixed = fixed, envPpr = my.env,
                 gFun = yuimaPpr@gFun, Kern = IntegKern,
                 intensityParm = yuimaPpr@Ppr@allparamPpr,
                 logLikelihood = TRUE, ...)
  }


  if(length(lower)>0 && length(upper)>0 && length(fixed)>0){
    out <- optim(par = param1, fn = aux.lambdaFromData,
      method = method, lower = lower, fixed = fixed, upper = upper,
      envPpr = my.env, gFun = yuimaPpr@gFun, Kern = IntegKern,
      intensityParm = yuimaPpr@Ppr@allparamPpr, logLikelihood = TRUE, ...)
  }


  if(is.null(out)){
    out <- optim(par = param1, fn = aux.lambdaFromData,
                 method = method, envPpr = my.env, gFun = yuimaPpr@gFun,
                 Kern = IntegKern, intensityParm = yuimaPpr@Ppr@allparamPpr,
                 logLikelihood = TRUE, ...)
  }


  return(out)

}

lambdaFromData <- function(yuimaPpr, PprData=NULL, parLambda=list()){
 if(is.null(PprData)){
   PprData<-yuimaPpr@data
 }else{
  # checklambdaFromData(yuimaPpr,PprData)
 }
 if(!any(names(parLambda) %in% yuimaPpr@Ppr@allparamPpr)){yuima.stop("1 ...")}
 if(!any(yuimaPpr@Ppr@allparamPpr %in% names(parLambda))){yuima.stop("2 ...")}
  Time <- index(yuimaPpr@data@zoo.data[[1]])
  envPpr <- list()
  dY <- paste0("d",yuimaPpr@Ppr@var.dx)
  for(i in (c(1:(length(Time)-1)))){
   envPpr[[i]]<-new.env()
    assign(yuimaPpr@gFun@param@time.var,rep(Time[i+1],i),envir=envPpr[[i]])
    assign(yuimaPpr@Ppr@var.dt,Time[1:i],envir=envPpr[[i]])
    if(length(yuimaPpr@Ppr@covariates)>0){
      for(j in c(1:length(yuimaPpr@Ppr@covariates))){
        cond<-colnames(yuimaPpr@data@original.data)%in%yuimaPpr@Ppr@covariates[[j]]
          assign(yuimaPpr@Ppr@covariates[[j]],
                as.numeric(yuimaPpr@data@original.data[1:(i+1),cond]),
                envir=envPpr[[i]])
      }
    }
    for(j in c(1:length(yuimaPpr@Ppr@counting.var))){
      cond<-colnames(yuimaPpr@data@original.data)%in%yuimaPpr@Ppr@counting.var[[j]]
      assign(yuimaPpr@Ppr@counting.var[[j]],
             as.numeric(yuimaPpr@data@original.data[1:(i+1),cond]),
             envir=envPpr[[i]])
    }

    for(j in c(1:length(yuimaPpr@Ppr@var.dx))){
      cond<-c(colnames(yuimaPpr@data@original.data),yuimaPpr@Ppr@var.dt)%in%c(yuimaPpr@Ppr@var.dx,yuimaPpr@Ppr@var.dt)[[j]]
      if(any(cond[-length(cond)])){
      assign(paste0("d",yuimaPpr@Ppr@var.dx[[j]]),
             diff(as.numeric(yuimaPpr@data@original.data[1:(i+1),cond[-length(cond)]])),
             envir=envPpr[[i]])
      }
      if(tail(cond,n=1L)){
        assign(paste0("d",yuimaPpr@Ppr@var.dx[[j]]),
               as.numeric(diff(Time[1:(i+1),cond[-length(cond)]])),
               envir=envPpr[[i]])
      }
    }
  }
  IntKernExpr<- function(Kern,dy){
    dum<-paste(Kern,dy,sep="*")
    dum<-paste0(dum, collapse = " + ")
    dum <- paste0("sum( ", dum, " )")
    return(parse(text = dum))
  }

  IntegKern <- lapply(yuimaPpr@Kernel@Integrand@IntegrandList,IntKernExpr,dY)
  Integrator <- t(as.matrix(eval(parse(text=dY[1]),envir=envPpr[[length(envPpr)]])))
  if(length(dY)>1){
    for(i in c(2:length(dY))){
      Integrator <- rbind(Integrator,
                          t(as.matrix(eval(parse(text=dY[1]),envir=envPpr[[length(envPpr)]]))))
    }
  }
  assign("Integrator",Integrator,envir=envPpr[[length(envPpr)]])
  assign("Nlamb",length(yuimaPpr@Ppr@counting.var),envir=envPpr[[length(envPpr)]])
  res<-aux.lambdaFromData(param = unlist(parLambda), gFun=yuimaPpr@gFun,
    Kern =IntegKern, intensityParm = yuimaPpr@Ppr@allparamPpr,
    envPpr)
  return(res)
}
# my.lapply <- function (X, FUN, ...){
# #  FUN <- match.fun(FUN)
#   .Internal(lapply(X, FUN))
# }
myfun3<-function(X,Kern){
  t(simplify2array(lapply(Kern, FUN=DumFun,
                             Y = X), higher = (TRUE == "array")))}
dumFun2<-function(X,Y){list2env(Y,envir=X)}

aux.lambdaFromData <-function(param, gFun, Kern, intensityParm, envPpr,logLikelihood = FALSE){
  lapply(envPpr,FUN=dumFun2,Y=as.list(param))
  lastEnv <- tail(envPpr,n=1L)[[1]]
  # gFunVect<- t(simplify2array(my.lapply(gFun@formula, FUN=DumFun,
  #            Y = lastEnv), higher = (TRUE == "array")))

  gFunVect<- matrix(unlist(lapply(gFun@formula, FUN=DumFun,
     Y = lastEnv)),nrow=lastEnv$Nlamb,byrow=TRUE)
  # IntKer<- simplify2array(my.lapply(envPpr,function(my.env){
  #   t(simplify2array(my.lapply(Kern, FUN=DumFun,
  #             Y = my.env), higher = (TRUE == "array")))}),
  #   higher = (TRUE == "array")
  #   )
  # IntKer<- simplify2array(my.lapply(envPpr,myfun3,Kern=Kern),
  #                         higher = (TRUE == "array")
  # )

  # IntKer<- matrix(unlist(my.lapply(envPpr,myfun3,Kern=Kern)),
  #        nrow=1,byrow=TRUE)
  IntKer<- matrix(unlist(lapply(envPpr,myfun3,Kern=Kern)),
              nrow=lastEnv$Nlamb)
  lambda <- gFunVect+cbind(0,IntKer)
  time <- (c(lastEnv$s,lastEnv$t[1]))
  if(!logLikelihood){
    Intensity <- zoo(t(lambda), order.by = time)
    return(Intensity)
  }
  dn <- dim(lambda)
  if(dn[1]==1){
    logLiklihood2 <- -sum(lambda[,-1]*diff(time)[1])
    logLiklihood1 <- sum(log(lambda[,-1])*lastEnv$Integrator)
  }else{
    logLiklihood2 <- -rowSums(lambda[,-1]*diff(time)[1])
    logLiklihood1 <- rowSums(log(lambda[,-1])*lastEnv$Integrator)
  }
  minusLoglik <- -sum(logLiklihood2+logLiklihood1)
  #cat(sprintf("\n%.5f",minusLoglik))
  return(minusLoglik)
}

DumFun<- function(X,Y){eval(X,envir=Y)}
