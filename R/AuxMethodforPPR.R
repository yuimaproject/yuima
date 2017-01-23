## Here we write all auxiliar functions for the Point Process
## Regression Model
is.Ppr <- function(yuimaPpr){is(yuimaPpr,"yuima.Ppr")}

Internal.LogLikPPR <- function(param,my.envd1=NULL,
                               my.envd2=NULL,my.envd3=NULL){
  param<-unlist(param)

  IntLambda<-InternalConstractionIntensity2(param,my.envd1,
                                           my.envd2,my.envd3)

  # IntLambda<-InternalConstractionIntensity(param,my.envd1,
  #                                          my.envd2,my.envd3)
  Index<-my.envd3$gridTime
  Integr1 <- -sum(IntLambda[-length(IntLambda)]*diff(Index),na.rm=TRUE)
  # if(is.nan(Integr1)){
  #   Integr1 <- -10^6
  # }
  cond2 <- diff(as.numeric(my.envd3$YUIMA.PPR@data@original.data))
  Integr2<- sum(log(IntLambda[-1][cond2>0]),na.rm=TRUE)
  # if(is.nan(Integr2)){
  #   Integr2 <- -10^6
  # }
  logLik <- Integr1+Integr2
  cat("\n ",logLik, param)
  return(-logLik)
}


quasiLogLik.Ppr <- function(yuimaPpr, parLambda=list(), method=method, fixed = list(),
                            lower, upper, call, ...){

  yuimaPpr->yuimaPPr
  parLambda->param
  gfun<-yuimaPPr@gFun@formula

  gfun<-yuimaPPr@gFun@formula

  dimIntegr <- length(yuimaPPr@Kernel@Integrand@IntegrandList)
  Integrand2 <- character(length=dimIntegr)
  for(i in c(1:dimIntegr)){
    Integrand1 <- as.character(yuimaPPr@Kernel@Integrand@IntegrandList[[i]])
    timeCond <- paste0(" * (",yuimaPPr@Kernel@variable.Integral@var.time," < ",yuimaPPr@Kernel@variable.Integral@upper.var,")")
    Integrand2[i] <-paste0(Integrand1,timeCond)
  }

  Integrand2<- matrix(Integrand2,yuimaPPr@Kernel@Integrand@dimIntegrand[1],yuimaPPr@Kernel@Integrand@dimIntegrand[2])


  for(j in c(1:yuimaPPr@Kernel@Integrand@dimIntegrand[2])){
    Integrand2[,j]<-paste0(Integrand2[,j]," * d",yuimaPPr@Kernel@variable.Integral@var.dx[j])
  }
  colnames(Integrand2) <- paste0("d",yuimaPPr@Kernel@variable.Integral@var.dx)
  NamesIntegrandExpr <- as.character(matrix(colnames(Integrand2), dim(Integrand2)[1],dim(Integrand2)[2], byrow = TRUE))
  Integrand2expr<- parse(text=Integrand2)

  gridTime <- time(yuimaPPr@data@original.data)

  yuimaPPr@Kernel@variable.Integral@var.dx
  if(any(yuimaPPr@Kernel@variable.Integral@var.dx %in% yuimaPPr@model@solve.variable)){
    my.envd1<-new.env()
    ExistdN<-TRUE
  }else{
    ExistdN<-FALSE
  }
  Univariate<-FALSE
  if(length(yuimaPPr@Ppr@counting.var)==1){
    Univariate<-TRUE
  }
  if(any(!(yuimaPPr@Kernel@variable.Integral@var.dx %in% yuimaPPr@model@solve.variable))){
    my.envd2<-new.env()
    ExistdX<-TRUE
  }else{
    my.envd2<-new.env()
    ExistdX<-FALSE
  }

  my.envd3 <- new.env()
  namesparam<-names(param)
  if(!(all(namesparam %in% yuimaPPr@Ppr@allparamPpr) && length(namesparam)==length(yuimaPPr@Ppr@allparamPpr))){
    return(NULL)
  }

  # construction my.envd1
  if(ExistdN){

    #CountingVariable
    for(i in c(1:length(yuimaPPr@Ppr@counting.var))){
      cond <- yuimaPPr@Ppr@counting.var[i] %in% yuimaPPr@model@solve.variable
      dummyData <-unique(yuimaPPr@data@original.data[,cond])[-1]
      assign(yuimaPPr@Ppr@counting.var[i], dummyData,envir=my.envd1)
    }
    # Names expression
    assign("NamesIntgra", NamesIntegrandExpr, envir=my.envd1)
    #dN
    namedX <-NULL
    namedJumpTimeX <- NULL
    for(i in c(1:length(yuimaPPr@Kernel@variable.Integral@var.dx))){
      if(yuimaPPr@Kernel@variable.Integral@var.dx[i] %in% yuimaPPr@Ppr@counting.var){
        cond <- yuimaPPr@model@solve.variable %in% yuimaPPr@Kernel@variable.Integral@var.dx[i]
        namedX<-c(namedX,paste0("d",yuimaPPr@Kernel@variable.Integral@var.dx[i]))
        namedJumpTimeX <-c(namedJumpTimeX,paste0("JumpTime.d",yuimaPPr@Kernel@variable.Integral@var.dx[i]))
        dummyData <- diff(as.numeric(yuimaPPr@data@original.data[,cond]))# We consider only Jump
        dummyJumpTime <- gridTime[-1][dummyData>0]
        dummyData2 <- diff(unique(cumsum(dummyData)))
        #dummyData3 <- zoo(dummyData2,order.by = dummyJumpTime)
        dummyData3 <- dummyData2
        JumpTime <- dummyJumpTime
        assign(paste0("d",yuimaPPr@Kernel@variable.Integral@var.dx[i]), dummyData3 ,envir=my.envd1)
        assign(paste0("JumpTime.d",yuimaPPr@Kernel@variable.Integral@var.dx[i]), dummyJumpTime ,envir=my.envd1)
      }
    }
    assign("namedX",namedX, envir = my.envd1)
    assign("namedJumpTimeX",namedJumpTimeX, envir = my.envd1)
    assign("var.time",yuimaPPr@Kernel@variable.Integral@var.time,envir=my.envd1)
    assign("t.time",yuimaPPr@Kernel@variable.Integral@upper.var,envir=my.envd1)

    # Covariates
    if(length(yuimaPPr@Ppr@covariates)>1){
      # Covariates should be identified at jump time
      return(NULL)
    }

  }
  # end coonstruction my.envd1

  # construction my.envd2
  if(ExistdX){
    #Covariate

    #CountingVariable
    for(i in c(1:length(yuimaPPr@Ppr@counting.var))){
      cond <- yuimaPPr@Ppr@counting.var[i] %in% yuimaPPr@model@solve.variable
      dummyData <-yuimaPPr@data@original.data[,cond]
      assign(yuimaPPr@Ppr@counting.var[i], dummyData,envir=my.envd1)
    }


  }else{
    assign("KerneldX",NULL,envir=my.envd2)
  }

  # end construction my.envd2

  # construction my.envd3

  #Covariate

  #CountingVariable
  for(i in c(1:length(yuimaPPr@Ppr@counting.var))){
    cond <- yuimaPPr@Ppr@counting.var[i] %in% yuimaPPr@model@solve.variable
    dummyData <-yuimaPPr@data@original.data[,cond]
    assign(yuimaPPr@Ppr@counting.var[i], dummyData,envir=my.envd3)
  }
  #time
  assign(yuimaPPr@model@time.variable, gridTime, my.envd3)

  #Model
  assign("YUIMA.PPR",yuimaPPr,envir=my.envd3)
  assign("namesparam",namesparam,envir=my.envd3)
  assign("gfun",gfun,envir=my.envd3)
  assign("Integrand2",Integrand2,envir=my.envd3)
  assign("Integrand2expr",Integrand2expr,envir=my.envd3)

  assign("gridTime",as.numeric(gridTime),envir=my.envd3)
  assign("Univariate",Univariate,envir=my.envd3)
  assign("ExistdN",ExistdN,envir=my.envd3)
  assign("ExistdX",ExistdX,envir=my.envd3)



  out<-NULL

  if(length(lower)==0 && length(upper)>0 && length(fixed)==0){
    out <- optim(par=param, fn=Internal.LogLikPPR,
                 my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
                 method = method, upper=upper, ...)

  }

  if(length(lower)==0 && length(upper)==0 && length(fixed)>0){
    out <- optim(par=param, fn=Internal.LogLikPPR,
                 my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
                 method = method, fixed = fixed, ...)

  }


  if(length(lower)>0 && length(upper)==0 && length(fixed)==0){
    out <- optim(par = param,  fn=Internal.LogLikPPR,
                 my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
                 method = method, lower=lower, ...)
  }

  if(length(lower)>0 && length(upper)>0 && length(fixed)==0){
    out <- optim(par=param, fn=Internal.LogLikPPR,
                 my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
                 method = method, upper = upper,
                 lower=lower, ...)
  }


  if(length(lower)==0 && length(upper)>0 && length(fixed)>0){
    out <- optim(par=param, fn=Internal.LogLikPPR,
                 my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
                 method = method, upper = upper,
                 fixed = fixed,  ...)
  }

  if(length(lower)>0 && length(upper)==0 && length(fixed)>0){
    out <- optim(par=param, fn=Internal.LogLikPPR,
                 my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
                 method = method, lower = lower,
                 fixed = fixed, ...)
  }


  if(length(lower)>0 && length(upper)>0 && length(fixed)>0){
    out <- optim(par=param, fn=Internal.LogLikPPR,
                 my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
                 method = method, lower = lower, fixed = fixed, upper = upper, ...)
  }


  if(is.null(out)){
    out <- optim(par=param, fn=Internal.LogLikPPR,
                 my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
                 method = method, ...)
  }


  return(out)

}


# quasiLogLik.Ppr <- function(yuimaPpr, parLambda=list(), method=method, fixed = list(),
#                             lower, upper, call, ...){
#
#   yuimaPpr->yuimaPPr
#   parLambda->param
#   gfun<-yuimaPPr@gFun@formula
#
#   dimIntegr <- length(yuimaPPr@Kernel@Integrand@IntegrandList)
#   Integrand2 <- character(length=dimIntegr)
#   for(i in c(1:dimIntegr)){
#     Integrand1 <- as.character(yuimaPPr@Kernel@Integrand@IntegrandList[[i]])
#     timeCond <- paste0(" * (",yuimaPPr@Kernel@variable.Integral@var.time," < ",yuimaPPr@Kernel@variable.Integral@upper.var,")")
#     Integrand2[i] <-paste0(Integrand1,timeCond)
#   }
#
#   Integrand2<- matrix(Integrand2,yuimaPPr@Kernel@Integrand@dimIntegrand[1],yuimaPPr@Kernel@Integrand@dimIntegrand[2])
#
#
#   for(j in c(1:yuimaPPr@Kernel@Integrand@dimIntegrand[2])){
#     Integrand2[,j]<-paste0(Integrand2[,j]," * d",yuimaPPr@Kernel@variable.Integral@var.dx[j])
#   }
#   colnames(Integrand2) <- paste0("d",yuimaPPr@Kernel@variable.Integral@var.dx)
#   NamesIntegrandExpr <- as.character(matrix(colnames(Integrand2), dim(Integrand2)[1],dim(Integrand2)[2], byrow = TRUE))
#   Integrand2expr<- parse(text=Integrand2)
#
#   gridTime <- time(yuimaPPr@data@original.data)
#
#   yuimaPPr@Kernel@variable.Integral@var.dx
#   if(any(yuimaPPr@Kernel@variable.Integral@var.dx %in% yuimaPPr@model@solve.variable)){
#     my.envd1<-new.env()
#     ExistdN<-TRUE
#   }else{
#     ExistdN<-FALSE
#   }
#   Univariate<-FALSE
#   if(length(yuimaPPr@Ppr@counting.var)==1){
#     Univariate<-TRUE
#   }
#   if(any(!(yuimaPPr@Kernel@variable.Integral@var.dx %in% yuimaPPr@model@solve.variable))){
#     my.envd2<-new.env()
#     ExistdX<-TRUE
#   }else{
#     my.envd2<-new.env()
#     ExistdX<-FALSE
#   }
#
#   my.envd3 <- new.env()
#   namesparam<-names(param)
#   if(!(all(namesparam %in% yuimaPPr@Ppr@allparamPpr) && length(namesparam)==length(yuimaPPr@Ppr@allparamPpr))){
#     return(NULL)
#   }
#
#   # construction my.envd1
#   if(ExistdN){
#
#     #CountingVariable
#     for(i in c(1:length(yuimaPPr@Ppr@counting.var))){
#       cond <- yuimaPPr@Ppr@counting.var[i] %in% yuimaPPr@model@solve.variable
#       dummyData <-unique(yuimaPPr@data@original.data[,cond])[-1]
#       assign(yuimaPPr@Ppr@counting.var[i], dummyData,envir=my.envd1)
#     }
#     # Names expression
#     assign("NamesIntgra", NamesIntegrandExpr, envir=my.envd1)
#     #dN
#     namedX <-NULL
#     for(i in c(1:length(yuimaPPr@Kernel@variable.Integral@var.dx))){
#       if(yuimaPPr@Kernel@variable.Integral@var.dx[i] %in% yuimaPPr@Ppr@counting.var){
#         cond <- yuimaPPr@model@solve.variable %in% yuimaPPr@Kernel@variable.Integral@var.dx[i]
#         namedX<-c(namedX,paste0("d",yuimaPPr@Kernel@variable.Integral@var.dx[i]))
#         dummyData <- diff(as.numeric(yuimaPPr@data@original.data[,cond]))# We consider only Jump
#         dummyJumpTime <- gridTime[-1][dummyData>0]
#         dummyData2 <- diff(unique(cumsum(dummyData)))
#         dummyData3 <- zoo(dummyData2,order.by = dummyJumpTime)
#         assign(paste0("d",yuimaPPr@Kernel@variable.Integral@var.dx[i]), dummyData3 ,envir=my.envd1)
#       }
#     }
#     assign("namedX",namedX, envir = my.envd1)
#     assign("var.time",yuimaPPr@Kernel@variable.Integral@var.time,envir=my.envd1)
#     assign("t.time",yuimaPPr@Kernel@variable.Integral@upper.var,envir=my.envd1)
#
#     # Covariates
#     if(length(yuimaPPr@Ppr@covariates)>1){
#       # Covariates should be identified at jump time
#       return(NULL)
#     }
#
#   }
#   # end coonstruction my.envd1
#
#   # construction my.envd2
#   if(ExistdX){
#     #Covariate
#
#     #CountingVariable
#     for(i in c(1:length(yuimaPPr@Ppr@counting.var))){
#       cond <- yuimaPPr@Ppr@counting.var[i] %in% yuimaPPr@model@solve.variable
#       dummyData <-yuimaPPr@data@original.data[,cond]
#       assign(yuimaPPr@Ppr@counting.var[i], dummyData,envir=my.envd1)
#     }
#
#
#   }else{
#     assign("KerneldX",NULL,envir=my.envd2)
#   }
#
#   # end construction my.envd2
#
#   # construction my.envd3
#
#   #Covariate
#
#   #CountingVariable
#   for(i in c(1:length(yuimaPPr@Ppr@counting.var))){
#     cond <- yuimaPPr@Ppr@counting.var[i] %in% yuimaPPr@model@solve.variable
#     dummyData <-yuimaPPr@data@original.data[,cond]
#     assign(yuimaPPr@Ppr@counting.var[i], dummyData,envir=my.envd3)
#   }
#   #time
#   assign(yuimaPPr@model@time.variable, gridTime, my.envd3)
#
#   #Model
#   assign("YUIMA.PPR",yuimaPPr,envir=my.envd3)
#   assign("namesparam",namesparam,envir=my.envd3)
#   assign("gfun",gfun,envir=my.envd3)
#   assign("Integrand2",Integrand2,envir=my.envd3)
#   assign("Integrand2expr",Integrand2expr,envir=my.envd3)
#
#   assign("gridTime",gridTime,envir=my.envd3)
#   assign("Univariate",Univariate,envir=my.envd3)
#   assign("ExistdN",ExistdN,envir=my.envd3)
#   assign("ExistdX",ExistdX,envir=my.envd3)
#   out<-NULL
#
#   if(length(lower)==0 && length(upper)>0 && length(fixed)==0){
#     out <- optim(par=param, fn=Internal.LogLikPPR,
#                  my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
#                  method = method, upper=upper, ...)
#
#   }
#
#   if(length(lower)==0 && length(upper)==0 && length(fixed)>0){
#     out <- optim(par=param, fn=Internal.LogLikPPR,
#                  my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
#                  method = method, fixed = fixed, ...)
#
#   }
#
#
#   if(length(lower)>0 && length(upper)==0 && length(fixed)==0){
#     out <- optim(par = param,  fn=Internal.LogLikPPR,
#                  my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
#                  method = method, lower=lower, ...)
#   }
#
#   if(length(lower)>0 && length(upper)>0 && length(fixed)==0){
#     out <- optim(par=param, fn=Internal.LogLikPPR,
#                  my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
#                  method = method, upper = upper,
#                  lower=lower, ...)
#   }
#
#
#   if(length(lower)==0 && length(upper)>0 && length(fixed)>0){
#     out <- optim(par=param, fn=Internal.LogLikPPR,
#                  my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
#                  method = method, upper = upper,
#                  fixed = fixed,  ...)
#   }
#
#   if(length(lower)>0 && length(upper)==0 && length(fixed)>0){
#     out <- optim(par=param, fn=Internal.LogLikPPR,
#                  my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
#                  method = method, lower = lower,
#                  fixed = fixed, ...)
#   }
#
#
#   if(length(lower)>0 && length(upper)>0 && length(fixed)>0){
#     out <- optim(par=param, fn=Internal.LogLikPPR,
#                  my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
#                  method = method, lower = lower, fixed = fixed, upper = upper, ...)
#   }
#
#
#   if(is.null(out)){
#     out <- optim(par=param, fn=Internal.LogLikPPR,
#                  my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3,
#                  method = method, ...)
#   }
#
#
#   return(out)
#
# }


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
  # lambda <- gFunVect+cbind(0,IntKer)
  lambda <- gFunVect+IntKer
  time <- (c(lastEnv$s,lastEnv$t[1]))
  if(!logLikelihood){
    Intensity <- zoo(t(lambda), order.by = time)
    return(Intensity)
  }
  dn <- dim(lambda)
  if(dn[1]==1){
    #logLiklihood2 <- -sum(lambda*diff(time)[1])
    logLiklihood2 <- -1/2*sum((lambda[1:(length(lambda)-1)]+lambda[2:(length(lambda))])*diff(time)[1])

    logLiklihood1 <- sum(log(lambda)*lastEnv$CountVar)

    # logLiklihood2 <- -sum(lambda*diff(time))
    # dummyLamb <- lambda[lastEnv$CountVar]
    # #logLiklihood1 <- sum(log(dummyLamb[-length(dummyLamb)]))
    # logLiklihood1 <- sum(log(dummyLamb))


    # newlamb <- unique(as.numeric(lambda))
    #
    # cond <- as.numeric(lastEnv$CountVar)!=0
    # timeJ <- time[-1][cond]
    # timeJ1 <- unique(c(0,timeJ,lastEnv$t[1]))
    # logLiklihood2 <- -sum(newlamb*diff(timeJ1))
    # InternCount<-c(0:length(timeJ))
    # if((length(timeJ)+2)==length(timeJ1)){
    #   InternCount <- c(InternCount,tail(InternCount, n=1L))
    # }
    #
    # logLiklihood1 <- sum(log(newlamb)*diff(InternCount))


  }else{
    #### NO Rewrite

    # logLiklihood2 <- -rowSums(lambda[,-1]*diff(time)[1])
    # logLiklihood1 <- rowSums(log(lambda[,-1])*lastEnv$Integrator)
    logLiklihood2 <- -rowSums(lambda*diff(time)[1])
    #cond <- t(apply(as.matrix(lastEnv$CountVar),FUN = "diff",MARGIN = 2))!=0
    logLiklihood1 <- rowSums(log(lambda)*lastEnv$CountVar)
  }
  if(is.nan(logLiklihood1)){
    logLiklihood1 <- -10^10
  }
  if(is.nan(logLiklihood2)){
    logLiklihood2 <- -10^10
  }
  minusLoglik <- -sum(logLiklihood2+logLiklihood1)
   cat(sprintf("\n%.5f",minusLoglik))
   cat(sprintf("\n%.5f",param))
  return(minusLoglik)
}

DumFun<- function(X,Y){eval(X,envir=Y)}
