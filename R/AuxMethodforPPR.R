## Here we write all auxiliar functions for the Point Process
## Regression Model
is.PPR <- function(yuimaPPR){is(yuimaPPR,"yuima.PPR")}

Internal.LogLikPPR <- function(param,my.envd1=NULL,
                               my.envd2=NULL,my.envd3=NULL){
  param<-unlist(param)

  IntLambda<-InternalConstractionIntensity2(param,my.envd1,
                                           my.envd2,my.envd3)

  # IntLambda<-InternalConstractionIntensity(param,my.envd1,
  #                                          my.envd2,my.envd3)
  Index<-my.envd3$gridTime
  Integr1a <- -sum(IntLambda[-length(IntLambda)]*diff(Index),na.rm=TRUE)
  Integr1b <- -sum(IntLambda[-1]*diff(Index),na.rm=TRUE)
  Integr1 <- (Integr1a+Integr1b)/2
  # if(is.nan(Integr1)){
  #   Integr1 <- -10^6
  # }
  if(length(my.envd3$YUIMA.PPR@PPR@counting.var)>0){
    cond1 <- my.envd3$YUIMA.PPR@model@solve.variable %in% my.envd3$YUIMA.PPR@PPR@counting.var
    cond2 <- diff(as.numeric(my.envd3$YUIMA.PPR@data@original.data[,cond1]))
    #Integr2<- sum(log(IntLambda[-1][cond2!=0]),na.rm=TRUE)
    Integr2 <- sum(log(IntLambda[cond2!=0]),na.rm=TRUE)
    #Integr2 <- (Integr2a+Integr2b)/2
  }else{
    yuima.stop("Spal")
  }
  # if(is.nan(Integr2)){
  #   Integr2 <- -10^6
  # }
  logLik <- Integr1+Integr2
  if(is.null(my.envd1$oldpar)){
    oldpar <- param
  }else{
    oldpar <- my.envd1$oldpar
  }
  ret <- -logLik/sum(cond2,na.rm=TRUE)#+sum((param-oldpar)^2*param^2)/2
  # line 40 necessary for the development of the cod
  #cat("\n ",logLik, param)
  
  #assign("oldpar",param,envir = my.envd1)
  
  return(ret)
}


quasiLogLik.PPR <- function(yuimaPPR, parLambda=list(), method=method, fixed = list(),
                            lower, upper, call, ...){

  yuimaPPR->yuimaPPR
  parLambda->param
  gfun<-yuimaPPR@gFun@formula

  gfun<-yuimaPPR@gFun@formula

  dimIntegr <- length(yuimaPPR@Kernel@Integrand@IntegrandList)
  Integrand2 <- character(length=dimIntegr)
  for(i in c(1:dimIntegr)){
    Integrand1 <- as.character(yuimaPPR@Kernel@Integrand@IntegrandList[[i]])
    timeCond <- paste0(" * (",yuimaPPR@Kernel@variable.Integral@var.time," < ",yuimaPPR@Kernel@variable.Integral@upper.var,")")
    Integrand2[i] <-paste0(Integrand1,timeCond)
  }

  Integrand2<- matrix(Integrand2,yuimaPPR@Kernel@Integrand@dimIntegrand[1],yuimaPPR@Kernel@Integrand@dimIntegrand[2])


  for(j in c(1:yuimaPPR@Kernel@Integrand@dimIntegrand[2])){
    Integrand2[,j]<-paste0(Integrand2[,j]," * d",yuimaPPR@Kernel@variable.Integral@var.dx[j])
  }
  colnames(Integrand2) <- paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx)
  NamesIntegrandExpr <- as.character(matrix(colnames(Integrand2), dim(Integrand2)[1],dim(Integrand2)[2], byrow = TRUE))
  Integrand2expr<- parse(text=Integrand2)

  gridTime <- time(yuimaPPR@data@original.data)

  yuimaPPR@Kernel@variable.Integral@var.dx
  if(any(yuimaPPR@Kernel@variable.Integral@var.dx %in% yuimaPPR@model@solve.variable)){
    my.envd1<-new.env()
    ExistdN<-TRUE
  }else{
    ExistdN<-FALSE
  }
  Univariate<-FALSE
  if(length(yuimaPPR@PPR@counting.var)==1){
    Univariate<-TRUE
  }
  if(any(!(yuimaPPR@Kernel@variable.Integral@var.dx %in% yuimaPPR@model@solve.variable))){
    my.envd2<-new.env()
    ExistdX<-TRUE
  }else{
    my.envd2<-new.env()
    ExistdX<-FALSE
  }

  my.envd3 <- new.env()
  namesparam<-names(param)
  if(!(all(namesparam %in% yuimaPPR@PPR@allparamPPR) && length(namesparam)==length(yuimaPPR@PPR@allparamPPR))){
    return(NULL)
  }

  # construction my.envd1
  if(ExistdN){

    #CountingVariable
    for(i in c(1:length(yuimaPPR@PPR@counting.var))){
      cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@counting.var[i]
      dummyData <-unique(yuimaPPR@data@original.data[,cond])[-1]
      assign(yuimaPPR@PPR@counting.var[i], rep(1,length(dummyData)),envir=my.envd1)
    }
    # Names expression
    assign("NamesIntgra", NamesIntegrandExpr, envir=my.envd1)
    #dN
    namedX <-NULL
    namedJumpTimeX <- NULL
    for(i in c(1:length(yuimaPPR@Kernel@variable.Integral@var.dx))){
      if(yuimaPPR@Kernel@variable.Integral@var.dx[i] %in% yuimaPPR@PPR@counting.var){
        cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@Kernel@variable.Integral@var.dx[i]
        namedX<-c(namedX,paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx[i]))
        namedJumpTimeX <-c(namedJumpTimeX,paste0("JumpTime.d",yuimaPPR@Kernel@variable.Integral@var.dx[i]))
        dummyData <- diff(as.numeric(yuimaPPR@data@original.data[,cond]))# We consider only Jump
        dummyJumpTime <- gridTime[-1][dummyData!=0]
        dummyData2 <- diff(unique(cumsum(dummyData)))
        #dummyData3 <- zoo(dummyData2,order.by = dummyJumpTime)
        dummyData3 <- rep(1,length(dummyData2))
        JumpTime <- dummyJumpTime
        assign(paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), dummyData3 ,envir=my.envd1)
        assign(paste0("JumpTime.d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), dummyJumpTime ,envir=my.envd1)
      }
    }
    assign("namedX",namedX, envir = my.envd1)
    assign("namedJumpTimeX",namedJumpTimeX, envir = my.envd1)
    assign("var.time",yuimaPPR@Kernel@variable.Integral@var.time,envir=my.envd1)
    assign("t.time",yuimaPPR@Kernel@variable.Integral@upper.var,envir=my.envd1)

    # Covariates
    if(length(yuimaPPR@PPR@covariates)>0){
      # Covariates should be identified at jump time
      for(i in c(1:length(yuimaPPR@PPR@covariates))){
        cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@covariates[i]
        condTime <- gridTime %in% my.envd1$JumpTime.dN
        assign(yuimaPPR@PPR@covariates[i],yuimaPPR@data@original.data[condTime,cond],envir = my.envd1)
      }
    }

  }
  # end coonstruction my.envd1

  # construction my.envd2
  if(ExistdX){
    #Covariate

    #CountingVariable
    for(i in c(1:length(yuimaPPR@PPR@counting.var))){
      cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@counting.var[i]
      dummyData <-yuimaPPR@data@original.data[,cond]
      assign(yuimaPPR@PPR@counting.var[i], dummyData,envir=my.envd1)
    }


  }else{
    assign("KerneldX",NULL,envir=my.envd2)
  }

  # end construction my.envd2

  # construction my.envd3

  #Covariate
  dimCov<-length(yuimaPPR@PPR@covariates)
  if(dimCov>0){
    for(i in c(1:dimCov)){
      cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@covariates[i] 
      dummyData <- yuimaPPR@data@original.data[,cond]
      assign(yuimaPPR@PPR@covariates[i], dummyData,envir=my.envd3)
    }
  }
  #CountingVariable
  for(i in c(1:length(yuimaPPR@PPR@counting.var))){
    cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@counting.var[i]
    dummyData <-cumsum(c(as.numeric(yuimaPPR@data@original.data[1,cond]!=0),as.numeric(diff(yuimaPPR@data@original.data[,cond])!=0)))
    assign(yuimaPPR@PPR@counting.var[i], dummyData,envir=my.envd3)
  }
  #time
  assign(yuimaPPR@model@time.variable, gridTime, my.envd3)

  #Model
  assign("YUIMA.PPR",yuimaPPR,envir=my.envd3)
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
   Hessian <- tryCatch(optimHess(as.list(out$par),
                 fn=Internal.LogLikPPR,
                 my.envd1=my.envd1,my.envd2=my.envd2,my.envd3=my.envd3),
                 error=function(){NULL})
   
   cond1 <- my.envd3$YUIMA.PPR@model@solve.variable %in% my.envd3$YUIMA.PPR@PPR@counting.var
   cond2 <- diff(as.numeric(my.envd3$YUIMA.PPR@data@original.data[,cond1]))
   N.jump <- sum(cond2,na.rm=TRUE)
   if(is.null(Hessian)){
      vcov <- matrix(NA,length(out$par),length(out$par))
   }else{
      vcov <- solve(Hessian)/N.jump
   }
   minuslog <- out$value*N.jump
   final_res<-new("yuima.PPR.qmle", call = call, coef = out$par, fullcoef = out$par,
                  vcov = vcov, min = minuslog, details = out, minuslogl = Internal.LogLikPPR,
                  method = method, nobs=as.integer(N.jump), model=my.envd3$YUIMA.PPR)
  return(final_res)

}


# quasiLogLik.PPR <- function(yuimaPPR, parLambda=list(), method=method, fixed = list(),
#                             lower, upper, call, ...){
#
#   yuimaPPR->yuimaPPR
#   parLambda->param
#   gfun<-yuimaPPR@gFun@formula
#
#   dimIntegr <- length(yuimaPPR@Kernel@Integrand@IntegrandList)
#   Integrand2 <- character(length=dimIntegr)
#   for(i in c(1:dimIntegr)){
#     Integrand1 <- as.character(yuimaPPR@Kernel@Integrand@IntegrandList[[i]])
#     timeCond <- paste0(" * (",yuimaPPR@Kernel@variable.Integral@var.time," < ",yuimaPPR@Kernel@variable.Integral@upper.var,")")
#     Integrand2[i] <-paste0(Integrand1,timeCond)
#   }
#
#   Integrand2<- matrix(Integrand2,yuimaPPR@Kernel@Integrand@dimIntegrand[1],yuimaPPR@Kernel@Integrand@dimIntegrand[2])
#
#
#   for(j in c(1:yuimaPPR@Kernel@Integrand@dimIntegrand[2])){
#     Integrand2[,j]<-paste0(Integrand2[,j]," * d",yuimaPPR@Kernel@variable.Integral@var.dx[j])
#   }
#   colnames(Integrand2) <- paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx)
#   NamesIntegrandExpr <- as.character(matrix(colnames(Integrand2), dim(Integrand2)[1],dim(Integrand2)[2], byrow = TRUE))
#   Integrand2expr<- parse(text=Integrand2)
#
#   gridTime <- time(yuimaPPR@data@original.data)
#
#   yuimaPPR@Kernel@variable.Integral@var.dx
#   if(any(yuimaPPR@Kernel@variable.Integral@var.dx %in% yuimaPPR@model@solve.variable)){
#     my.envd1<-new.env()
#     ExistdN<-TRUE
#   }else{
#     ExistdN<-FALSE
#   }
#   Univariate<-FALSE
#   if(length(yuimaPPR@PPR@counting.var)==1){
#     Univariate<-TRUE
#   }
#   if(any(!(yuimaPPR@Kernel@variable.Integral@var.dx %in% yuimaPPR@model@solve.variable))){
#     my.envd2<-new.env()
#     ExistdX<-TRUE
#   }else{
#     my.envd2<-new.env()
#     ExistdX<-FALSE
#   }
#
#   my.envd3 <- new.env()
#   namesparam<-names(param)
#   if(!(all(namesparam %in% yuimaPPR@PPR@allparamPPR) && length(namesparam)==length(yuimaPPR@PPR@allparamPPR))){
#     return(NULL)
#   }
#
#   # construction my.envd1
#   if(ExistdN){
#
#     #CountingVariable
#     for(i in c(1:length(yuimaPPR@PPR@counting.var))){
#       cond <- yuimaPPR@PPR@counting.var[i] %in% yuimaPPR@model@solve.variable
#       dummyData <-unique(yuimaPPR@data@original.data[,cond])[-1]
#       assign(yuimaPPR@PPR@counting.var[i], dummyData,envir=my.envd1)
#     }
#     # Names expression
#     assign("NamesIntgra", NamesIntegrandExpr, envir=my.envd1)
#     #dN
#     namedX <-NULL
#     for(i in c(1:length(yuimaPPR@Kernel@variable.Integral@var.dx))){
#       if(yuimaPPR@Kernel@variable.Integral@var.dx[i] %in% yuimaPPR@PPR@counting.var){
#         cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@Kernel@variable.Integral@var.dx[i]
#         namedX<-c(namedX,paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx[i]))
#         dummyData <- diff(as.numeric(yuimaPPR@data@original.data[,cond]))# We consider only Jump
#         dummyJumpTime <- gridTime[-1][dummyData>0]
#         dummyData2 <- diff(unique(cumsum(dummyData)))
#         dummyData3 <- zoo(dummyData2,order.by = dummyJumpTime)
#         assign(paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), dummyData3 ,envir=my.envd1)
#       }
#     }
#     assign("namedX",namedX, envir = my.envd1)
#     assign("var.time",yuimaPPR@Kernel@variable.Integral@var.time,envir=my.envd1)
#     assign("t.time",yuimaPPR@Kernel@variable.Integral@upper.var,envir=my.envd1)
#
#     # Covariates
#     if(length(yuimaPPR@PPR@covariates)>1){
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
#     for(i in c(1:length(yuimaPPR@PPR@counting.var))){
#       cond <- yuimaPPR@PPR@counting.var[i] %in% yuimaPPR@model@solve.variable
#       dummyData <-yuimaPPR@data@original.data[,cond]
#       assign(yuimaPPR@PPR@counting.var[i], dummyData,envir=my.envd1)
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
#   for(i in c(1:length(yuimaPPR@PPR@counting.var))){
#     cond <- yuimaPPR@PPR@counting.var[i] %in% yuimaPPR@model@solve.variable
#     dummyData <-yuimaPPR@data@original.data[,cond]
#     assign(yuimaPPR@PPR@counting.var[i], dummyData,envir=my.envd3)
#   }
#   #time
#   assign(yuimaPPR@model@time.variable, gridTime, my.envd3)
#
#   #Model
#   assign("YUIMA.PPR",yuimaPPR,envir=my.envd3)
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


lambdaFromData <- function(yuimaPPR, PPRData=NULL, parLambda=list()){
 if(is.null(PPRData)){
   PPRData<-yuimaPPR@data
 }else{
  # checklambdaFromData(yuimaPPR,PPRData)
 }
 if(!any(names(parLambda) %in% yuimaPPR@PPR@allparamPPR)){yuima.stop("1 ...")}
 if(!any(yuimaPPR@PPR@allparamPPR %in% names(parLambda))){yuima.stop("2 ...")}
  Time <- index(yuimaPPR@data@zoo.data[[1]])
  envPPR <- list()
  dY <- paste0("d",yuimaPPR@PPR@var.dx)
  for(i in (c(1:(length(Time)-1)))){
   envPPR[[i]]<-new.env()
    assign(yuimaPPR@gFun@param@time.var,rep(Time[i+1],i),envir=envPPR[[i]])
    assign(yuimaPPR@PPR@var.dt,Time[1:i],envir=envPPR[[i]])
    if(length(yuimaPPR@PPR@covariates)>0){
      for(j in c(1:length(yuimaPPR@PPR@covariates))){
        cond<-colnames(yuimaPPR@data@original.data)%in%yuimaPPR@PPR@covariates[[j]]
          assign(yuimaPPR@PPR@covariates[[j]],
                as.numeric(yuimaPPR@data@original.data[1:(i+1),cond]),
                envir=envPPR[[i]])
      }
    }
    for(j in c(1:length(yuimaPPR@PPR@counting.var))){
      cond<-colnames(yuimaPPR@data@original.data)%in%yuimaPPR@PPR@counting.var[[j]]
      assign(yuimaPPR@PPR@counting.var[[j]],
             as.numeric(yuimaPPR@data@original.data[1:(i+1),cond]),
             envir=envPPR[[i]])
    }

    for(j in c(1:length(yuimaPPR@PPR@var.dx))){
      cond<-c(colnames(yuimaPPR@data@original.data),yuimaPPR@PPR@var.dt)%in%c(yuimaPPR@PPR@var.dx,yuimaPPR@PPR@var.dt)[[j]]
      if(any(cond[-length(cond)])){
      assign(paste0("d",yuimaPPR@PPR@var.dx[[j]]),
             diff(as.numeric(yuimaPPR@data@original.data[1:(i+1),cond[-length(cond)]])),
             envir=envPPR[[i]])
      }
      if(tail(cond,n=1L)){
        assign(paste0("d",yuimaPPR@PPR@var.dx[[j]]),
               as.numeric(diff(Time[1:(i+1),cond[-length(cond)]])),
               envir=envPPR[[i]])
      }
    }
  }
  IntKernExpr<- function(Kern,dy){
    dum<-paste(Kern,dy,sep="*")
    dum<-paste0(dum, collapse = " + ")
    dum <- paste0("sum( ", dum, " )")
    return(parse(text = dum))
  }

  IntegKern <- lapply(yuimaPPR@Kernel@Integrand@IntegrandList,IntKernExpr,dY)
  Integrator <- t(as.matrix(eval(parse(text=dY[1]),envir=envPPR[[length(envPPR)]])))
  if(length(dY)>1){
    for(i in c(2:length(dY))){
      Integrator <- rbind(Integrator,
                          t(as.matrix(eval(parse(text=dY[1]),envir=envPPR[[length(envPPR)]]))))
    }
  }
  assign("Integrator",Integrator,envir=envPPR[[length(envPPR)]])
  assign("Nlamb",length(yuimaPPR@PPR@counting.var),envir=envPPR[[length(envPPR)]])
  res<-aux.lambdaFromData(param = unlist(parLambda), gFun=yuimaPPR@gFun,
    Kern =IntegKern, intensityParm = yuimaPPR@PPR@allparamPPR,
    envPPR)
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

aux.lambdaFromData <-function(param, gFun, Kern, intensityParm, envPPR,logLikelihood = FALSE){
  lapply(envPPR,FUN=dumFun2,Y=as.list(param))
  lastEnv <- tail(envPPR,n=1L)[[1]]
  # gFunVect<- t(simplify2array(my.lapply(gFun@formula, FUN=DumFun,
  #            Y = lastEnv), higher = (TRUE == "array")))

  gFunVect<- matrix(unlist(lapply(gFun@formula, FUN=DumFun,
     Y = lastEnv)),nrow=lastEnv$Nlamb,byrow=TRUE)
  # IntKer<- simplify2array(my.lapply(envPPR,function(my.env){
  #   t(simplify2array(my.lapply(Kern, FUN=DumFun,
  #             Y = my.env), higher = (TRUE == "array")))}),
  #   higher = (TRUE == "array")
  #   )
  # IntKer<- simplify2array(my.lapply(envPPR,myfun3,Kern=Kern),
  #                         higher = (TRUE == "array")
  # )

  # IntKer<- matrix(unlist(my.lapply(envPPR,myfun3,Kern=Kern)),
  #        nrow=1,byrow=TRUE)
  IntKer<- matrix(unlist(lapply(envPPR,myfun3,Kern=Kern)),
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



