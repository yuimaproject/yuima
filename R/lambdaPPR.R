# auxiliar function for the evaluation of g(t,X_t,N_t, theta)
internalGfunFromPPRModel <- function(gfun,my.envd3, univariate=TRUE){
  if(univariate){
    res<-as.numeric(eval(gfun, envir=my.envd3))
  }else{res<-NULL}
  return(res)
}

# auxiliar function for the evaluation of Kernel
InternalKernelFromPPRModel<-function(Integrand2,Integrand2expr,my.envd1=NULL,my.envd2=NULL,
                                     Univariate=TRUE, ExistdN, ExistdX, gridTime){
  if(Univariate){
    if(ExistdN){
      dimCol<- dim(Integrand2)[2]
      NameCol<-colnames(Integrand2)
      assign(my.envd1$t.time,gridTime, envir=my.envd1)
      IntegralKernel<- 0
      for(i in c(1:dimCol)){

        cond <- NameCol[i] %in% my.envd1$NamesIntgra
        assign(my.envd1$var.time, time(my.envd1[[my.envd1$namedX[cond]]]), my.envd1)
        # since it is just univariate we don't need a cycle for

        IntegralKernelDum<- sum(eval(Integrand2expr[cond], envir=my.envd1))
        IntegralKernel<-IntegralKernel+IntegralKernelDum
        #       cat("\n", IntegralKernel)
      }
    }

  }else{
    return(NULL)
  }

  return(IntegralKernel)
}

# auxiliar function for the evaluation of Intensity
InternalConstractionIntensity<-function(param,my.envd1=NULL,
                                        my.envd2=NULL,my.envd3=NULL){
  paramPPR <- my.envd3$YUIMA.PPR@PPR@allparamPPR
  namesparam <-my.envd3$namesparam


  gridTime  <-my.envd3$gridTime
  Univariate <-my.envd3$Univariate
  ExistdN <-my.envd3$ExistdN
  ExistdX <-my.envd3$ExistdX

  gfun<-my.envd3$gfun
  Integrand2<-my.envd3$Integrand2
  Integrand2expr<-my.envd3$Integrand2expr

  if(ExistdN){
    for(i in c(1:length(paramPPR))){
      cond<-namesparam %in% paramPPR[i]
      assign(paramPPR[i], param[cond], envir = my.envd1 )
    }
  }

  if(ExistdX){
    for(i in c(1:length(paramPPR))){
      cond<-namesparam %in% paramPPR[i]
      assign(paramPPR[i], param[cond], envir = my.envd2)
    }
  }

  #param
  for(i in c(1:length(paramPPR))){
    cond<-namesparam %in% paramPPR[i]
    assign(paramPPR[i], param[cond], envir = my.envd3)
  }


  KerneldN<- numeric(length=length(gridTime))
  for(i in c(1:length(gridTime))){
    KerneldN[i] <- InternalKernelFromPPRModel(Integrand2,Integrand2expr,my.envd1=my.envd1,my.envd2=my.envd2,
                                              Univariate=Univariate, ExistdN, ExistdX, gridTime=gridTime[i])
  }
  # KerneldN <- sapply(X=as.numeric(gridTime),FUN = InternalKernelFromPPRModel,
  #        Integrand2=Integrand2, Integrand2expr = Integrand2expr,my.envd1=my.envd1,my.envd2=my.envd2,
  #        Univariate=Univariate, ExistdN =ExistdN, ExistdX=ExistdX )
  KerneldCov<- numeric(length=length(gridTime))
  Evalgfun <- internalGfunFromPPRModel(gfun,my.envd3, univariate=Univariate)
  result<-KerneldN+KerneldCov+Evalgfun

}


InternalKernelFromPPRModel2<-function(Integrand2,Integrand2expr,my.envd1=NULL,my.envd2=NULL,
                                      Univariate=TRUE, ExistdN, ExistdX, gridTime){
  if(Univariate){
    if(ExistdN){
      dimCol<- dim(Integrand2)[2]
      NameCol<-colnames(Integrand2)
      assign(my.envd1$t.time,gridTime, envir=my.envd1)
      IntegralKernel<- 0
      for(i in c(1:dimCol)){

        # cond <- NameCol[i] %in% my.envd1$NamesIntgra
        # assign(my.envd1$var.time, time(my.envd1[[my.envd1$namedX[cond]]]), my.envd1)
        # since it is just univariate we don't need a cycle for
        cond <- paste0("JumpTime.",NameCol[i]) %in% my.envd1$namedJumpTimeX
        assign(my.envd1$var.time,my.envd1[[my.envd1$namedJumpTimeX[cond]]],envir=my.envd1)

        IntegralKernelDum<- sum(eval(Integrand2expr[cond], envir=my.envd1),na.rm = TRUE)
        IntegralKernel<-IntegralKernel+IntegralKernelDum
        #       cat("\n", IntegralKernel)
      }
    }

  }else{
    return(NULL)
  }

  return(IntegralKernel)
}


InternalConstractionIntensity2<-function(param,my.envd1=NULL,
                                         my.envd2=NULL,my.envd3=NULL){
  paramPPR <- my.envd3$YUIMA.PPR@PPR@allparamPPR
  namesparam <-my.envd3$namesparam


  gridTime  <-my.envd3$gridTime
  Univariate <-my.envd3$Univariate
  ExistdN <-my.envd3$ExistdN
  ExistdX <-my.envd3$ExistdX

  gfun<-my.envd3$gfun
  Integrand2<-my.envd3$Integrand2
  Integrand2expr<-my.envd3$Integrand2expr

  if(ExistdN){
    for(i in c(1:length(paramPPR))){
      cond<-namesparam %in% paramPPR[i]
      assign(paramPPR[i], param[cond], envir = my.envd1 )
    }
  }

  if(ExistdX){
    for(i in c(1:length(paramPPR))){
      cond<-namesparam %in% paramPPR[i]
      assign(paramPPR[i], param[cond], envir = my.envd2)
    }
  }

  #param
  for(i in c(1:length(paramPPR))){
    cond<-namesparam %in% paramPPR[i]
    assign(paramPPR[i], param[cond], envir = my.envd3)
  }


  KerneldN<- numeric(length=length(gridTime))
  # for(i in c(1:length(gridTime))){
  #   KerneldN[i] <- InternalKernelFromPPRModel(Integrand2,Integrand2expr,my.envd1=my.envd1,my.envd2=my.envd2,
  #                                             Univariate=Univariate, ExistdN, ExistdX, gridTime=gridTime[i])
  # }
  KerneldN <- sapply(X=gridTime,FUN = InternalKernelFromPPRModel2,
                     Integrand2=Integrand2, Integrand2expr = Integrand2expr,my.envd1=my.envd1,my.envd2=my.envd2,
                     Univariate=Univariate, ExistdN =ExistdN, ExistdX=ExistdX )
  KerneldCov<- numeric(length=length(gridTime))
  Evalgfun <- internalGfunFromPPRModel(gfun,my.envd3, univariate=Univariate)
  result<-KerneldN+KerneldCov+Evalgfun

}


Intensity.PPR <- function(yuimaPPR,param){
  # I need three envirnment
  # 1. my.envd1 is used when the counting variable is integrator
  # 2. my.envd2 is used when covariates or time are integrator
  # 3. my.envd3 for gfun

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
      condTime <- gridTime %in% my.envd1$JumpTime.dN
      dummyData <- yuimaPPR@data@original.data[condTime,cond]
      assign(yuimaPPR@PPR@counting.var[i], as.numeric(dummyData),envir=my.envd1)
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
        #dummyJumpTime <- gridTime[-1][dummyData>0]
        dummyJumpTime <- gridTime[-1][dummyData!=0]
        dummyData2 <- diff(unique(cumsum(dummyData)))
        #dummyData3 <- zoo(dummyData2,order.by = dummyJumpTime)
        dummyData3 <- dummyData2
        JumpTime <- dummyJumpTime
        assign(paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), as.numeric(dummyData3!=0) ,envir=my.envd1)
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
      # return(NULL)
      for(i in c(1:length(yuimaPPR@PPR@covariates))){
        cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@covariates[i]  
        dummyData <-yuimaPPR@data@original.data[,cond]
        assign(yuimaPPR@PPR@covariates[i], dummyData,envir=my.envd1)
      }
    }

  }
  # end coonstruction my.envd1

  # construction my.envd2
  if(ExistdX){
    #Covariate

    #CountingVariable
    for(i in c(1:length(yuimaPPR@PPR@counting.var))){
      cond <- yuimaPPR@PPR@counting.var[i] %in% yuimaPPR@model@solve.variable
      dummyData <-yuimaPPR@data@original.data[,cond]
      assign(yuimaPPR@PPR@counting.var[i], dummyData,envir=my.envd1)
    }


  }else{
    assign("KerneldX",NULL,envir=my.envd2)
  }

  # end construction my.envd2

  # construction my.envd3

  #Covariate
  dimCov <- length(yuimaPPR@PPR@covariates)
  if(dimCov>0){
    for(j in c(1:dimCov)){
      cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@covariates[j]
      dummyData <-yuimaPPR@data@original.data[,cond]
      assign(yuimaPPR@PPR@covariates[j], dummyData,envir=my.envd3)  
    }
  }

  #CountingVariable
  for(i in c(1:length(yuimaPPR@PPR@counting.var))){
    cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@counting.var[i]
    dummyData <-yuimaPPR@data@original.data[,cond]
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


  # end construction my.envd3





  ################################
  # Start Intensity construction #
  ################################
  #We define the parameters value for each environment

  param<-unlist(param)
  result<-InternalConstractionIntensity2(param,my.envd1,
                                           my.envd2,my.envd3)
  Int2<-zoo(as.matrix(result),order.by = gridTime)
  colnames(Int2)<-"lambda"
  res<-setData(Int2)
  return(res)


}





# Intensity.PPR <- function(yuimaPPR,param){
#   # I need three envirnment
#   # 1. my.envd1 is used when the counting variable is integrator
#   # 2. my.envd2 is used when covariates or time are integrator
#   # 3. my.envd3 for gfun
#
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
#
#
#   # end construction my.envd3
#
#
#
#
#
# ################################
# # Start Intensity construction #
# ################################
#   #We define the parameters value for each environment
#
#   param<-unlist(param)
#
#   result<-InternalConstractionIntensity(param,my.envd1,
#                                         my.envd2,my.envd3)
#
#   Int2<-zoo(as.matrix(result),order.by = gridTime)
#   colnames(Int2)<-"lambda"
#   res<-setData(Int2)
#   return(res)
#
#
# }
