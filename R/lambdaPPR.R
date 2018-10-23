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

InternalConstractionIntensityFeedBackIntegrand<-function(param,my.envd1,
                                               my.envd2,my.envd3){
  paramPPR <- my.envd3$YUIMA.PPR@PPR@allparamPPR
  namesparam <-my.envd3$namesparam
  
  
  gridTime  <-my.envd3$gridTime
  Univariate <-my.envd3$Univariate
  ExistdN <-my.envd3$ExistdN
  ExistdX <-my.envd3$ExistdX
  
  gfun<-my.envd3$gfun
  
  allVarsInG<- all.vars(gfun)
  CondIntFeedBacksToG <- my.envd3$YUIMA.PPR@PPR@additional.info %in% allVarsInG
  
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
  
  
  # for(i in c(1:length(gridTime))){
  #   KerneldN[i] <- InternalKernelFromPPRModel(Integrand2,Integrand2expr,my.envd1=my.envd1,my.envd2=my.envd2,
  #                                             Univariate=Univariate, ExistdN, ExistdX, gridTime=gridTime[i])
  # }
  
  if(Univariate){
    NameCol <- colnames(Integrand2)
    
    # Kernel <- sapply(X=gridTime,FUN = InternalKernelFromPPRModel3,
    #                  Integrand2=Integrand2, Integrand2expr = Integrand2expr,my.envd1=my.envd1,my.envd2=my.envd2,
    #                  my.envd3=my.envd3,
    #                  Univariate=Univariate, ExistdN =ExistdN, ExistdX=ExistdX,
    #                  dimCol=dim(Integrand2)[2], NameCol = NameCol,
    #                  JumpTimeName =paste0("JumpTime.",NameCol))
    # Kernel <- evalKernelCpp2(Integrand2,
    #                          Integrand2expr,
    #                          my.envd1, my.envd2, my.envd3$YUIMA.PPR@PPR@IntensWithCount,
    #                          my.envd3$YUIMA.PPR@PPR@counting.var,
    #                          my.envd3$YUIMA.PPR@PPR@covariates,
    #                          ExistdN, ExistdX,
    #                          gridTime, dimCol = dim(Integrand2)[2], NameCol = NameCol,
    #                          JumpTimeName =paste0("JumpTime.",NameCol))
    # Evalgfun <- internalGfunFromPPRModel(gfun[i],my.envd3, univariate=TRUE)
    # result<-Kernel+Evalgfun
    kernel <- numeric(length = length(gridTime))
    Intensity <- numeric(length = length(gridTime))
    JumpTimeName <- paste0("JumpTime.", NameCol)
    dimCol <- dim(Integrand2)[2] 
    IntensityatJumpTime<-as.list(numeric(length=length(gridTime)))
    differentialCounting <- paste0("d",
            c(my.envd3$YUIMA.PPR@PPR@counting.var,
             my.envd3$YUIMA.PPR@PPR@additional.info))
    if(!CondIntFeedBacksToG){
      EvalGFUN <- eval(gfun,envir=my.envd3)
      Intensity[1]<-EvalGFUN[1]
    }
    aaaa<- length(gridTime)
    CondallJump <- rep(FALSE,aaaa+1)
    
    for(t in c(2:aaaa)){
      assign(my.envd1$t.time,gridTime[[t]][1],envir=my.envd1)
      for(j in c(1:dimCol)){
        if(NameCol[j] %in% differentialCounting){
          # Intensity at Jump Time
          assign(my.envd3$YUIMA.PPR@PPR@additional.info,
                 IntensityatJumpTime[[t-1]],
                 envir=my.envd1)
          # Counting Var at Jump Time
          assign(my.envd3$YUIMA.PPR@PPR@counting.var,
                 my.envd1[[my.envd1$PosListCountingVariable]][[t]],
                 envir=my.envd1)
          # Jump time <= t
          assign(my.envd1$var.time,my.envd1[[JumpTimeName[j]]][[t]],envir=my.envd1)
          KerneldN <- sum(eval(Integrand2expr,envir=my.envd1)*my.envd1[[my.envd1$namedX[j]]][[t]])
          kernel[t-1] <- KerneldN
        }
        
      }
        # Evaluation gFun
      if(!CondIntFeedBacksToG){
        #EvalGFUN <- eval(gfun,envir=my.envd3)
        if(t<=aaaa){
          Intensity[t]<- EvalGFUN[t-1]+kernel[t-1]
        }
      }else{
        # Here we evaluate gFun time by time
        
      }
      
      for(j in c(1:dimCol)){
        if(t+1<=aaaa){
          if(NameCol[j] %in% differentialCounting){
            # 
              CondallJump[t] <-my.envd3$JumpTimeLogical[t] 
              IntensityatJumpTime[[t]]<- Intensity[CondallJump[-1]] 
        
          }
        }
      }
      if(t==77){
        ff<-2
      }
    }
    
  }else{
    n.row <- length(my.envd3$YUIMA.PPR@PPR@counting.var)
    n.col <- length(gridTime)
    result <- matrix(NA,n.row, n.col) 
    #Kernel<- numeric(length=n.col-1)
    
    dimCol<- dim(Integrand2)[2]
    NameCol<-colnames(Integrand2)
    JumpTimeName  <- paste0("JumpTime.",NameCol)
    
    for(i in c(1:n.row)){
      
      # Kernel <- sapply(X=gridTime,FUN = InternalKernelFromPPRModel2,
      #                  Integrand2=t(Integrand2[i,]), Integrand2expr = Integrand2expr[[i]],my.envd1=my.envd1,my.envd2=my.envd2,
      #                  Univariate=TRUE, ExistdN =ExistdN, ExistdX=ExistdX, dimCol=dimCol, NameCol = NameCol,
      #                  JumpTimeName =JumpTimeName)
      # Kernel <- evalKernelCpp2(t(Integrand2[i,]), Integrand2expr[[i]],
      #                          my.envd1, my.envd2, my.envd3$YUIMA.PPR@PPR@IntensWithCount, 
      #                          my.envd3$YUIMA.PPR@PPR@counting.var,
      #                          my.envd3$YUIMA.PPR@PPR@covariates,
      #                          ExistdN, ExistdX,
      #                          gridTime, dimCol = dim(Integrand2)[2], NameCol = NameCol,
      #                          JumpTimeName =paste0("JumpTime.",NameCol))
      
      # Evalgfun <- internalGfunFromPPRModel(gfun[i],my.envd3, univariate=TRUE)
      # result[i,]<-Kernel+Evalgfun
    }
  }
  return(Intensity)
}

InternalKernelFromPPRModel2<-function(Integrand2,Integrand2expr,my.envd1=NULL,my.envd2=NULL,
                                      Univariate=TRUE, ExistdN, ExistdX, gridTime, dimCol, NameCol,
                                      JumpTimeName){
  
  if(Univariate){
      # JumpTimeName  <- paste0("JumpTime.",NameCol[i])
      # dimCol<- dim(Integrand2)[2]
      # NameCol<-colnames(Integrand2)
      if(ExistdN){
        assign(my.envd1$t.time,gridTime, envir=my.envd1)
      }
      if(ExistdX){
        assign(my.envd2$t.time,gridTime, envir=my.envd2)
      }
      
      IntegralKernel<- 0
      for(i in c(1:dimCol)){

        # cond <- NameCol[i] %in% my.envd1$NamesIntgra
        # assign(my.envd1$var.time, time(my.envd1[[my.envd1$namedX[cond]]]), my.envd1)
        # since it is just univariate we don't need a cycle for
        if(ExistdN){  
          # cond <- paste0("JumpTime.",NameCol[i]) %in% my.envd1$namedJumpTimeX
          # cond <- my.envd1$namedJumpTimeX %in% paste0("JumpTime.",NameCol[i])
          cond <- my.envd1$namedJumpTimeX %in% JumpTimeName[i]
          
          if(any(cond)){
            assign(my.envd1$var.time,my.envd1[[my.envd1$namedJumpTimeX[cond]]],envir=my.envd1)
            # condpos <- NameCol %in% my.envd1$namedX
            condpos <- my.envd1$namedX %in% NameCol[i]  
            if(any(condpos)){
              IntegralKernelDum<- sum(eval(Integrand2expr[condpos], envir=my.envd1),na.rm=TRUE)
              IntegralKernel<-IntegralKernel+IntegralKernelDum
            }
          }
        }
        
        if(ExistdX){  
          # cond <- paste0("JumpTime.",NameCol[i]) %in% my.envd2$namedJumpTimeX
          # cond <- my.envd2$namedJumpTimeX %in% paste0("JumpTime.",NameCol[i]) 
          cond <- my.envd2$namedJumpTimeX %in% JumpTimeName[i]
          if(any(cond)){
            assign(my.envd2$var.time,my.envd2[[my.envd2$namedJumpTimeX[cond]]],envir=my.envd2)
            # condpos <- NameCol %in% my.envd2$namedX
            condpos <- my.envd2$namedX %in% NameCol[i]
            if(any(condpos)){
              IntegralKernelDum<- sum(eval(Integrand2expr[condpos], envir=my.envd2))
              IntegralKernel<-IntegralKernel+IntegralKernelDum
            }
          }
        }

        
    }

  }else{
    return(NULL)
  }

  return(IntegralKernel)
}

InternalKernelFromPPRModel3<-function(Integrand2,Integrand2expr,my.envd1=NULL,my.envd2=NULL,my.envd3=NULL,
                                      Univariate=TRUE, ExistdN, ExistdX, gridTime, dimCol, NameCol,
                                      JumpTimeName){
  
  if(Univariate){
    # JumpTimeName  <- paste0("JumpTime.",NameCol[i])
    # dimCol<- dim(Integrand2)[2]
    # NameCol<-colnames(Integrand2)
    if(ExistdN){
      #assign(my.envd1$t.time,gridTime[1], envir=my.envd1)
      my.envd1[[my.envd1$t.time]]<-gridTime[1]
    }
    if(ExistdX){
      assign(my.envd2$t.time,gridTime[1], envir=my.envd2)
    }
    if(my.envd3$YUIMA.PPR@PPR@IntensWithCount){
      for(i in c(1:length(my.envd3$YUIMA.PPR@PPR@counting.var))){
        # assign(my.envd3$YUIMA.PPR@PPR@counting.var[i],
        #   my.envd1[[my.envd1$PosListCountingVariable[i]]][[gridTime[2]]]
        #   ,envir = my.envd1)
        my.envd1[[my.envd3$YUIMA.PPR@PPR@counting.var[i]]]<-my.envd1[[my.envd1$PosListCountingVariable[i]]][[gridTime[2]]]
      }
      if(length(my.envd3$YUIMA.PPR@PPR@covariates)>0){
        for(i in c(1:length(my.envd3$YUIMA.PPR@PPR@covariates))){
          # assign(my.envd3$YUIMA.PPR@PPR@covariates[i],
          #      my.envd1[[my.envd1$PosListCovariates[i]]][[gridTime[2]]]
          #      ,envir = my.envd1)
          my.envd1[[my.envd3$YUIMA.PPR@PPR@covariates[i]]]<-my.envd1[[my.envd1$PosListCovariates[i]]][[gridTime[2]]]
        }
      }
    }
    IntegralKernel<- 0
    for(i in c(1:dimCol)){
      
      # cond <- NameCol[i] %in% my.envd1$NamesIntgra
      # assign(my.envd1$var.time, time(my.envd1[[my.envd1$namedX[cond]]]), my.envd1)
      # since it is just univariate we don't need a cycle for
      if(ExistdN){  
        # cond <- paste0("JumpTime.",NameCol[i]) %in% my.envd1$namedJumpTimeX
        # cond <- my.envd1$namedJumpTimeX %in% paste0("JumpTime.",NameCol[i])
        cond <- my.envd1$namedJumpTimeX %in% JumpTimeName[i]
        
        if(any(cond)){
          assign(my.envd1$var.time,my.envd1[[my.envd1$namedJumpTimeX[cond]]][[gridTime[2]]],envir=my.envd1)
          # condpos <- NameCol %in% my.envd1$namedX
          #condpos <- my.envd1$namedX %in% NameCol[i]  
          condpos <- NameCol %in% NameCol[i]
          if(any(condpos)){
            InterDum <- eval(Integrand2expr[condpos], envir=my.envd1)*my.envd1[[NameCol[i]]][[gridTime[2]]]
            IntegralKernelDum<- sum(InterDum,na.rm=TRUE)
            IntegralKernel<-IntegralKernel+IntegralKernelDum
          }
        }
      }
      
      if(ExistdX){  
        # cond <- paste0("JumpTime.",NameCol[i]) %in% my.envd2$namedJumpTimeX
        # cond <- my.envd2$namedJumpTimeX %in% paste0("JumpTime.",NameCol[i]) 
        cond <- my.envd2$namedJumpTimeX %in% JumpTimeName[i]
        if(any(cond)){
          #assign(my.envd2$var.time,my.envd2[[my.envd2$namedJumpTimeX[cond]]],envir=my.envd2)
          assign(my.envd2$var.time,my.envd2[[my.envd2$namedJumpTimeX[cond]]][1:gridTime[2]],envir=my.envd2)
          # condpos <- my.envd2$namedX %in% NameCol  
          condpos <- NameCol %in% NameCol[i]
          if(any(condpos)){
            IntegralKernelDum<- sum(eval(Integrand2expr[condpos], envir=my.envd2)*my.envd2[[NameCol[i]]][1:gridTime[2]] , na.rm=TRUE)
            IntegralKernel<-IntegralKernel+IntegralKernelDum
          }
        }
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

  
  # for(i in c(1:length(gridTime))){
  #   KerneldN[i] <- InternalKernelFromPPRModel(Integrand2,Integrand2expr,my.envd1=my.envd1,my.envd2=my.envd2,
  #                                             Univariate=Univariate, ExistdN, ExistdX, gridTime=gridTime[i])
  # }
  
  if(Univariate){
    # Kernel<- numeric(length=length(gridTime))
    NameCol <- colnames(Integrand2)
    # Kernel <- sapply(X=gridTime,FUN = InternalKernelFromPPRModel2,
    #                  Integrand2=Integrand2, Integrand2expr = Integrand2expr,my.envd1=my.envd1,my.envd2=my.envd2,
    #                  Univariate=Univariate, ExistdN =ExistdN, ExistdX=ExistdX, 
    #                  dimCol=dim(Integrand2)[2], NameCol = NameCol,
    #                  JumpTimeName =paste0("JumpTime.",NameCol))

    # NameCol <- colnames(Integrand2)
    # Kernel <- evalKernelCpp(Integrand2, Integrand2expr,my.envd1, my.envd2,
    #                         ExistdN, ExistdX, gridTime, dim(Integrand2)[2], NameCol,
    #                         paste0("JumpTime.",NameCol))
    
    # Kernel <- sapply(X=gridTime,FUN = InternalKernelFromPPRModel3,
    #                  Integrand2=Integrand2, Integrand2expr = Integrand2expr,my.envd1=my.envd1,my.envd2=my.envd2,
    #                  my.envd3=my.envd3,
    #                  Univariate=Univariate, ExistdN =ExistdN, ExistdX=ExistdX,
    #                  dimCol=dim(Integrand2)[2], NameCol = NameCol,
    #                  JumpTimeName =paste0("JumpTime.",NameCol))
    Kernel <- evalKernelCpp2(Integrand2,
                      Integrand2expr,
                      my.envd1, my.envd2, my.envd3$YUIMA.PPR@PPR@IntensWithCount,
                      my.envd3$YUIMA.PPR@PPR@counting.var,
                      my.envd3$YUIMA.PPR@PPR@covariates,
                      ExistdN, ExistdX,
                      gridTime, dimCol = dim(Integrand2)[2], NameCol = NameCol,
                      JumpTimeName =paste0("JumpTime.",NameCol))
  #KerneldCov<- numeric(length=length(gridTime))
    Evalgfun <- internalGfunFromPPRModel(gfun,my.envd3, univariate=Univariate)
    result<-Kernel+Evalgfun
  }else{
    n.row <- length(my.envd3$YUIMA.PPR@PPR@counting.var)
    n.col <- length(gridTime)
    result <- matrix(NA,n.row, n.col) 
    #Kernel<- numeric(length=n.col-1)
    
     dimCol<- dim(Integrand2)[2]
     NameCol<-colnames(Integrand2)
     JumpTimeName  <- paste0("JumpTime.",NameCol)
    
    for(i in c(1:n.row)){
      # Kernel <- pvec(v=gridTime,FUN = InternalKernelFromPPRModel2,
      #                  Integrand2=t(Integrand2[i,]), Integrand2expr = Integrand2expr[[i]],my.envd1=my.envd1,my.envd2=my.envd2,
      #                  Univariate=TRUE, ExistdN =ExistdN, ExistdX=ExistdX, dimCol=dimCol, NameCol = NameCol,
      #                  JumpTimeName =JumpTimeName)
      # Kernel <- sapply(X=gridTime,FUN = InternalKernelFromPPRModel2,
      #                  Integrand2=t(Integrand2[i,]), Integrand2expr = Integrand2expr[[i]],my.envd1=my.envd1,my.envd2=my.envd2,
      #                  Univariate=TRUE, ExistdN =ExistdN, ExistdX=ExistdX, dimCol=dimCol, NameCol = NameCol,
      #                  JumpTimeName =JumpTimeName)
      # Kernel <- sapply(X=gridTime,FUN = evalKernelCppPvec,
      #                   Integrand2=t(Integrand2[i,]), Integrand2expr = Integrand2expr[[i]],
      #                   myenvd1=my.envd1,myenvd2=my.envd2,
      #                   ExistdN =ExistdN, ExistdX=ExistdX,
      #                   dimCol=dimCol, NameCol = NameCol,
      #                   JumpTimeName =JumpTimeName)
      # Kernel <- evalKernelCpp(t(Integrand2[i,]), Integrand2expr[[i]],my.envd1, my.envd2,
      #                ExistdN, ExistdX, gridTime, dimCol, NameCol,
      #                JumpTimeName)
      Kernel <- evalKernelCpp2(t(Integrand2[i,]), Integrand2expr[[i]],
                               my.envd1, my.envd2, my.envd3$YUIMA.PPR@PPR@IntensWithCount, 
                               my.envd3$YUIMA.PPR@PPR@counting.var,
                               my.envd3$YUIMA.PPR@PPR@covariates,
                               ExistdN, ExistdX,
                               gridTime, dimCol = dim(Integrand2)[2], NameCol = NameCol,
                               JumpTimeName =paste0("JumpTime.",NameCol))
      Evalgfun <- internalGfunFromPPRModel(gfun[i],my.envd3, univariate=TRUE)
      result[i,]<-Kernel+Evalgfun
    }
  }
  return(result)
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
    #Integrand1 <- as.character(yuimaPPR@Kernel@Integrand@IntegrandList[[i]])
    #timeCond <- paste0(" * (",yuimaPPR@Kernel@variable.Integral@var.time," < ",yuimaPPR@Kernel@variable.Integral@upper.var,")")
    #Integrand2[i] <-paste0(Integrand1,timeCond)
    Integrand2[i] <- as.character(yuimaPPR@Kernel@Integrand@IntegrandList[[i]])
  }

  Integrand2<- matrix(Integrand2,yuimaPPR@Kernel@Integrand@dimIntegrand[1],yuimaPPR@Kernel@Integrand@dimIntegrand[2])


  # for(j in c(1:yuimaPPR@Kernel@Integrand@dimIntegrand[2])){
  #   Integrand2[,j]<-paste0(Integrand2[,j]," * d",yuimaPPR@Kernel@variable.Integral@var.dx[j])
  # }
  colnames(Integrand2) <- paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx)
  NamesIntegrandExpr <- as.character(matrix(colnames(Integrand2), dim(Integrand2)[1],dim(Integrand2)[2], byrow = TRUE))
  # if(yuimaPPR@Kernel@Integrand@dimIntegrand[2]==1 & yuimaPPR@Kernel@Integrand@dimIntegrand[1]==1)
  #   Integrand2expr<- parse(text=Integrand2)
  # 
  # if(yuimaPPR@Kernel@Integrand@dimIntegrand[2]>1 & yuimaPPR@Kernel@Integrand@dimIntegrand[1]==1){
  #   dum <- paste0(Integrand2,collapse=" + ")
  #   Integrand2expr <- parse(text=dum)
  # }
  if(yuimaPPR@Kernel@Integrand@dimIntegrand[1]==1){
    Integrand2expr <- parse(text=Integrand2)
  }else{
    Integrand2expr <- list()
    for(hh in c(1:yuimaPPR@Kernel@Integrand@dimIntegrand[1])){
      Integrand2expr[[hh]] <- parse(text=Integrand2[hh,])
    }
  }
  
  gridTime <- time(yuimaPPR@data@original.data)

  yuimaPPR@Kernel@variable.Integral@var.dx
  if(any(yuimaPPR@Kernel@variable.Integral@var.dx %in% yuimaPPR@PPR@counting.var)){
    my.envd1<-new.env()
    ExistdN<-TRUE
  }else{
    ExistdN<-FALSE
  }
  Univariate<-FALSE
  if(length(yuimaPPR@PPR@counting.var)==1){
    Univariate<-TRUE
  }
  if(any(yuimaPPR@Kernel@variable.Integral@var.dx %in% yuimaPPR@PPR@covariates)){
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
    # for(i in c(1:length(yuimaPPR@PPR@counting.var))){
    #   cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@counting.var[i] 
    #   condTime <- gridTime %in% my.envd1$JumpTime.dN
    #   dummyData <- yuimaPPR@data@original.data[condTime,cond]
    #   assign(yuimaPPR@PPR@counting.var[i], as.numeric(dummyData),envir=my.envd1)
    # }
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
        Jump <- lapply(X=as.numeric(gridTime), FUN = function(X,JumpT,Jump){Jump[JumpT<X]},
                       JumpT = dummyJumpTime, Jump = as.numeric(dummyData3!=0))
        assign(paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), 
               Jump ,
               envir=my.envd1)
        dummyJumpTimeNew <- lapply(X=as.numeric(gridTime), FUN = function(X,JumpT){JumpT[JumpT<X]},
                                   JumpT = dummyJumpTime)
        assign(paste0("JumpTime.d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), dummyJumpTimeNew ,envir=my.envd1)
      }
    }
    assign("namedX",namedX, envir = my.envd1)
    assign("namedJumpTimeX",namedJumpTimeX, envir = my.envd1)
    assign("var.time",yuimaPPR@Kernel@variable.Integral@var.time,envir=my.envd1)
    assign("t.time",yuimaPPR@Kernel@variable.Integral@upper.var,envir=my.envd1)

    #CountingVariable
    PosListCountingVariable <- NULL
    for(i in c(1:length(yuimaPPR@PPR@counting.var))){
      cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@counting.var[i]
      JUMPTIME <- tail(my.envd1[[paste0("JumpTime.d",yuimaPPR@Kernel@variable.Integral@var.dx[i])]],1L)[[1]]
      condTime <- gridTime %in% JUMPTIME 
      
      dummyData <- yuimaPPR@data@original.data[condTime,cond]
      dummyDataA <- lapply(X=as.numeric(gridTime), FUN = function(X,JumpT,Jump){Jump[JumpT<X]},
             JumpT = JUMPTIME, Jump = dummyData)
      dummyList <- paste0("List_",yuimaPPR@PPR@counting.var[i])
      PosListCountingVariable <- c(PosListCountingVariable,dummyList)
      assign(dummyList, dummyDataA, envir=my.envd1)
      assign(yuimaPPR@PPR@counting.var[i], numeric(length=0L), envir=my.envd1)
    }
    assign("PosListCountingVariable", PosListCountingVariable, envir=my.envd1)
    # Covariates
    if(length(yuimaPPR@PPR@covariates)>0){
      # Covariates should be identified at jump time
      # return(NULL)
      PosListCovariates <- NULL
      for(i in c(1:length(yuimaPPR@PPR@covariates))){
        cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@covariates[i]  
        #dummyData <-yuimaPPR@data@original.data[,cond]
        dummyData <- yuimaPPR@data@original.data[condTime, cond]
        dummyDataB <- lapply(X=as.numeric(gridTime), FUN = function(X,JumpT,Jump){Jump[JumpT<X]},
              JumpT = JUMPTIME, Jump = dummyData)
        dummyListCov <- paste0("List_",yuimaPPR@PPR@covariates[i])
        PosListCovariates <- c(PosListCovariates,dummyListCov)
        assign(dummyListCov, dummyDataB,envir=my.envd1)
        assign(yuimaPPR@PPR@covariates[i], numeric(length=0L),envir=my.envd1)
      }
      assign("PosListCovariates", PosListCovariates,envir=my.envd1)
    }

  }
  # end coonstruction my.envd1

  # construction my.envd2
  if(ExistdX){
    #Covariate
    dummyData<-NULL
    #CountingVariable
    for(i in c(1:length(yuimaPPR@PPR@counting.var))){
      cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@counting.var[i]  
      dummyData <-as.numeric(yuimaPPR@data@original.data[,cond])
     # assign(yuimaPPR@PPR@counting.var[i], dummyData[-length(dummyData)],envir=my.envd2)
      assign(yuimaPPR@PPR@counting.var[i], dummyData,envir=my.envd2)
    }
    namedX<-NULL
    namedJumpTimeX<-NULL
    for(i in c(1:length(yuimaPPR@Kernel@variable.Integral@var.dx))){
      if(yuimaPPR@Kernel@variable.Integral@var.dx[i] %in% yuimaPPR@PPR@covariates){
        cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@Kernel@variable.Integral@var.dx[i]
        namedX<-c(namedX,paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx[i]))
        namedJumpTimeX <-c(namedJumpTimeX,paste0("JumpTime.d",yuimaPPR@Kernel@variable.Integral@var.dx[i]))
        dummyData <- diff(as.numeric(yuimaPPR@data@original.data[,cond]))# We consider only Jump
        #dummyJumpTime <- gridTime[-1][dummyData>0]
        #assign(paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), dummyData ,envir=my.envd2)
        assign(paste0("d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), c(0,dummyData) ,envir=my.envd2)
        #assign(paste0("JumpTime.d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), gridTime[-1] ,envir=my.envd2)
        assign(paste0("JumpTime.d",yuimaPPR@Kernel@variable.Integral@var.dx[i]), as.numeric(gridTime) ,envir=my.envd2)
      }
    }
    
    assign("namedX",namedX, envir = my.envd2)
    assign("namedJumpTimeX",namedJumpTimeX, envir = my.envd2)
    assign("var.time",yuimaPPR@Kernel@variable.Integral@var.time,envir=my.envd2)
    assign("t.time",yuimaPPR@Kernel@variable.Integral@upper.var,envir=my.envd2)
    
    for(i in c(1:length(yuimaPPR@PPR@covariates))){
      cond <- yuimaPPR@model@solve.variable %in% yuimaPPR@PPR@covariates[i]  
      #dummyData <-yuimaPPR@data@original.data[,cond]
      dummyData <-as.numeric(yuimaPPR@data@original.data[, cond])
      #assign(yuimaPPR@PPR@covariates[i], dummyData[-length(dummyData)],envir=my.envd2)
      assign(yuimaPPR@PPR@covariates[i], dummyData,envir=my.envd2)
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

  l1 =as.list(as.numeric(gridTime))
  l2 = as.list(c(1:(length(l1))))
  l3 = mapply(c, l1, l2, SIMPLIFY=FALSE)
  
  assign("gridTime",l3,envir=my.envd3)
  assign("Univariate",Univariate,envir=my.envd3)
  assign("ExistdN",ExistdN,envir=my.envd3)
  assign("ExistdX",ExistdX,envir=my.envd3)
  
  assign("JumpTimeLogical",c(FALSE,as.integer(diff(my.envd3$N))!=0),envir=my.envd3)

  # end construction my.envd3





  ################################
  # Start Intensity construction #
  ################################
  #We define the parameters value for each environment

  param<-unlist(param)
  if(my.envd3$YUIMA.PPR@PPR@additional.info %in% all.vars(my.envd3$Integrand2expr)){
    result <- InternalConstractionIntensityFeedBackIntegrand(param,my.envd1,
                                   my.envd2,my.envd3)
  }else{
    result<-InternalConstractionIntensity2(param,my.envd1,
                                           my.envd2,my.envd3)
  }
  if(class(result)=="matrix"){
    my.matr <- t(result)
    colnames(my.matr) <-paste0("lambda",c(1:yuimaPPR@Kernel@Integrand@dimIntegrand[1]))
    Int2<-zoo(my.matr,order.by = gridTime)
  }else{
    Int2<-zoo(as.matrix(result),order.by = gridTime)
    colnames(Int2)<-"lambda"
  }
  
  
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
