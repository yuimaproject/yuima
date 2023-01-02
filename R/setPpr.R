setPPR <- function(yuima, counting.var="N", gFun, Kernel,
                   var.dx= "s", var.dt = "s", lambda.var = "lambda",
                   lower.var="0", upper.var = "t",
                   nrow =1 ,ncol=1){

  general <- TRUE
  ret <- aux.setPPR(yuima, counting.var, gFun,
    Kernel, var.dx, var.dt, lambda.var,
    lower.var="0", upper.var = "t",
    nrow =1 ,ncol=1,general =general)
  return(ret)
}


aux.setPPR <-function(yuima, counting.var="N", gFun, Kernel,
                      var.dx, var.dt = "s", lambda.var = "lambda",
                      lower.var="0", upper.var = "t",
                      nrow =1 ,ncol=1, general){
  g <- setMap(func = gFun, yuima = yuima,
              nrow = nrow, ncol = ncol)

  yuimadum <- yuima
  yuimadum@time.variable <- var.dt

  HawkesType <- FALSE
  if(all(counting.var %in% var.dx)){
    HawkesType <- TRUE
  }
  if(!HawkesType){
    Integral <- setIntegral(yuima=yuimadum,
                            integrand = Kernel, var.dx = var.dx,
                            lower.var = lower.var, upper.var = upper.var,
                            out.var = "", nrow = nrow, ncol = ncol)
  }else{
    Integral <- setIntegral(yuima=yuimadum,
                            integrand = Kernel, var.dx = var.dx,
                            lower.var = lower.var, upper.var = upper.var,
                            out.var = "", nrow = nrow, ncol = ncol)
  }
  if(g@Output@dimension[1]!=Integral@Integral@Integrand@dimIntegrand[1]){
    yuima.stop("dimension gFun and kernel mismatch")
  }


  allparam <- unique(c(yuima@parameter@all, g@Output@param@allparamMap,
                       Integral@Integral@param.Integral@Integrandparam))
  common <- unique(c(g@Output@param@common,
                     Integral@Integral@param.Integral@common))
  paramHawkes <- list(allparam = allparam, common = common,
                      gFun = g@Output@param@allparamMap,
                      Kern = Integral@Integral@param.Integral@Integrandparam)

  #   IntPPR<- yuima:::setIntegral(yuima=yuimadum,
  #     integrand y= Kernel, var.dx = "N",
  #     lower.var = lower.var, upper.var = upper.var,
  #     out.var = "", nrow = nrow, ncol = ncol)

  #   return(list(Count.Proc = counting.var,
  #     gFun = list(param=g@Output@param, output=g@Output),
  #     Kernel = Integral, paramHawkes = paramHawkes,
  #     model = yuima, SelfEx = HawkesType))
  yuima1 <- setYuima(model =yuima)
  type <- yuima@measure.type
  if(all(type == "code")){
    if(!(is(yuima@measure$df,"yuima.law")))
      measure <- list(df = as.character(yuima@measure$df$expr))
  }else{
    measure <- yuima@measure
    }
  if(all(type == "CP")){
    if(!(is(yuima@measure$df,"yuima.law")))
      measure <- list(intensity = as.character(yuima@measure$intensity),
                    df= as.character(yuima@measure$df$expr))
  }else{
    measure <- yuima@measure
  }
  IntensWithCount<-FALSE
  if(!is.list(g@Output@formula)){
    if(any(counting.var%in%all.vars(g@Output@formula)))
       IntensWithCount<-TRUE
  }else{
    ddd<- length(g@Output@formula)
    for(i in c(1:ddd)){
      if(any(counting.var%in%all.vars(g@Output@formula[[i]])))
        IntensWithCount<-TRUE
    }
  }
  if(any(counting.var%in%Integral@Integral@variable.Integral@var.dx))
     IntensWithCount<-TRUE

  if(!is.list(Integral@Integral@Integrand@IntegrandList)){
    if(any(counting.var%in%all.vars(Integral@Integral@Integrand@IntegrandList)))
      IntensWithCount<-TRUE
  }else{
    ddd<- length(Integral@Integral@Integrand@IntegrandList)
    for(i in c(1:ddd)){
      if(any(counting.var%in%all.vars(Integral@Integral@Integrand@IntegrandList[[i]])))
        IntensWithCount<-TRUE
    }
  }

  RegressWithCount <- FALSE

  if(general){
    covariates<-character()
    if(sum(!(counting.var==yuima@solve.variable))!=0){
      condCovariate<-!(counting.var==yuima@solve.variable)
      covariates<-yuima@solve.variable[condCovariate]
      if(length(covariates)>0){
        covariate.drift <- yuima@drift[condCovariate]
        covariate.diff <- yuima@diffusion[condCovariate]
        covariate.jump <- yuima@jump.coeff[condCovariate]
      }
      if(any(counting.var %in% all.vars(covariate.drift))){
        RegressWithCount <-TRUE
      }

      ddd.dif <-length(covariate.diff)
      if(length(covariate.diff)>0){
        for(i in c(1:ddd.dif)){
          if(any(counting.var %in% all.vars(covariate.diff[[i]]))){
            RegressWithCount <-TRUE
          }
        }
      }

      ddd.jump <-length(covariate.jump)
      if(length(covariate.jump)>0){
        for(i in c(1:ddd.jump)){
          if(any(counting.var %in% all.vars(covariate.jump[[i]]))){
            RegressWithCount <-TRUE
          }
        }
      }

    }
    PPR <- new("info.PPR",
               allparam = paramHawkes$allparam,
               allparamPPR = unique(c(paramHawkes$gFun,paramHawkes$Kern)),
               common = paramHawkes$common,
               counting.var = counting.var,
               var.dx = var.dx,
               upper.var = upper.var,
               lower.var = lower.var,
               covariates = covariates,
               var.dt = var.dt,
               additional.info = lambda.var,
               Info.measure = list(type=type,measure=measure),
               RegressWithCount=RegressWithCount,
               IntensWithCount=IntensWithCount)


    ret <- new("yuima.PPR", PPR = PPR,
               gFun = g@Output,
               Kernel = Integral@Integral,
               yuima = yuima1)
  }else{

    PPR <- new("info.PPR",
               allparam = paramHawkes$allparam,
               allparamPPR = unique(c(paramHawkes$gFun,paramHawkes$Kern)),
               common = paramHawkes$common,
               counting.var = counting.var,
               var.dx = var.dx,
               upper.var = upper.var,
               lower.var = lower.var,
               covariates = character(),
               var.dt = var.dt,
               additional.info = "Exponential Hawkes",
               Info.measure = list(type=type,measure=measure))
    ret <- new("yuima.Hawkes", PPR = PPR,
               gFun = g@Output,
               Kernel = Integral@Integral,
               yuima = yuima1)

  }
  return(ret)
}
