setPpr <- function(yuima, counting.var="N", gFun, Kernel,
                   var.dx, var.dt = "s", lambda.var = "lambda",
                   lower.var="0", upper.var = "t",
                   nrow =1 ,ncol=1){

  general <- TRUE
  ret <- aux.setPpr(yuima, counting.var, gFun,
    Kernel, var.dx, var.dt, lambda.var,
    lower.var="0", upper.var = "t",
    nrow =1 ,ncol=1,general =general)
  return(ret)
}


aux.setPpr <-function(yuima, counting.var="N", gFun, Kernel,
                      var.dx, var.dt = "s", lambda.var = "lambda",
                      lower.var="0", upper.var = "t",
                      nrow =1 ,ncol=1, general){
  g <- setMap(func = gFun, yuima = yuima,
              nrow = nrow, ncol = ncol)

  yuimadum <- yuima
  yuimadum@time.variable <- var.dt

  HawkesType <- FALSE
  if(counting.var %in% var.dx){
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

  #   IntPpr<- yuima:::setIntegral(yuima=yuimadum,
  #     integrand y= Kernel, var.dx = "N",
  #     lower.var = lower.var, upper.var = upper.var,
  #     out.var = "", nrow = nrow, ncol = ncol)

  #   return(list(Count.Proc = counting.var,
  #     gFun = list(param=g@Output@param, output=g@Output),
  #     Kernel = Integral, paramHawkes = paramHawkes,
  #     model = yuima, SelfEx = HawkesType))
  yuima1 <- setYuima(model =yuima)
  type <- yuima@measure.type
  if(type == "code"){
    measure <- list(df = as.character(yuima@measure$df$expr))
  }
  if(type == "CP"){
    measure <- list(intensity = as.character(yuima@measure$intensity),
                    df= as.character(yuima@measure$df$expr))
  }
  if(general){
    covariates<-character()
    if(sum(!(counting.var==yuima@solve.variable))!=0){
      covariates<-yuima@solve.variable[!(counting.var==yuima@solve.variable)]
    }
    Ppr <- new("info.Ppr",
               allparam = paramHawkes$allparam,
               allparamPpr = unique(c(paramHawkes$gFun,paramHawkes$Kern)),
               common = paramHawkes$common,
               counting.var = counting.var,
               var.dx = var.dx,
               upper.var = upper.var,
               lower.var = lower.var,
               covariates = covariates,
               var.dt = var.dt,
               additional.info = lambda.var,
               Info.measure = list(type=type,measure=measure))


    ret <- new("yuima.Ppr", Ppr = Ppr,
               gFun = g@Output,
               Kernel = Integral@Integral,
               yuima = yuima1)
  }else{

    Ppr <- new("info.Ppr",
               allparam = paramHawkes$allparam,
               allparamPpr = unique(c(paramHawkes$gFun,paramHawkes$Kern)),
               common = paramHawkes$common,
               counting.var = counting.var,
               var.dx = var.dx,
               upper.var = upper.var,
               lower.var = lower.var,
               covariates = character(),
               var.dt = var.dt,
               additional.info = "Exponential Hawkes",
               Info.measure = list(type=type,measure=measure))
    ret <- new("yuima.Hawkes", Ppr = Ppr,
               gFun = g@Output,
               Kernel = Integral@Integral,
               yuima = yuima1)

  }
  return(ret)
}
