setClass("yuima.th", slots=c(method = "character",
                                  up = "numeric", 
                                  low = "numeric",
                                  N = "numeric", 
                                  N_grid = "numeric", 
                                  regular_par = "ANY",
                                  h = "numeric",
                                  Dof = "character"),
         contains= "yuima.law")
setMethod("initialize", "yuima.th",
          function(.Object, h=1, method="LAG",  up=7, low=-7, 
                   N=180, N_grid=1000, regular_par=NULL){
            .Object@density<- dtLevy
            .Object@cdf<- ptLevy
            .Object@quantile<- qtLevy
            .Object@rng<- rtLevy
            .Object@method <- method
            .Object@up <- up
            .Object@low <- low
            .Object@N <- N
            .Object@N_grid <- N_grid
            .Object@regular_par <- regular_par
            .Object@h <- h
            .Object@Dof <- "nu"
            .Object@param.measure <- "nu"
           # validObject(.Object)
            return(.Object)
          }
)

setLaw_th <- function(h=1,method="LAG",  up=7, low=-7, 
                      N=180, N_grid=1000, regular_par=NULL, ...) {
  initialize(new("yuima.th"), h=h, method=method,  up=up, low=low, 
             N=N, N_grid=N_grid, regular_par=regular_par, ...)
}

setMethod("dens", "yuima.th", 
          function(object, x, param, log = FALSE, ...){
            if(is.list(param)){
              nu <- unlist(param)[object@param.measure]
            }else{
              nu <- param[object@param.measure]
            }
            res <- object@density(x = x, nu = nu, h = object@h, method = object@method,
                          up = object@up, low = object@low, N=object@N, 
                          N_grid = object@N_grid, regular_par = object@regular_par)
            if(log){
              res <- log(res)
            }
            return(res)
          } 
)

setMethod("cdf", "yuima.th", 
          function(object, q, param, ...){
            if(is.list(param)){
              nu <- unlist(param)[object@param.measure]
            }else{
              nu <- param[object@param.measure]
            }
            res <- object@cdf(q, nu = nu, h = object@h, method = object@method,
                                  up = object@up, low = object@low, N=object@N, 
                                  N_grid = object@N_grid, regular_par = object@regular_par)
            return(res)
          } 
)

setMethod("quant", "yuima.th", 
          function(object, p, param, ...){
            if(is.list(param)){
              nu <- unlist(param)[object@param.measure]
            }else{
              nu <- param[object@param.measure]
            }
            res <- object@quantile(p, nu = nu, h = object@h, method = object@method,
                              up = object@up, low = object@low, N=object@N, 
                              N_grid = object@N_grid, regular_par = object@regular_par)
            return(res)
          } 
)

setMethod("rand", "yuima.th", 
          function(object, n, param, ...){
            if(is.list(param)){
              nu <- unlist(param)[object@param.measure]
            }else{
              nu <- param[object@param.measure]
            }
            res <- numeric(length = n)
            res <- object@rng(n, nu = nu, h = object@h, method = object@method,
                                   up = object@up, low = object@low, N=object@N, 
                                   N_grid = object@N_grid, regular_par = object@regular_par)
            return(res)
          } 
)

setClass("yuima.LevyRM",
         slots=c(unit_Levy ="yuima.th",
                      regressors ="character",
                      LevyRM = "character",
                      paramRM = "character",
                      paramAll = "character",
                      solve.varRM = "character",
                      coeff = "character"),
         contains = "yuima")

setMethod("initialize", "yuima.LevyRM",
          function(.Object,  unit_Levy, yuima_regressors, LevyRM="Y", coeff=c("mu","sigma0"), 
                   data, sampling, characteristic, functional ){
            .Object@coeff <- c(paste0(coeff[1],1:length(yuima_regressors@solve.variable)),coeff[2])
            .Object@paramRM <- c(.Object@coeff,unit_Levy@param.measure)
            .Object@paramAll <- c(yuima_regressors@parameter@all, .Object@paramRM)
            .Object@solve.varRM <- c(yuima_regressors@solve.variable,LevyRM)
            .Object@LevyRM <- LevyRM
            .Object@regressors <- yuima_regressors@solve.variable
            .Object@unit_Levy <- unit_Levy
             test <- new("yuima", data = data, model = yuima_regressors, sampling = sampling, 
                   characteristic = characteristic, functional = functional)
            .Object@data <- test@data
            .Object@sampling  <- test@sampling
            .Object@model <- yuima_regressors
            .Object@characteristic <- test@characteristic
            .Object@functional <- test@functional 
            #validObject(.Object)
            return(.Object)
          }
)

setLRM <- function(unit_Levy, yuima_regressors, LevyRM="Y", coeff=c("mu","sigma0"), 
                   data=NULL, sampling = NULL, characteristic = NULL,
                   functional = NULL, ...) {
  new("yuima.LevyRM", unit_Levy=unit_Levy, yuima_regressors=yuima_regressors, 
             LevyRM = LevyRM, coeff = coeff, data = data, sampling = sampling, 
             characteristic = characteristic, functional = functional)
} 
#setMethod("simulate","yuima.LevyRM",
aux.simulateLevyRM <- function(object, nsim=1, seed=NULL, xinit, true.parameter, space.discretized = FALSE, 
          increment.W = NULL, increment.L = NULL, method = "euler", hurst, methodfGn = "WoodChan",
          sampling=sampling, subsampling=subsampling, ...){
            Samp <- sampling
            param <- true.parameter
            if(length(param[object@model@parameter@all])==0){
              traj <- simulate(object = object@model, nsim=nsim, seed=seed, 
                      xinit=xinit, space.discretized = space.discretized, 
                      increment.W = increment.W, increment.L = increment.L, 
                      method = method, hurst=hurst, methodfGn = methodfGn,
                      sampling=Samp, subsampling=subsampling)
            }else{
              traj <- simulate(object = object@model, nsim=nsim, seed=seed, 
                               true.parameter = param[object@model@parameter@all],
                               xinit=xinit, space.discretized = space.discretized, 
                               increment.W = increment.W, increment.L = increment.L, 
                               method = method, hurst=hurst, methodfGn = methodfGn,
                               sampling=Samp, subsampling=subsampling)
            }
            myint <- object@unit_Levy
            myint@h <- sampling@delta
            param <- unlist(param)
            # defaultW <- getOption("warn")
            # options(warn = -1)
            deltaLt <- rand(myint,n=sampling@n,param[myint@param.measure])
            process <- c(0,cumsum(deltaLt))
            test<- traj@data@original.data%*%param[object@paramRM[1:length(object@regressors)]]
            test<- test+param[object@paramRM[length(object@regressors)+1]]*process/sqrt(param[myint@param.measure])
            #plot(x=sampling@grid[[1]],y=as.numeric(test),type="l")
            res<-cbind(traj@data@original.data,test)
            colnames(res) <- object@solve.varRM 
            mydata <- setData(as.matrix(res),delta=sampling@delta, t0=sampling@Initial)
            object@data <- mydata
            object@unit_Levy@h<-1
            # res<-setLRM(object@unit_Levy, object@yuima_regressors, LevyRM="Y", coeff=c("mu","sigma0"), 
            #                    data=mydata, sampling = Samp, characteristic = characteristic,
            #                    functional = functional)
            # res<-new("yuima.LevyRM", unit_Levy=object@unit_Levy,
            #     yuima_regressors=object@yuima_regressors,
            #     LevyRM = object@LevyRM, coeff = c("mu","sigma"),
            #     data = mydata, sampling = Samp,
            #     characteristic = object@characteristic, functional = object@functional)
            return(object)
          }
#)


estimation_LRM <- function(start, model, data, upper, lower, PT=500, n_obs1=NULL){
  mydata <- data
  if(is.null(n_obs1)){
    n_obs1 <- floor((dim(mydata@original.data)[1]-1)/diff(range(index(mydata@original.data))))
  }
  Term <- tail(index(data@zoo.data[[1]]),1L)
  NofW <- n_obs1*Term # the number of the whole data
  m <- n_obs1*PT # the number of the data for the estimation of mu and sigma
  labY <- model@LevyRM
  Y <- as.numeric(mydata@original.data[,labY])
  regrlab <- model@regressors
  X <- as.matrix(mydata@original.data[,regrlab])
  
  dY <- Y[2:NofW] - Y[1:(NofW-1)]
  dX <- X[2:NofW,] - X[1:(NofW-1),]
  
  # first stage estimation
  # minus Cauchy quasi-likelihood 
  h<-1/n_obs1
  day <- n_obs1
  Term <- floor((dim(mydata@original.data)[1]-1)*h)
  if(length(regrlab)==1){
    mcql <- function(par, model, h, m){
      mu <- par[model@paramRM[1:length(model@regressors)]]
      sigma <- par[model@paramRM[length(model@regressors)+1]]
      #sum(log(sigma) + log(1 + (dY[1:m] - mu1*dX1[1:m] - mu2*dX2[1:m])^2/(h*sigma)^2))
      sum(log(sigma) + log(1 + (dY[1:m] - dX[1:m]*mu)^2/(h*sigma)^2))
    }
  }else{
    mcql <- function(par, model, h, m){
      mu <- par[model@paramRM[1:length(model@regressors)]]
      sigma <- par[model@paramRM[length(model@regressors)+1]]
      #sum(log(sigma) + log(1 + (dY[1:m] - mu1*dX1[1:m] - mu2*dX2[1:m])^2/(h*sigma)^2))
      sum(log(sigma) + log(1 + (dY[1:m] - dX[1:m,]%*%mu)^2/(h*sigma)^2))
    }
  }
  # estimation of mu and sigma (partial data)
  fres <- optim(par = unlist(start), fn = mcql, lower = lower, upper = upper, 
                method = "L-BFGS-B", model=model, h=h, m=m)
  esig <-  fres$par[model@paramRM[length(model@regressors)+1]]
  emu <- fres$par[model@paramRM[1:length(model@regressors)]]
  
  # second stage estimation
  # unit-time increments
  # if(is.null(data)){
  #   Term <- floor(tail(index(model@data@zoo.data[[1]]),1L))
  # }else{
  #   Term <- floor(tail(index(data@zoo.data[[1]]),1L))
  # }
  duY <- numeric(Term)
  duX <- matrix(NA, Term, length(model@regressors))
  #duX2 <- numeric(Term)
  duY[1] <- Y[day]
  duX[1,] <- X[day,]
  for(i in 2:Term){
    duY[i] <- Y[i*day] - Y[(i-1)*day + 1]
    duX[i,] <- X[i*day,] - X[(i-1)*day + 1,]
  }
  # unit-time residuals
  ures <- numeric(Term)
  for(i in 1:Term){
    ures[i] <- esig^(-1)*(duY[i] - duX[i,]%*%emu ) 
  }
  # hres <- numeric(m)
  # for(i in 1:m){
  #   hres[i] <- esig^(-1)/h*(dY[i] - dX[i,]%*%emu) 
  # }
  # minus (scaled) student-t likelihood function (unit time)
  stl <- function(nu){
    sum(-log(gamma((nu + 1)/2)) + log(gamma(nu/2)) + (nu + 1)/2*log(1 + (ures)^2))
    #-sum(log(dt(ures, nu)))
  }
  #sres <- optimize(stl, c(0.01, 10))
  sres<-optim(1,fn=stl,method="Brent",upper=100,lower=0.00001)
  enu <- sres$par
  names(enu)<-"nu"
  est<-c(emu, esig, enu)
  dFPos <- length(model@paramRM)
  # 
  scalepos <- length(model@paramRM)-1
  mycoef_a <- est[model@paramRM[-c(scalepos,dFPos)]]
  X_data <- data@original.data[,model@regressors]
  if(length(model@regressors)==1){
    VarX_data<- as.matrix(diff(X_data))
  }else{
    VarX_data<-apply(X_data,2,"diff")
  }
  #Nn <- dim(VarX_data)[1]
  Sn <- as.matrix(VarX_data[1,])%*%t(as.matrix(VarX_data[1,]))
  #m <- n_obs1*PT
  for(i in c(2:m)){
    Sn <- Sn+as.matrix(VarX_data[i,])%*%t(as.matrix(VarX_data[i,]))
  }
  Sn <- 1/(h^2*m)*Sn
  esig <- est[model@paramRM[scalepos]]
  GAM_a <- matrix(0, dim(Sn)[1]+1,dim(Sn)[2]+1)
  GAM_a[1:dim(Sn)[1],1:dim(Sn)[1]]<-Sn/(2*esig^2)
  GAM_a[dim(Sn)[1]+1,dim(Sn)[2]+1] <- 1/(2*esig^2)
  Vcov1 <- 1/m*solve(GAM_a)
  enu <- est[model@paramRM[dFPos]]
  GAM_nu <- 1/4*(trigamma(enu/2)-trigamma((enu+1)/2))
  Vcov2 <- 1/Term*solve(GAM_nu)
  Vcov <- matrix(0 , dim(Vcov1)[1]+1,dim(Vcov1)[2]+1)
  Vcov[1:dim(Vcov1)[1],1:dim(Vcov1)[1]]<- Vcov1
  Vcov[dim(Vcov1)[1]+1,dim(Vcov1)[1]+1]<- Vcov2
  call <- match.call()
  colnames(Vcov)<-names(est)
  rownames(Vcov)<-names(est)
  min <- c(fres$value,sres$value)
  final_res <- new("yuima.qmle", call = call, coef = est, fullcoef = est,
                 vcov = Vcov, min = min, details = list(), minuslogl = function(){NULL},
                 method = "L-BFGS-B", nobs=as.integer(NofW), model=model@model)
  #return(list(est=est,vcov=Vcov))
  return(final_res)
}