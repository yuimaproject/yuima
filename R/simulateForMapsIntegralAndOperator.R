# Method for Map
setMethod("simulate", "yuima.Map",
                    function(object, nsim=1, seed=NULL, xinit, true.parameter,
                             space.discretized=FALSE, increment.W=NULL, increment.L=NULL, method="euler",
                             hurst, methodfGn="WoodChan",
                             sampling, subsampling,
                             #Initial = 0, Terminal = 1, n = 100, delta,
                             #	grid, random = FALSE, sdelta=as.numeric(NULL),
                             #	sgrid=as.numeric(NULL), interpolation="none"
                             ...){
                      res <- aux.simulatOutput(object, nsim = nsim, seed = seed,
                        xinit = xinit, true.parameter = true.parameter,
                        space.discretized = space.discretized, increment.W = increment.W,
                        increment.L = increment.L, method = method, hurst = hurst,
                        methodfGn = methodfGn, sampling = sampling, subsampling = subsampling)

                      return(res)
                    }
          )

#aux.simulatOutput<-function(yuima.output, param=list(), sampling){
 # param <- list(a=1,b=1,s1=0.1,s2=0.2,a1=0.1,a2=0.1)
aux.simulatOutput<-function(object, nsim, seed, xinit,
  true.parameter, space.discretized, increment.W, increment.L, method, hurst,
  methodfGn, sampling, subsampling){
  mod <- object@model

  if(missing(sampling)){
       sampling <- setSampling()
  }

  sim.Inputs <- simulate(mod, nsim, seed, xinit,
    true.parameter, space.discretized,
    increment.W, increment.L, method, hurst,
    methodfGn, sampling, subsampling)

  for(i in c(1:object@model@equation.number)){
    assign(object@Output@param@Input.var[[i]],
           get.zoo.data(sim.Inputs)[[i]])
  }
  assign(object@Output@param@time.var,
         sim.Inputs@sampling@grid[[1]])
  par <- unlist(true.parameter)
  nam <- names(par)
  for(i in c(1:length(nam))){
    assign(nam[i],par[i])
  }
  my.data<-eval(object@Output@formula[[1]])
  if(prod(object@Output@dimension)>1){
    for(i in c(2:prod(object@Output@dimension))){
      my.data<-cbind(my.data,
        eval(object@Output@formula[[i]]))
    }
  }
  names(my.data)<-object@Output@param@out.var

  data1 <- setData(my.data)
  object@data <- data1
  return(object)

}

# Method for Map
setMethod("simulate", "yuima.Integral",
          function(object, nsim=1, seed=NULL, xinit, true.parameter,
                   space.discretized=FALSE, increment.W=NULL, increment.L=NULL, method="euler",
                   hurst, methodfGn="WoodChan",
                   sampling, subsampling,
                   #Initial = 0, Terminal = 1, n = 100, delta,
                   #	grid, random = FALSE, sdelta=as.numeric(NULL),
                   #	sgrid=as.numeric(NULL), interpolation="none"
                   ...){
            res <- aux.simulatIntegral(object, nsim = nsim, seed = seed,
                                     xinit = xinit, true.parameter = true.parameter,
                                     space.discretized = space.discretized, increment.W = increment.W,
                                     increment.L = increment.L, method = method, hurst = hurst,
                                     methodfGn = methodfGn, sampling = sampling, subsampling = subsampling)

            return(res)
          }
)

aux.simulatIntegral <- function(object, nsim = nsim, seed = seed,
  xinit = xinit, true.parameter = true.parameter, space.discretized = space.discretized,
  increment.W = increment.W, increment.L = increment.L, method = method, hurst = hurst,
  methodfGn = methodfGn, sampling = sampling, subsampling = subsampling){

  if(missing(sampling)){
    sampling <- setSampling()
  }

  param <- unlist(true.parameter)
  info.par <- object@Integral@param.Integral
  info.var <- object@Integral@variable.Integral
  info.int <- object@Integral@Integrand

  mod1 <- object@model
  labmod.par <- mod1@parameter@all

  nm <- names(param)
  CondModPar <- nm%in%labmod.par
  ValModPar <- param[CondModPar]
  IntModPar <- param[!CondModPar]

  #Simulation Internal trajectories
  sim.Inputs <- simulate(mod1, nsim, seed, xinit, true.parameter,
    space.discretized, increment.W, increment.L, method, hurst,
    methodfGn, sampling)

  # Data of underlying SDE

  Data <- get.zoo.data(sim.Inputs)
  time <- index(sim.Inputs@data@original.data)
  my.env <- new.env()
  assign(info.var@var.time,time,envir=my.env)
  for(i in c(1:mod1@equation.number)){
    assign(mod1@solve.variable[i],
           as.numeric(Data[[i]]), envir = my.env)
  }
  df <- character(length=info.int@dimIntegrand[2])
  M.dX <- matrix(0,
     nrow = info.int@dimIntegrand[2],
     ncol = sim.Inputs@sampling@n)

  for(i in c(1:info.int@dimIntegrand[2])){
    df[i] <- paste0("diff(as.numeric(",info.var@var.dx[i],"))")
    M.dX[i,] <- eval(parse(text = df[i]), envir = my.env)
  }
  for(i in c(1:length(info.par@Integrandparam))){
    cond <- nm%in%info.par@Integrandparam[i]
    if(any(cond))
      assign(nm[cond],param[nm[cond]], envir = my.env)
  }

  #assign(info.var@var.time,time[-length(time)],envir=my.env)



#   matrInt <-matrix(0, nrow = info.int@dimIntegrand[1],
#                    ncol = info.int@dimIntegrand[2])
  res <- NULL
  PosInMatr <- matrix(c(1:(info.int@dimIntegrand[2]*info.int@dimIntegrand[1])),
    nrow = info.int@dimIntegrand[1], ncol = info.int@dimIntegrand[2])
  my.fun <- function(my.list, my.env, i){
    dum <- eval(my.list,envir = my.env)
    #if(length(dum)==1){
     # return(rep(dum,i))
    #}else{
      return(dum[1:i])
    #}
  }
 # res<-matrix(0,info.int@dimIntegrand[1],(length(time)-1))
 #  element <- matrix(0, nrow =info.int@dimIntegrand[1], ncol = 1)
 #
 #  for(i in c(1:(length(time)-1))){
 #    assign(info.var@upper.var,time[i+1],envir=my.env)
 #    Inter2 <- lapply(info.int@IntegrandList,
 #      FUN = my.fun, my.env = my.env, i = i)
 #    for(j in c(1:info.int@dimIntegrand[1])){
 #      element[j,] <- M.dX[,c(1:i)]%*%matrix(unlist(Inter2[PosInMatr[j,]]),
 #                                            ncol = info.int@dimIntegrand[2])
 #    }
 #    res[,i] <-  element
 #  }
  LengTime<-(length(time)-1)
  my.listenv<-as.list(c(1:LengTime))
  names(my.listenv)<- rep("i",LengTime)
  globalMyEnv <-new.env()
  globalMyEnv$time <- time
  globalMyEnv$my.env <- my.env
  element <- matrix(0, nrow =info.int@dimIntegrand[1], ncol = 1)
  my.listenv2<-lapply(X=my.listenv,
                      globalMyEnv=globalMyEnv,
                      FUN = function(X,globalMyEnv){
                        assign(info.var@upper.var,globalMyEnv$time[X+1],
                               envir= globalMyEnv$my.env)
                        assign(object@model@time.variable,globalMyEnv$time[c(1:X)],
                               envir= globalMyEnv$my.env)
                        Inter2 <- lapply(info.int@IntegrandList,
                                         FUN = my.fun, my.env = globalMyEnv$my.env,
                                         i = X)
                        for(j in c(1:info.int@dimIntegrand[1])){
                          element[j,] <- M.dX[,c(1:X)]%*%matrix(unlist(Inter2[PosInMatr[j,]]),
                                                                ncol = info.int@dimIntegrand[2])
                        }
                        list(element)
                      })
  res<-as.numeric(unlist(my.listenv2))
  res<-matrix(res,info.int@dimIntegrand[1],(length(time)-1))
  res <- cbind(0,res)
  rownames(res) <- object@Integral@variable.Integral@out.var
  my.data <- zoo(x = t(res), order.by = time)
  data1 <- setData(my.data)
  object@data <- data1
  if(missing(subsampling))
    return(object)
  subsampling(object, subsampling)
}
