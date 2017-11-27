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
  res <- aux.simulatPpr(object, nsim = nsim, seed = seed,
                        xinit = xinit, true.parameter = true.parameter,
                        space.discretized = space.discretized, increment.W = increment.W,
                        increment.L = increment.L, method = method, hurst = hurst,
                        methodfGn = methodfGn, sampling = sampling, subsampling = subsampling)

  return(res)

#   object@Kernel@param.Integral@allparam
#   simOzaki.aux(gFun=object@gFun@formula,a,cCoeff, Time, numJump)
}

setMethod("simulate", "yuima.Ppr",
          function(object, nsim=1, seed=NULL, xinit, true.parameter,
                   space.discretized=FALSE, increment.W=NULL, increment.L=NULL, method="euler",
                   hurst, methodfGn="WoodChan",
                   sampling, subsampling,
                   #Initial = 0, Terminal = 1, n = 100, delta,
                   #	grid, random = FALSE, sdelta=as.numeric(NULL),
                   #	sgrid=as.numeric(NULL), interpolation="none"
                   ...){
            res <- aux.simulatPpr(object, nsim = nsim, seed = seed,
                                     xinit = xinit, true.parameter = true.parameter,
                                     space.discretized = space.discretized, increment.W = increment.W,
                                     increment.L = increment.L, method = method, hurst = hurst,
                                     methodfGn = methodfGn, sampling = sampling, subsampling = subsampling)

            return(res)
          }
)

constHazIntPr <- function(g.Fun , Kern.Fun){
  numb.Int <- length(g.Fun)
  Int.Intens <- list()
  for(i in c(1:numb.Int)){
    dum.g <- as.character(g.Fun[i])
    dum.g <- paste("tail(",dum.g,", n=1L)")
    dum.Ker <- as.character(Kern.Fun@Integrand@IntegrandList[[i]])
    dif.dx <- paste("d",Kern.Fun@variable.Integral@var.dx, sep="")
    dum.Ker <- paste(dum.Ker,dif.dx, sep = "*")
    if(length(dum.Ker)>1){
      dum.Ker <- paste(dum.Ker,collapse = "+")
    }
    dum.Ker <- paste("(",dum.Ker,") * (")
    cond.Sup <- paste(Kern.Fun@variable.Integral@var.time, "<", Kern.Fun@variable.Integral@upper.var)
    dum.Ker <- paste(dum.Ker, cond.Sup, ")")
    dum.Ker <- paste0("sum(",dum.Ker,")")
    Int.Intens[[i]] <- parse(text = paste(dum.g, dum.Ker, sep = " + "))
  }
  res <- list(Intens = Int.Intens)
}

aux.simulatPpr<- function(object, nsim = nsim, seed = seed,
               xinit = xinit, true.parameter = true.parameter,
               space.discretized = space.discretized, increment.W = increment.W,
               increment.L = increment.L, method = method, hurst = hurst,
               methodfGn = methodfGn, sampling = sampling,
               subsampling = subsampling){
  ROLDVER<-!(is(object@model@measure$df,"yuima.law"))
  if(ROLDVER){
    object <- aux.simulatPprROldVersion(object, nsim = nsim, seed = seed,
                                          xinit = xinit, true.parameter = true.parameter,
                                          space.discretized = space.discretized, increment.W = increment.W,
                                          increment.L = increment.L, method = method, hurst = hurst,
                                          methodfGn = methodfGn, sampling = sampling,
                                          subsampling = subsampling)
  }else{
    object <- aux.simulatPprROldNew(object, nsim = nsim, seed = seed,
                                        xinit = xinit, true.parameter = true.parameter,
                                        space.discretized = space.discretized, increment.W = increment.W,
                                        increment.L = increment.L, method = method, hurst = hurst,
                                        methodfGn = methodfGn, sampling = sampling,
                                        subsampling = subsampling)
  }
  return(object)
}

aux.simulatPprROldNew<-function(object, nsim = nsim, seed = seed,
                          xinit = xinit, true.parameter = true.parameter,
                          space.discretized = space.discretized, increment.W = increment.W,
                          increment.L = increment.L, method = method, hurst = 0.5,
                          methodfGn = methodfGn, sampling = sampling,
                          subsampling = subsampling){

  # We need an expression for the evaluation of the hazard
  if(missing(hurst)){
    hurst<-0.5
  }
  samp <- sampling
  Model <- object@model
  gFun <- object@gFun
  Kern <- object@Kernel

  if(missing(xinit)){
    if(object@Ppr@RegressWithCount){

      yuima.warn("Counting Variables are also covariates.
                 In this case, the algorthim will be implemented
                 as soon as possible.")
      return(NULL)
    }
  }else{
    if(object@Ppr@RegressWithCount){
      yuima.warn("Counting Variables are also covariates.
                 In this case, the algorthim will be implemented
                 as soon as possible.")
      return(NULL)
    }
  }
  if(!object@Ppr@RegressWithCount && !object@Ppr@IntensWithCount){
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
      cod <-Model@solve.variable%in%object@Ppr@counting.var
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
    if(!object@Ppr@RegressWithCount && object@Ppr@IntensWithCount){
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
        Kern.Fun = object@Kernel)$Intens

      if(length(ExprHaz)==1){

        Time <- samp@Initial

        my.env <- new.env()
        for(i in c(1:length(object@Ppr@allparam))){
          assign(object@Ppr@allparam[i],
            as.numeric(true.parameter[object@Ppr@allparam[i]]),
            envir = my.env)
        }

        dimCov <- length(object@Ppr@covariates)
        if (dimCov>0){
          for(j in c(1:dimCov)){
          assign(object@Ppr@covariates[j],
                   as.numeric(simMod@data@original.data[1,object@Ppr@covariates[j]]),
                   envir = my.env)
          }
        }
        
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
            Lambda <- eval(ExprHaz[[1]], envir=my.env)
            return(Lambda)
          }else{
            return(NULL)
          }
        }

        IntensityProc <- 0
        #set.seed(1)
        dN<-0

        myhawkesP <- function(simMod, Kern,
                              samp, Model, my.env, ExprHaz,
                              Time){
          noExit<-TRUE
          const <- -log(runif(1))
          delta <- samp@delta
          grid <- samp@grid[[1]]
          while(const<delta){
            const <- -log(runif(1))
          }
          jumpT<-0
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
              dimCov <- length(object@Ppr@covariates)
                
              if (dimCov>0){
                  for(j in c(1:dimCov)){
                    assign(object@Ppr@covariates[j],
                           as.numeric(simMod@data@original.data[1:i,object@Ppr@covariates[j]]),
                           envir = my.env)
                  }
                }  
            

              
              cat("\n ", i, grid[i])
              }
            }
            if(i<dim(simMod@data@original.data)[1]){ 
              jumpT<-c(jumpT,grid[i])
              # if(i==7001){
              #   cat("\n",noExit)
              # }
              dN<-c(dN,1)
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

        prova1<-myhawkesP(simMod, Kern,
                          samp, Model, my.env, ExprHaz,
                          Time)

        Time<-prova1$jumpT
      #  return(Time)
      #   cost <- runif(1)
      #   while(-log(cost)<samp@delta)
      #   {cost <- runif(1)}
      #
      #   compErrHazR <- function(simMod, Kern,
      #                           samp, Model, my.env, ExprHaz,
      #                           cost, Time, dN, intensityProc){
      #     dummyLambda <- numeric(length=(samp@n+1))
      #     if(length(Kern@variable.Integral@var.dx)==1){
      #       MyPos <- sum(samp@grid[[1]]<=tail(Time,n=1L))
      #       assign(Kern@variable.Integral@var.time, Time, envir = my.env)
      #       i<- MyPos-1
      #       cond <- -log(cost)-sum(dummyLambda)*samp@delta
      #       while(cond>0 && i<=(samp@n) ){
      #         i<-i+1
      #         LastTime <- samp@grid[[1]][-1][i]
      #         assign(Model@time.variable, LastTime, envir = my.env)
      #         assign(paste0("d",Kern@variable.Integral@var.dx), dN, envir =my.env)
      #         dummyLambda[i] <- eval(ExprHaz[[1]], envir=my.env)
      #         cond <- -log(cost)-sum(dummyLambda)*samp@delta
      #
      #       }
      #
      #       if(i<=samp@n){
      #         NoExit<-TRUE
      #         Time <- c(Time,LastTime)
      #         dummydN <-  simMod@data@original.data[i+1,Kern@variable.Integral@var.dx]-simMod@data@original.data[i,Kern@variable.Integral@var.dx]
      #         dN <- c(dN,as.numeric(dummydN))
      #         intensityProc<-c(intensityProc,dummyLambda[MyPos:i])
      #       }else{
      #         NoExit<-FALSE
      #       }
      #       cat(sprintf("\n%.5f ", tail(Time,n=1L)))
      #       res <- list(Time = Time, dN = dN,
      #                   intensityProc = intensityProc, NoExit=NoExit)
      #       return(res)
      #     }else{
      #       return(NULL)
      #     }
      #   }
      #
      #   NoExit<-TRUE
      #   dN<-0
      #   intensityProc <- NULL
      #
      #   while(NoExit){
      #     ParRes <- compErrHazR(simMod, Kern,
      #                           samp, Model, my.env, ExprHaz,
      #                           cost, Time, dN, intensityProc)
      #     Time <- ParRes$Time
      #     dN <- ParRes$dN
      #     intensityProc <- ParRes$intensityProc
      #     NoExit <- ParRes$NoExit
      #     cost <- runif(1)
      #     while(-log(cost)<samp@delta)
      #     {cost <- runif(1)}
      #     # plot(y=intensityProc,x=samp@grid[[1]][-1][samp@grid[[1]][-1]<=tail(Time,n=1L)],type="l")
      #     # Time
      #     # NoExit
      #
      #   }
      # {
      # #   posLeft <- 1
      # #   posRight <- samp@n+1
      # #
      # #   posMid <- floor((posLeft+posRight)/2)
      # #   solveLeft <- -log(cost)
      # #   solveRight <- NULL
      # #
      # #   exit <- FALSE
      # #   prova <- NULL
      # #
      # #   globEx <- FALSE
      # #
      # #   while(tail(Time,n=1L)<(samp@Terminal-samp@delta) && !globEx){
      # #     while(!exit){
      # #       if((posMid-posLeft)==1){
      # #         posMid
      # #       }
      # #       oldprova <- prova
      # #       prova <- SolvePpr(posMid, posLeft, posRight, solveLeft, solveRight,
      # #                         cost, Kern, simMod, samp, Model, ExprHaz,
      # #                          my.env, Time, IntensityProc)
      # #       if(length(prova$left)==0){
      # #         globEx <- TRUE
      # #       }else{
      # #       if(prova$left){
      # #         posMid <- floor((prova$posLeft+prova$posRight)/2)
      # #         posLeft <- prova$posLeft
      # #         posRight <- prova$posRight
      # #         solveRight <- prova$solveRight
      # #         solveLeft <- prova$solveLeft
      # #       }else{
      # #         posMid <- floor((prova$posLeft+prova$posRight)/2)
      # #         posLeft <- prova$posLeft
      # #         posRight <- prova$posRight
      # #         solveRight <- prova$solveRight
      # #         solveLeft <- prova$solveLeft
      # #       }
      # #        exit<-prova$exit
      # #     }
      # #         if(globEx){
      # #           exit <- TRUE
      # #           oldprova -> prova
      # #         }
      # #
      # #
      # #     }
      # #
      # #     cost <- runif(1)
      # #     while(-log(cost)<samp@delta)
      # #     {cost <- runif(1)}
      # #     exit<- FALSE
      # #     posRight <- samp@n+1
      # #     solveLeft <- -log(cost)
      # #      posMid <- floor((posLeft+posRight)/2)
      # #     # if(length(prova$VeryExit)!=0){
      # #     #   if(prova$VeryExit)
      # #     #     globEx <- TRUE
      # #     # }
      # #     if(!globEx)
      # #       Time <- prova$Time
      # #     IntensityProc <- prova$IntensityProc
      # #     #cat(tail(prova$Time,n=1L))
      # #     #cat(sprintf("\n%.5f ", posMid))
      # #     cat(sprintf("\n%.5f ", tail(Time,n=1L)))
      # # }
      # }
        cond <- samp@grid[[1]][-1] %in% Time[-1]
        countVar <- Model@solve.variable %in%  object@Ppr@counting.var
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

      }
    }
  }
  return(NULL)
}

# SolvePpr <- function(posMid, posLeft, posRight, solveLeft = NULL, solveRight = NULL,
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


# SolvePpr <- function(TopposInGridIn, OldTimePoint, solveLambdaInOld,
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




aux.simulatPprROldVersion <- function(object, nsim = nsim, seed = seed,
                                      xinit = xinit, true.parameter = true.parameter,
                                      space.discretized = space.discretized, increment.W = increment.W,
                                      increment.L = increment.L, method = method, hurst = hurst,
                                      methodfGn = methodfGn, sampling = sampling,
                                      subsampling = subsampling){
  Time <- sampling@Terminal
  numbVardx <- length(object@Ppr@var.dx)
  numbCountVar <- length(object@Ppr@counting.var)
  U <- runif(numbCountVar)

  my.env<- new.env()

  true.parameter <- unlist(true.parameter)

  if(!all(names(true.parameter)==object@Ppr@allparam)){
    yuima.stop("true.parameters mismatch the model parameters")
  }
  for(i in c(1:length(object@Ppr@allparam))){
    assign(object@Ppr@allparam[i],true.parameter[object@Ppr@allparam[i]], envir = my.env)
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
    # if(length(object@Ppr@covariates)!=0)
    #    Condcovariate <- TRUE
    cond <- !(object@model@solve.variable %in% object@Ppr@counting.var)
    if(any(cond)){
      Condcovariate <- TRUE
    }
    dimMd <- length(object@model@solve.variable)
    dumMod <- setModel(drift = rep("0",dimMd),
                       diffusion = matrix("0",dimMd,1),
                       jump.coeff = diag("1",dimMd,dimMd),
                       measure = object@Ppr@Info.measure$measure,
                       measure.type = object@Ppr@Info.measure$type,
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
      for(i in c(1:length(object@Ppr@covariates))){
        assign(object@Ppr@covariates[i],
               eval(yuima@xinit[yuima@solve.variable==object@Ppr@covariates[i]],
                    envir = env), envir = env)
      }
      if(u!=0){
        # Mat<-matrix(0,length(yuima@solve.variable),length(env$t)+1)
        # for(i in c(1:length(yuima@solve.variable))){
        #   Mat[i,1] = eval(yuima@xinit[i],envir = env)
        # }
        Linc <- env$L[,c(1:(length(part)-1))]
        # Linc[yuima@solve.variable!=object@Ppr@covariates,]<-matrix(0,
        #   sum(yuima@solve.variable!=object@Ppr@covariates), dim(Linc)[2])
        Linc[yuima@solve.variable!=object@Ppr@covariates,] <- 0
        DumUnderlMod <- simulate(yuima, true.parameter = true.parameter,
                                 increment.L = env$L[,c(1:(length(part)-1))],
                                 sampling = setSampling(Terminal = u, n= (length(part)-1)))


        for(i in c(1:length(object@Ppr@covariates))){
          VariableDum <- DumUnderlMod@data@original.data[,yuima@solve.variable==object@Ppr@covariates[i]]
          assign(object@Ppr@covariates[i], as.numeric(VariableDum), envir = env)
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
      LevIncr[object@Ppr@counting.var,]<-0
      yuima<- object@model
      xinit<- numeric(length(object@Ppr@covariates))
      names(xinit)<- object@Ppr@covariates
      for(i in c(1:length(object@Ppr@covariates))){
        xinit[i] <- env[[object@Ppr@covariates[i]]]
      }

      xinitCount <- numeric(length(object@Ppr@counting.var))
      names(xinitCount) <- object@Ppr@counting.var
      for(i in c(1:length(xinitCount))){
        xinitCount[i] <- tail(env[[object@Ppr@counting.var[i]]],n = 1)
      }
      xinit <- c(xinit,xinitCount)
      if(part[length(part)]-part[1]!=0){
        DumVarCov  <- simulate(yuima,
                               true.parameter = true.parameter,
                               increment.L = LevIncr,
                               sampling =  setSampling(Terminal = (part[length(part)]-part[1]),
                                                       n = dim(LevIncr)[2]),
                               xinit=xinit[yuima@solve.variable])
        for(i in c(1:length(object@Ppr@covariates))){
          VariableDum <- DumVarCov@data@original.data[,yuima@solve.variable==object@Ppr@covariates[i]]
          assign(object@Ppr@covariates[i], as.numeric(VariableDum), envir = env)
        }
      }else{
        for(i in c(1:length(object@Ppr@covariates))){
          VariableDum <- xinit[yuima@solve.variable==object@Ppr@covariates[i]]
          assign(object@Ppr@covariates[i], as.numeric(VariableDum), envir = env)
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
  names(u) <- object@Ppr@counting.var
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

  rownames(JUMP)<- object@Ppr@counting.var
  pos<-sum(sampling@grid[[1]][-1]<=t_1)
  t_1 <- sampling@grid[[1]][-1][pos]
  recordTime<-c(0,t_1)
  pos0<-0

  JUMP[namesContVarJump, pos] <- L[namesContVarJump, pos]
  ntot <- sampling@n
  dL <- L
  dL[object@Ppr@counting.var,c((pos0+1):pos)]<-JUMP[object@Ppr@counting.var,c((pos0+1):pos)]

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
      dL[object@Ppr@counting.var,c((pos0+1):pos)]<-JUMP[object@Ppr@counting.var,c((pos0+1):pos)]
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
  N.VarPPr<-length(solve.variable)

  dummy.mod <- setModel(drift=rep("0",N.VarPPr),
                        diffusion = NULL, jump.coeff = diag(rep("1",N.VarPPr)),
                        measure = object@Ppr@Info.measure$measure,
                        measure.type = object@Ppr@Info.measure$type,
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
  dummyObject@variable.Integral@out.var <-object@Ppr@additional.info
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
  colnames(Alldata)<-c(solve.variable,object@Ppr@additional.info)
  # for(i in c(1:N.VarPPr)){
  #   assign(solve.variable[i],interRes@data@original.data[,i],envir=my.env)
  # }
  # dummy<-NULL
  # for(t in c(1:length(object@Ppr@additional.info))){
  #   dummy <-eval(object@gFun)
  #   assign(object@Ppr@additional.info[[]])
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
