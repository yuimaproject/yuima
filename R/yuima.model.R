# setMethod("initialize", "model.parameter",
#           function(.Object,
#                    all,
#                    common,
#                    diffusion,
#                    drift,
#                    jump,
#                    measure,
#                    xinit){
#             .Object@all <- all
#             .Object@common <- common
#             .Object@diffusion <- diffusion
#             .Object@drift <- drift
#             .Object@jump <- jump
#             .Object@measure <- measure
#             .Object@xinit <- xinit
#             return(.Object)
#           })

setMethod("initialize", "model.parameter",
          function(.Object,
                   all = character(),
                   common = character(),
                   diffusion = character(),
                   drift = character(),
                   jump = character(),
                   measure = character(),
                   xinit = character()){
            .Object@all <- all
            .Object@common <- common
            .Object@diffusion <- diffusion
            .Object@drift <- drift
            .Object@jump <- jump
            .Object@measure <- measure
            .Object@xinit <- xinit
            return(.Object)
          })

# setMethod("initialize", "yuima.model",
#           function(.Object,
#                    drift ,
#                    diffusion,
# 				           hurst,
#                    jump.coeff,
#                    measure,
#                    measure.type,
#                    parameter,
#                    state.variable,
#                    jump.variable,
#                    time.variable,
#                    noise.number,
#                    equation.number,
#                    dimension,
#                    solve.variable,
#                    xinit,
#                    J.flag){
#             .Object@drift <- drift
#             .Object@diffusion <- diffusion
# 			.Object@hurst <- hurst
#             .Object@jump.coeff <- jump.coeff
#             .Object@measure <- measure
#             .Object@measure.type <- measure.type
#             .Object@parameter <- parameter
#             .Object@state.variable <- state.variable
#             .Object@jump.variable <- jump.variable
#             .Object@time.variable <- time.variable
#             .Object@noise.number <- noise.number
#             .Object@equation.number <- equation.number
#             .Object@dimension <- dimension
#             .Object@solve.variable <- solve.variable
#             .Object@xinit <- xinit
#             .Object@J.flag <- J.flag
#             return(.Object)
#           })

# 23/11 we need to provide the default values for the yuima.model object class
# in order to construct a new class that inherits from yuima.model class

setMethod("initialize", "yuima.model",
          function(.Object,
                   drift = expression() ,
                   diffusion = list() ,
                   hurst = 0.5,
                   jump.coeff = list(),
#jump.coeff = expression(),
                   measure=list(),
                   measure.type=character(),
                   parameter = new("model.parameter"),
                   state.variable = "x",
                   jump.variable = "z",
                   time.variable = "t",
                   noise.number = numeric(),
                   equation.number = numeric(),
                   dimension = numeric(),
                   solve.variable = character(),
                   xinit = expression(),
                   J.flag = logical()){
            .Object@drift <- drift
            .Object@diffusion <- diffusion
            .Object@hurst <- hurst
            .Object@jump.coeff <- jump.coeff
            .Object@measure <- measure
            .Object@measure.type <- measure.type
            .Object@parameter <- parameter
            .Object@state.variable <- state.variable
            .Object@jump.variable <- jump.variable
            .Object@time.variable <- time.variable
            .Object@noise.number <- noise.number
            .Object@equation.number <- equation.number
            .Object@dimension <- dimension
            .Object@solve.variable <- solve.variable
            .Object@xinit <- xinit
            .Object@J.flag <- J.flag
            return(.Object)
          })


## setModel
## setter of class 'yuima.model'
## set yuima model from SDE
setModel <- function(drift=NULL,
                     diffusion=NULL,
                     hurst=0.5,
                     jump.coeff=NULL,
                     measure=list(),
                     measure.type=character(),
                     state.variable="x",
                     jump.variable="z",
                     time.variable="t",
                     solve.variable,
                     xinit=NULL){
  ## we need a temp env for simplifications
  mylengdumMeas<-length(measure)
  if(mylengdumMeas>0){
    for(i in c(1:mylengdumMeas)){
      if(is(measure[[i]],"yuima.law")){
        res<- aux.setModelLaw(drift,diffusion,
          hurst, jump.coeff, measure, measure.type,
          state.variable, jump.variable, time.variable,
          solve.variable, xinit, posyuimalaw=i)
        res@measure[[i]]<-measure[[i]]
        return(res)
      }
    }
  }

  yuimaENV <- new.env()
  ##::measure and jump term #####################################

  ##::initialize objects ########
  MEASURE <- list()

  ##::end initialize objects ########

  ##::make type function list ####
  CPlist <- c("dnorm", "dgamma", "dexp", "dconst")
  codelist <- c("rIG", "rNIG", "rgamma", "rbgamma", "rvgamma", "rstable","rpts","rnts","yuima.law")
  ## added "yuima.law",dconst LM (2024/09/10)
  ## added "rpts" and "rnts" by YU (2016/10/4)
  ##::end make type function list ####

  ## Multivariate YUIMA model


  if(!is.null(jump.coeff)){
    if(is.matrix(jump.coeff)){

      if(dim(jump.coeff)[2]!=1){
        intensity <- NULL
        #if(is.null(names(measure)) || names(measure)=="df"){
        if(is.null(names(measure)) || all(names(measure)%in%"df")){
          names(measure) <- "df"
        }

        df <- as.list(measure[["df"]])
        if(any(measure.type=="CP")){
          intensity <- measure[["intensity"]]
        }
        my.cond <- TRUE
        tmp <- regexpr("\\(", measure$df)[1]
        measurefunc <- substring(measure$df, 1, tmp-1)
        if(!is.na(match(measurefunc, codelist))){
          my.cond <- FALSE
        }
        if(my.cond){
          res <- setMultiModel(drift = drift, diffusion = diffusion,
                               hurst = hurst, jump.coeff = jump.coeff,
                               intensity = intensity, df = df,
                               measure.type = measure.type, state.variable = state.variable,
                               jump.variable = jump.variable, time.variable = time.variable,
                               solve.variable = solve.variable, xinit= xinit)
          return(res)
        }
      }
    }
  }


  if(!length(measure.type)){
    if( length(jump.coeff) || length(measure) ){
      yuima.warn("measure type does not match with jump term.")
      return(NULL)
    }
    jump.variable <- character()
    measure.par <- character()
  }else{
    if( !length(jump.coeff) || !length(measure) ){
      yuima.warn("measure type isn't matched with jump term.")
      return(NULL)
   # }else
      #       if(length(jump.coeff)!=1){
      #        yuima.warn("multi dimentional jump term is not supported yet.")
      #
      #         return(NULL)
      #     }

    } else if(measure.type=="CP"){ ##::CP
        #        if(length(measure)!=2){
        # yuima.warn(paste("length of measure must be two on type", measure.type, "."))
        # return(NULL)
        #}
        if(!is.list(measure)){
          measure <- list(intensity=measure[1], df=measure[2],dimension=measure[3])
        } else {
        #if(length(measure[[1]])!=1 || length(measure[[2]])!=1){
        #   yuima.warn("multi dimentional jump term is not supported yet.")
        #   return(NULL)
        # }
          ##::naming measure list ########
          tmpc <- names(measure)
          if(is.null(tmpc)){
            names(measure) <- c("intensity", "df","dimension")
          }else{
            whichint <- match("intensity", tmpc)
            whichdf <- match("df", tmpc)
            if(!is.na(whichint)){
              if(names(measure)[-whichint]=="" || names(measure)[-whichint]=="df"){
                names(measure)[-whichint] <- "df"
              }else{
                yuima.warn("names of measure are incorrect.")
                return(NULL)
              }
            }else if(!is.na(whichdf)){
              if(names(measure)[-whichdf]=="" || names(measure)[-whichdf]=="intensity"){
                names(measure)[-whichdf] <- "intensity"
              }else{
                yuima.warn("names of measure are incorrect.")
                return(NULL)
              }
            }else{
              yuima.warn("names of measure are incorrect.")
              return(NULL)
            }
          }
          ##::end naming measure list ########
        }

        ##::check df name ####################
        tmp <- regexpr("\\(", measure$df)[1]
        measurefunc <- substring(measure$df, 1, tmp-1)
        if(!is.na(match(measurefunc, codelist))){
          yuima.warn(paste("distribution function", measurefunc, "should be defined as type code."))
          return(NULL)
        }#else if(is.na(match(measurefunc, CPlist))){
        # warning(paste("\ndistribution function", measurefunc, "is not officialy supported as type CP.\n"))
        #}
        MEASURE$df$func <- eval(parse(text=measurefunc)) #LM 15/05/2017
        MEASURE$df$expr <- parse(text=measure$df)
        MEASURE$intensity <- parse(text=measure$intensity)

        measure.par <- unique( c( all.vars(MEASURE$intensity), all.vars(MEASURE$df$expr) ) )
        ##measure.par$intensity <- unique(all.vars(MEASURE$intensity))
        ##::end check df name ####################
        ##::end CP

      } else if(measure.type=="code"){ ##::code
        if(length(measure)!=1){
          yuima.warn(paste("length of measure must be one on type", measure.type, "."))
          return(NULL)
        }

        if(!is.list(measure)){
          measure <- list(df=measure)
        }else{
          if(length(measure[[1]])!=1){
            yuima.warn("multi dimentional jump term is not supported yet.")
            return(NULL)
          }
          ##::naming measure list #############
          if(is.null(names(measure)) || names(measure)=="df"){
            names(measure) <- "df"
          }else{
            yuima.warn("name of measure is incorrect.")
            return(NULL)
          }
          ##::end naming measure list #############
        }

        ##::check df name ####################
        tmp <- regexpr("\\(", measure$df)[1]
        measurefunc <- substring(measure$df, 1, tmp-1)
        if(!is.na(match(measurefunc, CPlist))){
          yuima.warn(paste("\ndistribution function", measurefunc, "should be defined as type CP."))
          return(NULL)
        }else if(is.na(match(measurefunc, codelist))){
          warning(paste("\ndistribution function", measurefunc, "is not officialy supported as type code.\n"))
        }
        ##MEASURE$df$func <- eval(parse(text=measurefunc))
        MEASURE$df$expr <- parse(text=measure$df)

        measure.par <- unique(all.vars(MEASURE$df$expr))
        ##::end check df name ####################
        ##::end code
      }else if(measure.type=="density"){ ##::density
        if(length(measure)!=1){
          yuima.warn(paste("length of measure must be one on type", measure.type, "."))
          return(NULL)
        }

        if(!is.list(measure)){
          measure <- list(df=measure)
        }else{
          if(length(measure[[1]])!=1){
            yuima.warn("multi dimentional jump term is not supported yet.")
            return(NULL)
          }

          ##::naming measure list #############
          if(is.null(names(measure))){
            names(measure) <- "df"
          }else if(names(measure)!="density" && names(measure)!="df"){
            yuima.warn("name of measure is incorrect.")
            return(NULL)
          }
          ##::end naming measure list #############
        }

        ##::check df name ####################
        tmp <- regexpr("\\(", measure[[names(measure)]])[1]
        measurefunc <- substring(measure[[names(measure)]], 1, tmp-1)
        if(!is.na(match(measurefunc, CPlist))){
          yuima.warn(paste("distribution function", measurefunc, "should be defined as type CP."))
          return(NULL)
        }else if(!is.na(match(measurefunc, codelist))){
          yuima.warn(paste("distribution function", measurefunc, "should be defined as type code."))
          return(NULL)
        }
        MEASURE[[names(measure)]]$func <- eval(parse(text=measurefunc))
        MEASURE[[names(measure)]]$expr <- parse(text=measure[[names(measure)]])

        measure.par <- unique(all.vars(MEASURE[[names(measure)]]$expr))
        ##::end check df name ####################
        ##::end density
      }else{ ##::else
        yuima.warn(paste("measure type", measure.type, "isn't supported."))
        return(NULL)
      }
    n.eqn3 <- 1
    n.jump <- 1
  }




  ##::end measure and jump term #####################################

  ##:: check for errors and reform values
  if(any(time.variable %in% state.variable)){
    yuima.warn("time and state(s) variable must be different.")
    return(NULL)
  }
  if(is.null(dim(drift))){ # this is a vector
    n.eqn1 <- length(drift)
    n.drf <- 1
  }else{ # it is a matrix
    n.eqn1 <- dim(drift)[1]
    n.drf <- dim(drift)[2]
  }

  if(is.null(dim(diffusion))){ # this is a vector
    n.eqn2 <- length(diffusion)
    n.noise <- 1
  }else{ # it is a matrix
    n.eqn2 <- dim(diffusion)[1]
    n.noise <- dim(diffusion)[2]
  }

  if(is.null(diffusion)){
    diffusion <- rep("0", n.eqn1)
    n.eqn2 <- n.eqn1
    n.noise <- 1
  }

  ## TBC
  n.eqn3 <- n.eqn1

  if(!length(measure)){
    n.eqn3 <- n.eqn1
  }

  if(n.eqn1 != n.eqn2 || n.eqn1 != n.eqn3){
    yuima.warn("Malformed model, number of equations in the drift and diffusion do not match.")
    return(NULL)
  }
  n.eqn <- n.eqn1

  if(is.null(xinit)){
    # xinit <- numeric(n.eqn)
    xinit <- character(n.eqn)
  }else if(length(xinit) != n.eqn){
    if(length(xinit)==1){
      xinit <- rep(xinit, n.eqn)
    }else{
      yuima.warn("Dimension of xinit variables missmatch.")
      return(NULL)
    }
  }

  if(missing(solve.variable)){
    yuima.warn("Solution variable (lhs) not specified. Trying to use state variables.")
    solve.variable <- state.variable
  }
  if(n.eqn != length(solve.variable)){
    yuima.warn("Malformed model, number of solution variables (lhs) do no match number of equations (rhs).")
    return(NULL)
  }

  loc.drift <- matrix(drift, n.eqn, n.drf)
  loc.diffusion <- matrix(diffusion, n.eqn, n.noise)
  # Modification starting point 6/11
  loc.xinit<-matrix(xinit,n.eqn,n.drf)

  ##:: allocate vectors
  DRIFT <- vector(n.eqn, mode="expression")
  DIFFUSION <- vector(n.eqn, mode="list")
  # Modification starting point 6/11
  XINIT<-vector(n.eqn, mode = "expression")

  ##:: function to make expression from drift characters
  pre.proc <- function(x){
    for(i in 1:length(x)){
      if(length(parse(text=x[i]))==0){
        x[i] <- "0"
      }
    }
    parse(text=paste(sprintf("(%s)", x), collapse="+"))
  }
  ##22/11:: function to simplify expression in drift, diffusion, jump and xinit characters
  yuima.Simplifyobj<-function(x){
    dummy<-yuima.Simplify(x, yuima.env=yuimaENV)
    dummy1<-yuima.Simplify(dummy, yuima.env=yuimaENV)
    dummy2<-as.character(dummy1)
    res<-parse(text=paste0("(",dummy2,")",collapse=NULL))
    return(res)
  }


  ##:: make expressions of drifts and diffusions and jump
  for(i in 1:n.eqn){
    DRIFT[i] <- pre.proc(loc.drift[i,])
    # 22/11 Simplify expressions
    DRIFT[i] <- yuima.Simplifyobj(DRIFT[i])
    # Modification starting point 6/11
    XINIT[i]<-pre.proc(loc.xinit[i, ])
    XINIT[i]<- yuima.Simplifyobj(XINIT[i])
    for(j in 1:n.noise){
      expr <- parse(text=loc.diffusion[i,j])
      if(length(expr)==0){
        expr <- expression(0)  # expr must have something
      }
#       DIFFUSION[[i]][j] <- expr
      #22/11
      DIFFUSION[[i]][j] <- yuima.Simplifyobj(expr)
    }
#22/11

#if (length(JUMP)>0){
#    JUMP[i] <- parse(text=jump.coeff[i])
#    JUMP[i] <- yuima.Simplifyobj(JUMP[i])
#}

  }



#print(length(jump.coeff))
#if (length(jump.coeff)==0){
#    JUMP <- list(parse(text=jump.coeff))
#}else{
#    #    JUMP <- vector(n.eqn, mode="expression")
#    JUMP <- vector(n.eqn, mode="list")
#}

if(length(jump.coeff)==0){
    JUMP <- list()
} else {
    if(length(jump.coeff)==1 & !is.matrix(jump.coeff)){ # is a scalar
        expr <- parse(text=jump.coeff)
        if(length(expr)==0){
            expr <- expression(0)  # expr must have something
        }
        JUMP <- list(yuima.Simplifyobj(expr))
    } else { # must be matrix, n.col = dimension of Levy noise
        jump.coeff <- as.matrix(jump.coeff)
        c.j <- ncol(jump.coeff)
        r.j <- nrow(jump.coeff)
        #print(c.j)
        #print(r.j)
        #print(jump.coeff)
        JUMP <- vector(r.j, mode="list")
        for(i in 1:r.j){
            for(j in 1:c.j){
                #cat(sprintf("\ni=%d,j=%d\n",i,j))
                expr <- parse(text=jump.coeff[i,j])
                if(length(expr)==0){
                    expr <- expression(0)  # expr must have something
                }
                JUMP[[i]][j] <- yuima.Simplifyobj(expr)
            }
        }
    }
}
#print(str(JUMP))

 #

  ##:: get parameters in drift expression
  drift.par <- unique(all.vars(DRIFT))
  # Modification starting point 6/11
  xinit.par <- unique(all.vars(XINIT))

  drift.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), drift.par)))
  if(length(drift.idx)>0){
    drift.par <- drift.par[-drift.idx]
  }

  ##:: get parameters in diffusion expression
  diff.par <- unique(unlist(lapply(DIFFUSION, all.vars)))
  diff.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), diff.par)))
  if(length(diff.idx)>0){
    diff.par <- diff.par[-diff.idx]
  }

  ##:: get parameters in jump expression
  J.flag <- FALSE
  #  jump.par <- unique(all.vars(JUMP))
  jump.par <- unlist(lapply(JUMP,all.vars))
  if(is.null(jump.par))
   jump.par <- character()
  if(length(na.omit(match(jump.par, jump.variable)))){
    J.flag <- TRUE
  }
  jump.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), jump.par)))
  if(length(jump.idx)>0){
    jump.par <- jump.par[-jump.idx]
  }

  ##:: get parameters in measure expression
  measure.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), measure.par)))
  if(length(measure.idx)>0){
    measure.par <- measure.par[-measure.idx]
  }

  ##:: order parameters for 'yuima.pars'
  ##id1 <- which(diff.par %in% drift.par)
  ##id2 <- which(drift.par %in% diff.par)
  ##common <- unique(c(diff.par[id1], drift.par[id2]))
  common <- c(drift.par, diff.par)
  common <- common[duplicated(common)]

  common1<-common
  # modification 06/11 common1 contains only
  # parameters that appear in both drift and diffusion terms.

  # Modification 06/11 common contains only parameters that appear
  # in drift, diff, Jump and xinit
  if (length(xinit)) {
    common <- c(common, xinit.par)
    common <- common[duplicated(common)]
    common <- c(common, xinit.par)
    common <- common[duplicated(common)]
  }


  if(length(measure)){
    common <- c(common, jump.par)
    common <- common[duplicated(common)]
    common <- c(common, measure.par)
    common <- common[duplicated(common)]
  }
  #   all.par <- unique(c(drift.par, diff.par, jump.par, measure.par))
  all.par <- unique(c(drift.par, diff.par, jump.par, measure.par, xinit.par))

  ##:: instanciate class
  tmppar <- new("model.parameter",
                all= all.par,
                #                 common= common,
                common= common1,
                diffusion= diff.par,
                drift= drift.par,
                jump= jump.par,
                measure= measure.par,
                xinit=xinit.par)
  tmp <- new("yuima.model",
             drift= DRIFT,
             diffusion= DIFFUSION,
             hurst=as.numeric(hurst),
             jump.coeff=JUMP,
             measure= MEASURE,
             measure.type= measure.type,
             parameter= tmppar,
             state.variable= state.variable,
             jump.variable= jump.variable,
             time.variable= time.variable,
             noise.number= n.noise,
             equation.number= n.eqn,
             dimension= c(
               length(tmppar@all),
               length(tmppar@common),
               length(tmppar@diffusion),
               length(tmppar@drift),
               length(tmppar@jump),
               length(tmppar@measure)
             ),
             solve.variable= solve.variable,
             xinit= XINIT,
             J.flag <- J.flag)
  return(tmp)
}

aux.setModelLaw <- function(drift,diffusion,
                      hurst, jump.coeff, measure, measure.type,
                      state.variable, jump.variable, time.variable,
                      solve.variable, xinit, posyuimalaw){

  dummyMeasure <- paste0(c("yuima.law(",
                    paste0(measure[[posyuimalaw]]@param.measure,collapse=", ")
                    ,")"), collapse="")
  auxmeasure <- measure
  auxmeasure[[posyuimalaw]]<-dummyMeasure
  names(auxmeasure[posyuimalaw]) <- "df"
  setModel(drift = drift,diffusion = diffusion,
           hurst = hurst, jump.coeff = jump.coeff, measure = auxmeasure,
           measure.type = measure.type,
           state.variable = state.variable,
           jump.variable = jump.variable, time.variable,
           solve.variable, xinit)
}

# yuima.model rbind

# setGeneric("rbind.yuima",
#            function(x, ...)
#              standardGeneric("rbind.yuima")
# )

# setMethod("cbind.yuima", signature(x="yuima"),
#           function(x, ...){
#             ##:: init
#             y.list <- list(x, ...)
#             y.num <- length(y.list)
#
#             ##:: bind yuima.data in yuima
#
#             ##:: return result
#             return(NULL)
#           }
# )

# setMethod("rbind.yuima", signature(x="yuima.model"),
#           function(x, ...){
#             y.list <- list(x, ...)
#             y.num <- length(y.list)
#             res <- aux.rbind.model(y.list,y.num)
#             return(res)
#           }
# )

rbind.yuima.model <- function(x, ...){

  y.list <- list(x, ...)
#  y.list1 <- lapply(y.list, FUN = only.yuima.model)
  y.num <- length(y.list)
  new.list <- list()
  for(i in (1:y.num)){
    if(is(y.list[[i]],"yuima.model"))
          new.list[i] <- y.list[[i]]
  }
  new.y.num <- length(new.list)
  res <- aux.rbind.model(y.list = new.list,
        y.num = new.y.num, mycall = y.list)
  return(res)
}

aux.rbind.model<-function(y.list,y.num, mycall=list()){
  lapply(y.list, FUN = check.yuima.model)
  check.lev <- lapply(y.list, FUN = check.yuima.levy)
  check.lev <- unlist(check.lev)
  drift <- lapply(y.list, FUN = extract.model, type = "drift")
  diffusion <- lapply(y.list, FUN = extract.model, type = "diffusion")
  solve.variable <- lapply(y.list, FUN = extract.model, type = "solve.variable")
  state.variable <- lapply(y.list, FUN = extract.model, type = "state.variable")
  xinit <- lapply(y.list, FUN = extract.model, type = "xinit")
  noise.number <- lapply(y.list, FUN = extract.model, type = "noise.number")
  equation.number <- lapply(y.list, FUN = extract.model, type = "equation.number")
  #Until Here only diffusion process
  drift <- lapply(drift, FUN = ExpToString)
  drift <- unlist(drift)
  # drift
  nrow.diff <- sum(unlist(equation.number))
  ncol.diff <- sum(unlist(noise.number))
  matr.diff <- matrix("0", nrow = nrow.diff, ncol = ncol.diff)
  extrinf <- 1
  extrsup <- noise.number[[1]]
  j <- 1
  cond.eq <- equation.number[[1]]
  cond.eq1 <- 0
  for(i in c(1:nrow.diff)){
     if(i <= cond.eq){
        dum <- ExpToString(diffusion[[j]][[i-cond.eq1]])
        matr.diff[i,extrinf:extrsup] <- dum
        if(i == equation.number[[j]]){
          extrinf <- extrsup+1
          j <- j+1
          if(j <= nrow.diff){
            extrsup <- extrsup +  equation.number[[j]]
            cond.eq1 <- i
            cond.eq <- cond.eq +  equation.number[[j]]
          }
        }
     }
  }
  solve.variable <- lapply(solve.variable, FUN = ExpToString, cond = FALSE)
  solve.variable <- unlist(solve.variable)
  state.variable <- lapply(state.variable, FUN = ExpToString, cond = FALSE)
  state.variable <- unlist(state.variable)
  xinit <- lapply(xinit, FUN = ExpToString, cond = FALSE)
  xinit <- unlist(xinit)
  if(!any(check.lev)){
    mod <- setModel(drift = drift, diffusion = matr.diff,
      solve.variable = solve.variable, state.variable = state.variable,
      xinit = xinit)
  }else{
    MultiLevy <- y.list[check.lev]
    jump.coeff <- lapply(MultiLevy,
      FUN = extract.model, type = "jump.coeff")
    ncol.jump <- lapply(jump.coeff, FUN = numb.jump)
    dum.ncolj <- unlist(ncol.jump)
    ncol.jump <- sum(unlist(dum.ncolj))
    jump.coeff <- lapply(y.list,
      FUN = extract.model, type = "jump.coeff")
    #ncol.jump1 <- lapply(jump.coeff, FUN = numb.jump)
    matr.jump <- matrix("0",nrow = nrow.diff,
      ncol = ncol.jump)
    j <- 1
    h <- 0
    cond.eqa <- equation.number[[j]]
    cond.eqb <- 0
    extrinf <- 1
    extrsup <- 1
    if(check.lev[j])
      extrsup <- dum.ncolj[j]
    else{
      h <- h+1
    }
    for(i in c(1:nrow.diff)){
      if(i <= cond.eqa){
        if(check.lev[j]){
          dum <- ExpToString(jump.coeff[[j]][[i-cond.eqb]])
          matr.jump[i, extrinf:extrsup] <- dum
        }else{
#          matr.jump[i,] <- matr.jump[i,]
        }
        if(i == cond.eqa){
          cond.eqb <- i
          j <- j+1
          if(j<=length(equation.number))
                cond.eqa <- cond.eqa + equation.number[[j]]
          if(check.lev[j-1]){
            extrinf <- extrsup + 1
            extrsup <- extrsup + dum.ncolj[j-h]
          }else{
            extrinf <- extrinf
            extrsup <- extrsup
            h <- h+1
          }
        }
      }
    }

    # mod <- matr.jump
#     measure <- lapply(y.list,
#        FUN = extract.model, type = "measure")
#     measure
    df <- NULL
  if("df" %in% names(mycall))
      df <- mycall$df
    measure.type <- NULL
  if("measure.type" %in% names(mycall))
      measure.type <- mycall$measure.type
    intensity <-NULL
    if("intensity" %in% names(mycall))
      intensity <- mycall$intensity
    time.variable <- "t"
    if("time.variable" %in% names(mycall))
      time.variable <- mycall$time.variable
  mod <- setMultiModel(drift=drift, diffusion = matr.diff,
    jump.coeff =  matr.jump, solve.variable = solve.variable,
    xinit = xinit, time.variable = time.variable, df= df,
    intensity = intensity, measure.type = measure.type)

  }
  return(mod)
}
# only.yuima.model<- function(y.list){
#   if(is(y.list,"yuima.model")){
#     return(y.list)
#   }else{
#     NULL
#   }
# }
numb.jump <- function(x){length(x[[1]])}

check.yuima.levy <- function(x){
  Levy <- FALSE
  if(length(x@measure.type)>0){
    if(!is(x, "yuima.model")){
      yuima.stop("the Levy model have to belong to the yuima.multimodel class")
    }
    Levy <- TRUE
  }
  return(Levy)
}

ExpToString <- function(x, cond = TRUE){
  dum <- unlist(strsplit(toString(x),split=", "))
  if(cond)
    dum <- substr(dum, 2, nchar(dum)-1)
  return(dum)
}

extract.model <- function(x, type = "drift"){
  res<- slot(x,type)
  return(res)
}

check.yuima.model <- function(x){
  if(is.CARMA(x)){
    yuima.warn("The cbind for CARMA will be implemented as soon as possible")
    return(NULL)
  }
  if(is.COGARCH(x)){
    yuima.warn("The cbind for COGARCH will be implemented as soon as possible")
    return(NULL)
  }
  if(is.Poisson(x)){
    yuima.warn("The cbind for Poisson will be implemented as soon as possible")
    return(NULL)
  }
}

# setMethod("initialize", "model.parameter",
#           function(.Object,
#                    all,
#                    common,
#                    diffusion,
#                    drift,
#                    jump,
#                    measure,
#                    xinit){
#             .Object@all <- all
#             .Object@common <- common
#             .Object@diffusion <- diffusion
#             .Object@drift <- drift
#             .Object@jump <- jump
#             .Object@measure <- measure
#             .Object@xinit <- xinit
#             return(.Object)
#           })

setMethod(
  "initialize", "model.parameter",
  function(.Object,
           all = character(),
           common = character(),
           diffusion = character(),
           drift = character(),
           jump = character(),
           measure = character(),
           xinit = character()) {
    .Object@all <- all
    .Object@common <- common
    .Object@diffusion <- diffusion
    .Object@drift <- drift
    .Object@jump <- jump
    .Object@measure <- measure
    .Object@xinit <- xinit
    return(.Object)
  }
)

# setMethod("initialize", "yuima.model",
#           function(.Object,
#                    drift ,
#                    diffusion,
# 				           hurst,
#                    jump.coeff,
#                    measure,
#                    measure.type,
#                    parameter,
#                    state.variable,
#                    jump.variable,
#                    time.variable,
#                    noise.number,
#                    equation.number,
#                    dimension,
#                    solve.variable,
#                    xinit,
#                    J.flag){
#             .Object@drift <- drift
#             .Object@diffusion <- diffusion
# 			.Object@hurst <- hurst
#             .Object@jump.coeff <- jump.coeff
#             .Object@measure <- measure
#             .Object@measure.type <- measure.type
#             .Object@parameter <- parameter
#             .Object@state.variable <- state.variable
#             .Object@jump.variable <- jump.variable
#             .Object@time.variable <- time.variable
#             .Object@noise.number <- noise.number
#             .Object@equation.number <- equation.number
#             .Object@dimension <- dimension
#             .Object@solve.variable <- solve.variable
#             .Object@xinit <- xinit
#             .Object@J.flag <- J.flag
#             return(.Object)
#           })

# 23/11 we need to provide the default values for the yuima.model object class
# in order to construct a new class that inherits from yuima.model class

setMethod(
  "initialize", "yuima.model",
  function(.Object,
           drift = expression(),
           diffusion = list(),
           hurst = 0.5,
           jump.coeff = list(),
           # jump.coeff = expression(),
           measure = list(),
           measure.type = character(),
           parameter = new("model.parameter"),
           state.variable = "x",
           jump.variable = "z",
           time.variable = "t",
           noise.number = numeric(),
           equation.number = numeric(),
           dimension = numeric(),
           solve.variable = character(),
           xinit = expression(),
           J.flag = logical()) {
    .Object@drift <- drift
    .Object@diffusion <- diffusion
    .Object@hurst <- hurst
    .Object@jump.coeff <- jump.coeff
    .Object@measure <- measure
    .Object@measure.type <- measure.type
    .Object@parameter <- parameter
    .Object@state.variable <- state.variable
    .Object@jump.variable <- jump.variable
    .Object@time.variable <- time.variable
    .Object@noise.number <- noise.number
    .Object@equation.number <- equation.number
    .Object@dimension <- dimension
    .Object@solve.variable <- solve.variable
    .Object@xinit <- xinit
    .Object@J.flag <- J.flag
    return(.Object)
  }
)


## setModel
## setter of class 'yuima.model'
## set yuima model from SDE
setModel <- function(drift = NULL,
                     diffusion = NULL,
                     hurst = 0.5,
                     jump.coeff = NULL,
                     measure = list(),
                     measure.type = character(),
                     state.variable = "x",
                     jump.variable = "z",
                     time.variable = "t",
                     solve.variable,
                     xinit = NULL,
                     model.class = NULL,
                     observed.variable = NULL,
                     unobserved.variable = NULL) {
  ## if model.class is NULL, decide model.class automatically
  if (is.null(model.class)) {
    if (!is.null(observed.variable) || !is.null(unobserved.variable)) {
      model.class <- "linearStateSpaceModel"

      model <- setStateSpaceModel(
        drift = drift,
        diffusion = diffusion,
        hurst = hurst,
        jump.coeff = jump.coeff,
        measure = measure,
        measure.type = measure.type,
        state.variable = state.variable,
        jump.variable = jump.variable,
        time.variable = time.variable,
        solve.variable = solve.variable,
        xinit = xinit,
        observed.variable = observed.variable,
        unobserved.variable = unobserved.variable
      )

      unobserved.variable <- model@state.variable[!model@is.observed]

      ## check if the model is linear
      drift.is.linear.with.unobserved <- tryCatch(
        is.linear(model@drift, unobserved.variable),
        warning = function(w) {
          return(TRUE)
        }
      )
      if (!drift.is.linear.with.unobserved) {
        model.class <- "stateSpaceModel"
      }
      # check that any of unobserved.variable is not in diffusion
      for (i in 1:length(model@diffusion)) {
        if (any(params_in_expr(params = unobserved.variable, expr = model@diffusion[[i]]))) {
          model.class <- "stateSpaceModel"
        }
      }
    } else {
      model.class <- "model"
    }
  }



  ## state space cases
  if (model.class == "stateSpaceModel") {
    model <- setStateSpaceModel(
      drift = drift,
      diffusion = diffusion,
      hurst = hurst,
      jump.coeff = jump.coeff,
      measure = measure,
      measure.type = measure.type,
      state.variable = state.variable,
      jump.variable = jump.variable,
      time.variable = time.variable,
      solve.variable = solve.variable,
      xinit = xinit,
      observed.variable = observed.variable,
      unobserved.variable = unobserved.variable
    )
    return(model)
  } else if (model.class == "linearStateSpaceModel") {
    model <- setLinearStateSpaceModel(
      drift = drift,
      diffusion = diffusion,
      hurst = hurst,
      jump.coeff = jump.coeff,
      measure = measure,
      measure.type = measure.type,
      state.variable = state.variable,
      jump.variable = jump.variable,
      time.variable = time.variable,
      solve.variable = solve.variable,
      xinit = xinit,
      observed.variable = observed.variable,
      unobserved.variable = unobserved.variable
    )
    return(model)
  } else if (model.class == "model") {
    # normal model case
    ## we need a temp env for simplifications
    mylengdumMeas <- length(measure)
    if (mylengdumMeas > 0) {
      for (i in c(1:mylengdumMeas)) {
        if (is(measure[[i]], "yuima.law")) {
          res <- aux.setModelLaw(drift, diffusion,
            hurst, jump.coeff, measure, measure.type,
            state.variable, jump.variable, time.variable,
            solve.variable, xinit,
            posyuimalaw = i
          )
          res@measure[[i]] <- measure[[i]]
          return(res)
        }
      }
    }

    yuimaENV <- new.env()
    ## ::measure and jump term #####################################

    ## ::initialize objects ########
    MEASURE <- list()

    ## ::end initialize objects ########

    ## ::make type function list ####
    CPlist <- c("dnorm", "dgamma", "dexp", "dconst")
    codelist <- c("rIG", "rNIG", "rgamma", "rbgamma", "rvgamma", "rstable", "rpts", "rnts")
    ## added "rpts" and "rnts" by YU (2016/10/4)
    ## ::end make type function list ####

    ## Multivariate YUIMA model


    if (!is.null(jump.coeff)) {
      if (is.matrix(jump.coeff)) {
        if (dim(jump.coeff)[2] != 1) {
          intensity <- NULL
          # if(is.null(names(measure)) || names(measure)=="df"){
          if (is.null(names(measure)) || all(names(measure) %in% "df")) {
            names(measure) <- "df"
          }

          df <- as.list(measure[["df"]])
          if (any(measure.type == "CP")) {
            intensity <- measure[["intensity"]]
          }
          my.cond <- TRUE
          tmp <- regexpr("\\(", measure$df)[1]
          measurefunc <- substring(measure$df, 1, tmp - 1)
          if (!is.na(match(measurefunc, codelist))) {
            my.cond <- FALSE
          }
          if (my.cond) {
            res <- setMultiModel(
              drift = drift, diffusion = diffusion,
              hurst = hurst, jump.coeff = jump.coeff,
              intensity = intensity, df = df,
              measure.type = measure.type, state.variable = state.variable,
              jump.variable = jump.variable, time.variable = time.variable,
              solve.variable = solve.variable, xinit = xinit
            )
            return(res)
          }
        }
      }
    }

    if (!length(measure.type)) {
      if (length(jump.coeff) || length(measure)) {
        yuima.warn("measure type does not match with jump term.")
        return(NULL)
      }
      jump.variable <- character()
      measure.par <- character()
    } else {
      if (!length(jump.coeff) || !length(measure)) {
        yuima.warn("measure type isn't matched with jump term.")
        return(NULL)
        # }else
        #       if(length(jump.coeff)!=1){
        #        yuima.warn("multi dimentional jump term is not supported yet.")
        #
        #         return(NULL)
        #     }
      } else if (measure.type == "CP") { ## ::CP
        #        if(length(measure)!=2){
        # yuima.warn(paste("length of measure must be two on type", measure.type, "."))
        # return(NULL)
        # }
        if (!is.list(measure)) {
          measure <- list(intensity = measure[1], df = measure[2], dimension = measure[3])
        } else {
          # if(length(measure[[1]])!=1 || length(measure[[2]])!=1){
          #   yuima.warn("multi dimentional jump term is not supported yet.")
          #   return(NULL)
          # }
          ## ::naming measure list ########
          tmpc <- names(measure)
          if (is.null(tmpc)) {
            names(measure) <- c("intensity", "df", "dimension")
          } else {
            whichint <- match("intensity", tmpc)
            whichdf <- match("df", tmpc)
            if (!is.na(whichint)) {
              if (names(measure)[-whichint] == "" || names(measure)[-whichint] == "df") {
                names(measure)[-whichint] <- "df"
              } else {
                yuima.warn("names of measure are incorrect.")
                return(NULL)
              }
            } else if (!is.na(whichdf)) {
              if (names(measure)[-whichdf] == "" || names(measure)[-whichdf] == "intensity") {
                names(measure)[-whichdf] <- "intensity"
              } else {
                yuima.warn("names of measure are incorrect.")
                return(NULL)
              }
            } else {
              yuima.warn("names of measure are incorrect.")
              return(NULL)
            }
          }
          ## ::end naming measure list ########
        }

        ## ::check df name ####################
        tmp <- regexpr("\\(", measure$df)[1]
        measurefunc <- substring(measure$df, 1, tmp - 1)
        if (!is.na(match(measurefunc, codelist))) {
          yuima.warn(paste("distribution function", measurefunc, "should be defined as type code."))
          return(NULL)
        } # else if(is.na(match(measurefunc, CPlist))){
        # warning(paste("\ndistribution function", measurefunc, "is not officialy supported as type CP.\n"))
        # }
        MEASURE$df$func <- eval(parse(text = measurefunc)) # LM 15/05/2017
        MEASURE$df$expr <- parse(text = measure$df)
        MEASURE$intensity <- parse(text = measure$intensity)

        measure.par <- unique(c(all.vars(MEASURE$intensity), all.vars(MEASURE$df$expr)))
        ## measure.par$intensity <- unique(all.vars(MEASURE$intensity))
        ## ::end check df name ####################
        ## ::end CP
      } else if (measure.type == "code") { ## ::code
        if (length(measure) != 1) {
          yuima.warn(paste("length of measure must be one on type", measure.type, "."))
          return(NULL)
        }

        if (!is.list(measure)) {
          measure <- list(df = measure)
        } else {
          if (length(measure[[1]]) != 1) {
            yuima.warn("multi dimentional jump term is not supported yet.")
            return(NULL)
          }
          ## ::naming measure list #############
          if (is.null(names(measure)) || names(measure) == "df") {
            names(measure) <- "df"
          } else {
            yuima.warn("name of measure is incorrect.")
            return(NULL)
          }
          ## ::end naming measure list #############
        }

        ## ::check df name ####################
        tmp <- regexpr("\\(", measure$df)[1]
        measurefunc <- substring(measure$df, 1, tmp - 1)
        if (!is.na(match(measurefunc, CPlist))) {
          yuima.warn(paste("\ndistribution function", measurefunc, "should be defined as type CP."))
          return(NULL)
        } else if (is.na(match(measurefunc, codelist))) {
          warning(paste("\ndistribution function", measurefunc, "is not officialy supported as type code.\n"))
        }
        ## MEASURE$df$func <- eval(parse(text=measurefunc))
        MEASURE$df$expr <- parse(text = measure$df)

        measure.par <- unique(all.vars(MEASURE$df$expr))
        ## ::end check df name ####################
        ## ::end code
      } else if (measure.type == "density") { ## ::density
        if (length(measure) != 1) {
          yuima.warn(paste("length of measure must be one on type", measure.type, "."))
          return(NULL)
        }

        if (!is.list(measure)) {
          measure <- list(df = measure)
        } else {
          if (length(measure[[1]]) != 1) {
            yuima.warn("multi dimentional jump term is not supported yet.")
            return(NULL)
          }

          ## ::naming measure list #############
          if (is.null(names(measure))) {
            names(measure) <- "df"
          } else if (names(measure) != "density" && names(measure) != "df") {
            yuima.warn("name of measure is incorrect.")
            return(NULL)
          }
          ## ::end naming measure list #############
        }

        ## ::check df name ####################
        tmp <- regexpr("\\(", measure[[names(measure)]])[1]
        measurefunc <- substring(measure[[names(measure)]], 1, tmp - 1)
        if (!is.na(match(measurefunc, CPlist))) {
          yuima.warn(paste("distribution function", measurefunc, "should be defined as type CP."))
          return(NULL)
        } else if (!is.na(match(measurefunc, codelist))) {
          yuima.warn(paste("distribution function", measurefunc, "should be defined as type code."))
          return(NULL)
        }
        MEASURE[[names(measure)]]$func <- eval(parse(text = measurefunc))
        MEASURE[[names(measure)]]$expr <- parse(text = measure[[names(measure)]])

        measure.par <- unique(all.vars(MEASURE[[names(measure)]]$expr))
        ## ::end check df name ####################
        ## ::end density
      } else { ## ::else
        yuima.warn(paste("measure type", measure.type, "isn't supported."))
        return(NULL)
      }
      n.eqn3 <- 1
      n.jump <- 1
    }




    ## ::end measure and jump term #####################################

    ## :: check for errors and reform values
    if (any(time.variable %in% state.variable)) {
      yuima.warn("time and state(s) variable must be different.")
      return(NULL)
    }
    if (is.null(dim(drift))) { # this is a vector
      n.eqn1 <- length(drift)
      n.drf <- 1
    } else { # it is a matrix
      n.eqn1 <- dim(drift)[1]
      n.drf <- dim(drift)[2]
    }

    if (is.null(dim(diffusion))) { # this is a vector
      n.eqn2 <- length(diffusion)
      n.noise <- 1
    } else { # it is a matrix
      n.eqn2 <- dim(diffusion)[1]
      n.noise <- dim(diffusion)[2]
    }

    if (is.null(diffusion)) {
      diffusion <- rep("0", n.eqn1)
      n.eqn2 <- n.eqn1
      n.noise <- 1
    }

    ## TBC
    n.eqn3 <- n.eqn1

    if (!length(measure)) {
      n.eqn3 <- n.eqn1
    }

    if (n.eqn1 != n.eqn2 || n.eqn1 != n.eqn3) {
      yuima.warn("Malformed model, number of equations in the drift and diffusion do not match.")
      return(NULL)
    }
    n.eqn <- n.eqn1

    if (is.null(xinit)) {
      # xinit <- numeric(n.eqn)
      xinit <- character(n.eqn)
    } else if (length(xinit) != n.eqn) {
      if (length(xinit) == 1) {
        xinit <- rep(xinit, n.eqn)
      } else {
        yuima.warn("Dimension of xinit variables missmatch.")
        return(NULL)
      }
    }

    if (missing(solve.variable)) {
      yuima.warn("Solution variable (lhs) not specified. Trying to use state variables.")
      solve.variable <- state.variable
    }
    if (n.eqn != length(solve.variable)) {
      yuima.warn("Malformed model, number of solution variables (lhs) do no match number of equations (rhs).")
      return(NULL)
    }

    loc.drift <- matrix(drift, n.eqn, n.drf)
    loc.diffusion <- matrix(diffusion, n.eqn, n.noise)
    # Modification starting point 6/11
    loc.xinit <- matrix(xinit, n.eqn, n.drf)

    ## :: allocate vectors
    DRIFT <- vector(n.eqn, mode = "expression")
    DIFFUSION <- vector(n.eqn, mode = "list")
    # Modification starting point 6/11
    XINIT <- vector(n.eqn, mode = "expression")

    ## :: function to make expression from drift characters
    pre.proc <- function(x) {
      for (i in 1:length(x)) {
        if (length(parse(text = x[i])) == 0) {
          x[i] <- "0"
        }
      }
      parse(text = paste(sprintf("(%s)", x), collapse = "+"))
    }
    ## 22/11:: function to simplify expression in drift, diffusion, jump and xinit characters
    yuima.Simplifyobj <- function(x) {
      dummy <- yuima.Simplify(x, yuima.env = yuimaENV)
      dummy1 <- yuima.Simplify(dummy, yuima.env = yuimaENV)
      dummy2 <- as.character(dummy1)
      res <- parse(text = paste0("(", dummy2, ")", collapse = NULL))
      return(res)
    }


    ## :: make expressions of drifts and diffusions and jump
    for (i in 1:n.eqn) {
      DRIFT[i] <- pre.proc(loc.drift[i, ])
      # 22/11 Simplify expressions
      DRIFT[i] <- yuima.Simplifyobj(DRIFT[i])
      # Modification starting point 6/11
      XINIT[i] <- pre.proc(loc.xinit[i, ])
      XINIT[i] <- yuima.Simplifyobj(XINIT[i])
      for (j in 1:n.noise) {
        expr <- parse(text = loc.diffusion[i, j])
        if (length(expr) == 0) {
          expr <- expression(0) # expr must have something
        }
        #       DIFFUSION[[i]][j] <- expr
        # 22/11
        DIFFUSION[[i]][j] <- yuima.Simplifyobj(expr)
      }
      # 22/11

      # if (length(JUMP)>0){
      #    JUMP[i] <- parse(text=jump.coeff[i])
      #    JUMP[i] <- yuima.Simplifyobj(JUMP[i])
      # }
    }



    # print(length(jump.coeff))
    # if (length(jump.coeff)==0){
    #    JUMP <- list(parse(text=jump.coeff))
    # }else{
    #    #    JUMP <- vector(n.eqn, mode="expression")
    #    JUMP <- vector(n.eqn, mode="list")
    # }

    if (length(jump.coeff) == 0) {
      JUMP <- list()
    } else {
      if (length(jump.coeff) == 1 & !is.matrix(jump.coeff)) { # is a scalar
        expr <- parse(text = jump.coeff)
        if (length(expr) == 0) {
          expr <- expression(0) # expr must have something
        }
        JUMP <- list(yuima.Simplifyobj(expr))
      } else { # must be matrix, n.col = dimension of Levy noise
        jump.coeff <- as.matrix(jump.coeff)
        c.j <- ncol(jump.coeff)
        r.j <- nrow(jump.coeff)
        # print(c.j)
        # print(r.j)
        # print(jump.coeff)
        JUMP <- vector(r.j, mode = "list")
        for (i in 1:r.j) {
          for (j in 1:c.j) {
            # cat(sprintf("\ni=%d,j=%d\n",i,j))
            expr <- parse(text = jump.coeff[i, j])
            if (length(expr) == 0) {
              expr <- expression(0) # expr must have something
            }
            JUMP[[i]][j] <- yuima.Simplifyobj(expr)
          }
        }
      }
    }
    # print(str(JUMP))

    #

    ## :: get parameters in drift expression
    drift.par <- unique(all.vars(DRIFT))
    # Modification starting point 6/11
    xinit.par <- unique(all.vars(XINIT))

    drift.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), drift.par)))
    if (length(drift.idx) > 0) {
      drift.par <- drift.par[-drift.idx]
    }

    ## :: get parameters in diffusion expression
    diff.par <- unique(unlist(lapply(DIFFUSION, all.vars)))
    diff.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), diff.par)))
    if (length(diff.idx) > 0) {
      diff.par <- diff.par[-diff.idx]
    }

    ## :: get parameters in jump expression
    J.flag <- FALSE
    #  jump.par <- unique(all.vars(JUMP))
    jump.par <- unlist(lapply(JUMP, all.vars))
    if (is.null(jump.par)) {
      jump.par <- character()
    }
    if (length(na.omit(match(jump.par, jump.variable)))) {
      J.flag <- TRUE
    }
    jump.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), jump.par)))
    if (length(jump.idx) > 0) {
      jump.par <- jump.par[-jump.idx]
    }

    ## :: get parameters in measure expression
    measure.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), measure.par)))
    if (length(measure.idx) > 0) {
      measure.par <- measure.par[-measure.idx]
    }

    ## :: order parameters for 'yuima.pars'
    ## id1 <- which(diff.par %in% drift.par)
    ## id2 <- which(drift.par %in% diff.par)
    ## common <- unique(c(diff.par[id1], drift.par[id2]))
    common <- c(drift.par, diff.par)
    common <- common[duplicated(common)]

    common1 <- common
    # modification 06/11 common1 contains only
    # parameters that appear in both drift and diffusion terms.

    # Modification 06/11 common contains only parameters that appear
    # in drift, diff, Jump and xinit
    if (length(xinit)) {
      common <- c(common, xinit.par)
      common <- common[duplicated(common)]
      common <- c(common, xinit.par)
      common <- common[duplicated(common)]
    }


    if (length(measure)) {
      common <- c(common, jump.par)
      common <- common[duplicated(common)]
      common <- c(common, measure.par)
      common <- common[duplicated(common)]
    }
    #   all.par <- unique(c(drift.par, diff.par, jump.par, measure.par))
    all.par <- unique(c(drift.par, diff.par, jump.par, measure.par, xinit.par))

    ## :: instanciate class
    tmppar <- new("model.parameter",
      all = all.par,
      #                 common= common,
      common = common1,
      diffusion = diff.par,
      drift = drift.par,
      jump = jump.par,
      measure = measure.par,
      xinit = xinit.par
    )
    tmp <- new("yuima.model",
      drift = DRIFT,
      diffusion = DIFFUSION,
      hurst = as.numeric(hurst),
      jump.coeff = JUMP,
      measure = MEASURE,
      measure.type = measure.type,
      parameter = tmppar,
      state.variable = state.variable,
      jump.variable = jump.variable,
      time.variable = time.variable,
      noise.number = n.noise,
      equation.number = n.eqn,
      dimension = c(
        length(tmppar@all),
        length(tmppar@common),
        length(tmppar@diffusion),
        length(tmppar@drift),
        length(tmppar@jump),
        length(tmppar@measure)
      ),
      solve.variable = solve.variable,
      xinit = XINIT,
      J.flag <- J.flag
    )
    return(tmp)
  }
}

aux.setModelLaw <- function(drift, diffusion,
                            hurst, jump.coeff, measure, measure.type,
                            state.variable, jump.variable, time.variable,
                            solve.variable, xinit, posyuimalaw) {
  dummyMeasure <- paste0(c(
    "yuima.law(",
    paste0(measure[[posyuimalaw]]@param.measure, collapse = ", "),
    ")"
  ), collapse = "")
  auxmeasure <- measure
  auxmeasure[[posyuimalaw]] <- dummyMeasure
  names(auxmeasure[posyuimalaw]) <- "df"
  setModel(
    drift = drift, diffusion = diffusion,
    hurst = hurst, jump.coeff = jump.coeff, measure = auxmeasure,
    measure.type = measure.type,
    state.variable = state.variable,
    jump.variable = jump.variable, time.variable,
    solve.variable, xinit
  )
}

# yuima.model rbind

# setGeneric("rbind.yuima",
#            function(x, ...)
#              standardGeneric("rbind.yuima")
# )

# setMethod("cbind.yuima", signature(x="yuima"),
#           function(x, ...){
#             ##:: init
#             y.list <- list(x, ...)
#             y.num <- length(y.list)
#
#             ##:: bind yuima.data in yuima
#
#             ##:: return result
#             return(NULL)
#           }
# )

# setMethod("rbind.yuima", signature(x="yuima.model"),
#           function(x, ...){
#             y.list <- list(x, ...)
#             y.num <- length(y.list)
#             res <- aux.rbind.model(y.list,y.num)
#             return(res)
#           }
# )

rbind.yuima.model <- function(x, ...) {
  y.list <- list(x, ...)
  #  y.list1 <- lapply(y.list, FUN = only.yuima.model)
  y.num <- length(y.list)
  new.list <- list()
  for (i in (1:y.num)) {
    if (is(y.list[[i]], "yuima.model")) {
      new.list[i] <- y.list[[i]]
    }
  }
  new.y.num <- length(new.list)
  res <- aux.rbind.model(
    y.list = new.list,
    y.num = new.y.num, mycall = y.list
  )
  return(res)
}

aux.rbind.model <- function(y.list, y.num, mycall = list()) {
  lapply(y.list, FUN = check.yuima.model)
  check.lev <- lapply(y.list, FUN = check.yuima.levy)
  check.lev <- unlist(check.lev)
  drift <- lapply(y.list, FUN = extract.model, type = "drift")
  diffusion <- lapply(y.list, FUN = extract.model, type = "diffusion")
  solve.variable <- lapply(y.list, FUN = extract.model, type = "solve.variable")
  state.variable <- lapply(y.list, FUN = extract.model, type = "state.variable")
  xinit <- lapply(y.list, FUN = extract.model, type = "xinit")
  noise.number <- lapply(y.list, FUN = extract.model, type = "noise.number")
  equation.number <- lapply(y.list, FUN = extract.model, type = "equation.number")
  # Until Here only diffusion process
  drift <- lapply(drift, FUN = ExpToString)
  drift <- unlist(drift)
  # drift
  nrow.diff <- sum(unlist(equation.number))
  ncol.diff <- sum(unlist(noise.number))
  matr.diff <- matrix("0", nrow = nrow.diff, ncol = ncol.diff)
  extrinf <- 1
  extrsup <- noise.number[[1]]
  j <- 1
  cond.eq <- equation.number[[1]]
  cond.eq1 <- 0
  for (i in c(1:nrow.diff)) {
    if (i <= cond.eq) {
      dum <- ExpToString(diffusion[[j]][[i - cond.eq1]])
      matr.diff[i, extrinf:extrsup] <- dum
      if (i == equation.number[[j]]) {
        extrinf <- extrsup + 1
        j <- j + 1
        if (j <= nrow.diff) {
          extrsup <- extrsup + equation.number[[j]]
          cond.eq1 <- i
          cond.eq <- cond.eq + equation.number[[j]]
        }
      }
    }
  }
  solve.variable <- lapply(solve.variable, FUN = ExpToString, cond = FALSE)
  solve.variable <- unlist(solve.variable)
  state.variable <- lapply(state.variable, FUN = ExpToString, cond = FALSE)
  state.variable <- unlist(state.variable)
  xinit <- lapply(xinit, FUN = ExpToString, cond = FALSE)
  xinit <- unlist(xinit)
  if (!any(check.lev)) {
    mod <- setModel(
      drift = drift, diffusion = matr.diff,
      solve.variable = solve.variable, state.variable = state.variable,
      xinit = xinit
    )
  } else {
    MultiLevy <- y.list[check.lev]
    jump.coeff <- lapply(MultiLevy,
      FUN = extract.model, type = "jump.coeff"
    )
    ncol.jump <- lapply(jump.coeff, FUN = numb.jump)
    dum.ncolj <- unlist(ncol.jump)
    ncol.jump <- sum(unlist(dum.ncolj))
    jump.coeff <- lapply(y.list,
      FUN = extract.model, type = "jump.coeff"
    )
    # ncol.jump1 <- lapply(jump.coeff, FUN = numb.jump)
    matr.jump <- matrix("0",
      nrow = nrow.diff,
      ncol = ncol.jump
    )
    j <- 1
    h <- 0
    cond.eqa <- equation.number[[j]]
    cond.eqb <- 0
    extrinf <- 1
    extrsup <- 1
    if (check.lev[j]) {
      extrsup <- dum.ncolj[j]
    } else {
      h <- h + 1
    }
    for (i in c(1:nrow.diff)) {
      if (i <= cond.eqa) {
        if (check.lev[j]) {
          dum <- ExpToString(jump.coeff[[j]][[i - cond.eqb]])
          matr.jump[i, extrinf:extrsup] <- dum
        } else {
          #          matr.jump[i,] <- matr.jump[i,]
        }
        if (i == cond.eqa) {
          cond.eqb <- i
          j <- j + 1
          if (j <= length(equation.number)) {
            cond.eqa <- cond.eqa + equation.number[[j]]
          }
          if (check.lev[j - 1]) {
            extrinf <- extrsup + 1
            extrsup <- extrsup + dum.ncolj[j - h]
          } else {
            extrinf <- extrinf
            extrsup <- extrsup
            h <- h + 1
          }
        }
      }
    }

    # mod <- matr.jump
    #     measure <- lapply(y.list,
    #        FUN = extract.model, type = "measure")
    #     measure
    df <- NULL
    if ("df" %in% names(mycall)) {
      df <- mycall$df
    }
    measure.type <- NULL
    if ("measure.type" %in% names(mycall)) {
      measure.type <- mycall$measure.type
    }
    intensity <- NULL
    if ("intensity" %in% names(mycall)) {
      intensity <- mycall$intensity
    }
    time.variable <- "t"
    if ("time.variable" %in% names(mycall)) {
      time.variable <- mycall$time.variable
    }
    mod <- setMultiModel(
      drift = drift, diffusion = matr.diff,
      jump.coeff = matr.jump, solve.variable = solve.variable,
      xinit = xinit, time.variable = time.variable, df = df,
      intensity = intensity, measure.type = measure.type
    )
  }
  return(mod)
}
# only.yuima.model<- function(y.list){
#   if(is(y.list,"yuima.model")){
#     return(y.list)
#   }else{
#     NULL
#   }
# }
numb.jump <- function(x) {
  length(x[[1]])
}

check.yuima.levy <- function(x) {
  Levy <- FALSE
  if (length(x@measure.type) > 0) {
    if (!is(x, "yuima.model")) {
      yuima.stop("the Levy model have to belong to the yuima.multimodel class")
    }
    Levy <- TRUE
  }
  return(Levy)
}

ExpToString <- function(x, cond = TRUE) {
  dum <- unlist(strsplit(toString(x), split = ", "))
  if (cond) {
    dum <- substr(dum, 2, nchar(dum) - 1)
  }
  return(dum)
}

extract.model <- function(x, type = "drift") {
  res <- slot(x, type)
  return(res)
}

check.yuima.model <- function(x) {
  if (is.CARMA(x)) {
    yuima.warn("The cbind for CARMA will be implemented as soon as possible")
    return(NULL)
  }
  if (is.COGARCH(x)) {
    yuima.warn("The cbind for COGARCH will be implemented as soon as possible")
    return(NULL)
  }
  if (is.Poisson(x)) {
    yuima.warn("The cbind for Poisson will be implemented as soon as possible")
    return(NULL)
  }
}
