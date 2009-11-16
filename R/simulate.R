##:: function for approximation of function G
gfunc <- function(x){
  c0 <- 10
  c1 <- 10
  ret <- rep(0, length(x))
  idx <- which(x < 1/c0)
  ret[idx] <- 1
  
  idx <- which(1/c0 <= x)
  ret[idx] <- 1-pnorm(x[idx])
  for(i in 1:length(idx)){
    n <- 1:floor(c1/x[idx[i]])
    ret[idx[i]] <- 4 * (ret[idx[i]] - sum( pnorm((4*n+1)*x[idx[i]]) - pnorm((4*n-1)*x[idx[i]]) ))
  }
  
  idx <- which(1 < ret)
  ret[idx] <- 1
  return(ret)
}


##:: simulate
##:: solves SDE and returns result
setGeneric("simulate",
           function(yuima, xinit, true.parameter, space.discretized=FALSE, increment.W=NULL, increment.L=NULL)
           standardGeneric("simulate")
           )
setMethod("simulate", "yuima", function(yuima, xinit, true.parameter, space.discretized=FALSE, increment.W=NULL, increment.L=NULL){
  ##:: error checks
  if(missing(yuima)){
    cat("\nyuima object is missing.\n")
    return(NULL)
  }
  
  sdeModel <- yuima@model
  Terminal <- yuima@sampling@Terminal[1]
  division <- yuima@sampling@division[1]
  
  ##   str(sdeModel)
  ##   print(Terminal)
  ##   print(yuima@sampling@division)
  ##   readline()
  
  if(missing(true.parameter)){
    true.parameter <- numeric(length(sdeModel@parameter@all))
  }
  
  r.size <- sdeModel@noise.number
  d.size <- sdeModel@equation.number
  
  ##:: error check
  if(missing(xinit)){
    xinit <- sdeModel@xinit
  }else if(length(xinit) != d.size){
    if(length(xinit)==1){
      xinit <- rep(xinit, d.size)
    }else{
      cat("\nDimension of xinit variables missmuch.\n")
      return(NULL)
    }
  }
  
  if(space.discretized){
    if(r.size>1){
      warning("Space-discretized EM cannot be used for multi-dimentional models. Use standard method.")
      space.discretized <- FALSE
    }
    if(length(sdeModel@jump.coeff)){
      warning("Space-discretized EM is for only Wiener Proc. Use standard method.")
      space.discretized <- FALSE
    }
  }

  ##:: Error check for increment specified version.
  if(!missing(increment.W)){
    if(space.discretized == TRUE){
      cat("\nParameter increment must be invalid if space.discretized=TRUE.\n")
      return(NULL)
    }else if(dim(increment.W)[1] != r.size){
      cat("\nLength of increment's row must be same as yuima@model@noise.number.\n")
      return(NULL)
    }else if(dim(increment.W)[2] != division){
      cat("\nLength of increment's column must be same as yuima@sampling@division[1].\n")
      return(NULL)
    }
  }

    ##:: Error check for increment specified version.
  if(!missing(increment.L)){
    if(space.discretized == TRUE){
      cat("\nParameter increment must be invalid if space.discretized=TRUE.\n")
      return(NULL)
    }else if(dim(increment.L)[1] != r.size){
      cat("\nLength of increment's row must be same as yuima@model@noise.number.\n")
      return(NULL)
    }else if(dim(increment.L)[2] != division){
      cat("\nLength of increment's column must be same as yuima@sampling@division[1].\n")
      return(NULL)
    }
  }
  
  
  ##:: check if DRIFT and/or DIFFUSION has values
  has.drift <- sum(as.character(sdeModel@drift) != "(0)")
  var.in.diff <- is.logical(any(match(unlist(lapply(sdeModel@diffusion, all.vars)), sdeModel@state.variable)))
  
  ##:: set variables
  modelstate <- sdeModel@solve.variable
  modeltime <- sdeModel@time.variable
  V0 <- sdeModel@drift
  V <- sdeModel@diffusion
  
  par.len <- length(sdeModel@parameter@all)
  if(par.len>0){
    for(i in 1:par.len){
      pars <- sdeModel@parameter@all[i]
      assign(pars, true.parameter[i])
    }
  }
  
  ##:: Initialization
  ##:: set time step
  delta <- Terminal/division
  ##:: initialize state variables
  dX <- xinit
  
  if(space.discretized){   ##:: using Space-discretized Euler-Maruyama method
    ## if(r.size > 1){
    ##   cat("\nSpace-discretized Euler-Maruyama method cannot be used for multi-dimentional models.\n")
    ##   return(NULL)
    ## }
    ##:: set numbers for making random time steps
    dxx <- 0.0001
    xx <- seq(0, 1.7, dxx)
    
    ##:: approximate function G(gg)
    gg <- gfunc(xx)
    appfunc <- suppressWarnings( approxfun(gg, xx) )
    
    ##:: calculate inverse of G
    unif.a <- runif(division*2)
    inv.a <- pmin(qnorm(1 - unif.a/4), appfunc(unif.a), na.rm=TRUE)
    
    ##:: make random time steps
    ep <- sqrt(delta)
    dTW <- (ep/inv.a)^2
    time_idx <- cumsum(dTW) ##:: time index should be attached            
    div_sd <- min(which(time_idx > Terminal)) ##:: cut by time=1
    time_idx <- time_idx[1:div_sd]
    
    ##:: add diffusion term
    dTW <- rbind(dTW[1:div_sd],
                 t(matrix( (rbinom(div_sd*r.size, 1, 0.5)*2-1) * ep,
                          nrow=div_sd,
                          ncol=r.size)
                   )
                 )
    
    X_mat <- matrix(0, d.size, div_sd+1)              
    X_mat[,1] <- dX
    
    ##:: function to calculate coefficients of dTW
    p.b <- function(t, X=numeric(d.size)){
      ##:: assign names of variables
      for(i in 1:length(modelstate)){
        assign(modelstate[i], X[i])
      }
      assign(modeltime,t)
      tmp <- matrix(0, d.size, r.size+1)
      for(i in 1:d.size){
        tmp[i,1] <- eval(V0[i])
        for(j in 1:r.size){
          tmp[i,j+1] <- eval(V[[i]][j])
        }
      }
      return(tmp)
    }
    ##:: calcurate difference equation
    for(i in 1:div_sd){
      dX <- dX + p.b(t=time_idx[i], X=dX) %*% dTW[,i]
      X_mat[,i+1] <- dX
    }
    ##tsX <- ts(data=t(X_mat), deltat=delta , start=0)
    ##:: output zoo data
    zooX <- zoo(x=t(X_mat), order.by=c(0, time_idx))
    yuimaData <- setData(original.data=zooX)
    yuima@data <- yuimaData
    return(yuima)
  }
  ##:: function to calculate coefficients of dW(including drift term)
  ##:: common used in Wiener and CP
  p.b <- function(t, X=numeric(d.size)){
    ##:: assign names of variables
    for(i in 1:length(modelstate)){
      assign(modelstate[i], X[i])
    }
    assign(modeltime, t)
    ##:: solve diffusion term
    if(has.drift){
      tmp <- matrix(0, d.size, r.size+1)
      for(i in 1:d.size){
        tmp[i,1] <- eval(V0[i])
        for(j in 1:r.size){
          tmp[i,j+1] <- eval(V[[i]][j])
        }
      }
    }else{  ##:: no drift term (faster)
      tmp <- matrix(0, d.size, r.size)
      for(i in 1:d.size){
        for(j in 1:r.size){
          tmp[i,j] <- eval(V[[i]][j])
        }
      }
    }
    return(tmp)
  }
  
  X_mat <- matrix(0, d.size, division+1)
  X_mat[,1] <- dX
  
  ##:: using Euler-Maruyama method
  
  ##:: Diffusion terms
  if( missing(increment.W)){
    dW <- rnorm(division * r.size, 0, sqrt(delta))
    dW <- t(matrix(dW, nrow=division, ncol=r.size))
  }else{
    dW <- increment.W
  }
  
  ## [TBC] Levy incrementが指定された場合にも, シミュレーションをincrementを用いて行うように.
  
  if(has.drift){  ##:: consider drift term to be one of the diffusion term(dW=1) 
    dW <- rbind( rep(1, division)*delta , dW)
  }

  if(!length(sdeModel@jump.coeff)){ ##:: Wiener Proc
    ##:: using Euler-Maruyama method
        
    if(var.in.diff){  ##:: diffusions have state variables
      ##:: calcurate difference eq.    
      for( i in 1:division){
        dX <- dX + p.b(t=i*delta, X=dX) %*% dW[, i]
        X_mat[,i+1] <- dX
      }
    }else{  ##:: diffusions have no state variables (not use p.b(). faster)
      sde.tics <- seq(0, Terminal, length=(division+1))
      sde.tics <- sde.tics[2:length(sde.tics)]
      
      X_mat[, 1] <- dX
      
      ##:: assign names of variables
      for(i in 1:length(modelstate)){
        assign(modelstate[i], dX[i])
      }
      assign(modeltime, sde.tics)
      t.size <- length(sde.tics)
      
      ##:: solve diffusion term
      if(has.drift){
        pbdata <- matrix(0, d.size*(r.size+1), t.size)
        for(i in 1:d.size){
          pbdata[(i-1)*(r.size+1)+1, ] <- eval(V0[i])
          for(j in 1:r.size){
            pbdata[(i-1)*(r.size+1)+j+1, ] <- eval(V[[i]][j])
          }
        }
        dim(pbdata)<-(c(r.size+1, d.size*t.size))
      }else{
        pbdata <- matrix(0, d.size*r.size, t.size)
        for(i in 1:d.size){
          for(j in 1:r.size){
            pbdata[(i-1)*r.size+j, ] <- eval(V[[i]][j])
          }
        }
        dim(pbdata)<-(c(r.size, d.size*t.size))
      }
    
      pbdata <- t(pbdata)
      
      ##:: calcurate difference eq.
      for( i in 1:division){
        dX <- dX + pbdata[((i-1)*d.size+1):(i*d.size), ] %*% dW[, i]
        X_mat[, i+1] <- dX
      }
    }
    tsX <- ts(data=t(X_mat), deltat=delta , start=0)
    
  }else{ ##:: Levy
    JP <- sdeModel@jump.coeff
    mu.size <- length(JP)
    
    ##:: function to solve c(x,z)
    p.b.j <- function(t, X=numeric(d.size)){
      for(i in 1:length(modelstate)){
        assign(modelstate[i], X[i])
      }
      assign(modeltime, t)
      tmp <- numeric(d.size)
      for(i in 1:d.size){
        tmp[i] <-  eval(JP[i])
      }
      return(tmp)
    }
	  
    if(sdeModel@measure.type == "CP"){ ##:: Compound-Poisson type
      eta0 <- eval(sdeModel@measure$intensity)
      ##:: get lambda from nu()
      lambda <- integrate(sdeModel@measure$df$func, 0, Inf)$value * eta0
      
      ##:: lambda = nu() (p6)
      N_sharp <- rpois(1,Terminal*eta0)	##:: Po(Ne)
      if(N_sharp == 0){
        JAMP <- FALSE
      }else{
        JAMP <- TRUE
        Uj <- sort( runif(N_sharp, 0, Terminal) )
        ij <- NULL
        for(i in 1:length(Uj)){
          Min <- min(which(c(1:division)*delta > Uj[i]))
          ij <- c(ij, Min)
        }
      }
      
      ##:: make expression to create iid rand J
      if(grep("^[dexp|dnorm|dgamma]", sdeModel@measure$df$expr)) {
        ##:: e.g. dnorm(z,1,1) -> rnorm(mu.size*N_sharp,1,1)
        F <- suppressWarnings(parse(text=gsub("^d(.+?)\\(.+?,", "r\\1(mu.size*N_sharp,", sdeModel@measure$df$expr, perl=TRUE)))
      }else{
        stop("Sorry. CP only supports dexp, dnorm and dgamma yet.")
      }
      randJ <- eval(F)
      j <- 1
      for(i in 1:division){
        if(JAMP==FALSE || sum(i==ij)==0){
          Pi <- 0
        }else{
          J <- eta0*randJ[j]/lambda
          j <- j+1
          ##cat(paste(J,"\n"))
          ##Pi <- zeta(dX,J)
          assign(sdeModel@jump.variable, J)
          ##Pi <- p.b.j(t=i*delta,X=dX) %*% J
          Pi <- p.b.j(t=i*delta, X=dX)
        }
        dX <- dX + p.b(t=i*delta, X=dX) %*% dW[, i] + Pi
        X_mat[, i+1] <- dX
      }
      tsX <- ts(data=t(X_mat), deltat=delta, start=0)
      
    }else if(sdeModel@measure.type=="code"){  ##:: code type
      ##:: Jump terms
      code <- suppressWarnings(sub("^(.+?)\\(.+", "\\1", sdeModel@measure$df$expr, perl=TRUE))
      args <- unlist(strsplit(suppressWarnings(sub("^.+?\\((.+)\\)", "\\1", sdeModel@measure$df$expr, perl=TRUE)), ","))
      dZ <- switch(code,
                   rNIG=paste("rNIG(division, ", args[2], ", ", args[3], ", ", args[4], "*delta, ", args[5], "*delta)"),
                   rIG=paste("rIG(division,", args[2], "*delta, ", args[3], ")"),
                   rgamma=paste("rgamma(division, ", args[2], "*delta, ", args[3], ")"),
                   rbgamma=paste("rbgamma(division, ", args[2], "*delta, ", args[3], ", ", args[4], "*delta, ", args[5], ")"),
##                   rngamma=paste("rngamma(division, ", args[2], "*delta, ", args[3], ", ", args[4], ", ", args[5], "*delta, ", args[6], ")"),
                   rngamma=paste("rngamma(division, ", args[2], "*delta, ", args[3], ", ", args[4], ", ", args[5], "*delta)"),
##                   rstable=paste("rstable(division, ", args[2], ", ", args[3], ", ", args[4], ", ", args[5], ", ", args[6], ")")
                   rstable=paste("rstable(division, ", args[2], ", ", args[3], ", ", args[4], "*delta^(1/",args[2],"), ", args[5], "*delta)")
                   )
      
      if(is.null(dZ)){  ##:: "otherwise"
        cat(paste("Code \"", code, "\" not supported yet.\n", sep=""))
        return(NULL)
      }
      dZ <- eval(parse(text=dZ))
      ##:: calcurate difference eq.
      
      for(i in 1:division){
        assign(sdeModel@jump.variable, dZ[i])
        dX <- dX + p.b(t=i*delta, X=dX) %*% dW[, i] +p.b.j(t=i*delta, X=dX) * dZ[i]
        X_mat[, i+1] <- dX
      }
      tsX <- ts(data=t(X_mat), deltat=delta, start=0)
	  
    }else{
      cat(paste("Type \"", sdeModel@measure.type, "\" not supported yet.\n", sep=""))
      return(NULL)
    }
  }
  
  yuimaData <- setData(original.data=tsX)
  yuima@data <- yuimaData
  return(yuima)
})
