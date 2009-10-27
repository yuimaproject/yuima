# quasi-likelihood function

# extract diffusion term from yuima
#para: parameter of diffusion term (theta1)
"calc.diffusion" <- function(yuima,para){
  r.size <- yuima@model@noise.number
  d.size <- yuima@model@equation.number
  modelstate <- yuima@model@state.variable
  modelpara <- yuima@model@parameter@diffusion
  DIFFUSION <- yuima@model@diffusion
  division <- length(yuima)[1]
  diff <- array(0,dim=c(d.size,r.size,division))
  X <- as.matrix(onezoo(yuima))
  for(i in 1:length(modelpara)){
    assign(modelpara[i],para[i])
  }
  for(r in 1:r.size){
    for(t in 1:division){
      Xt <- X[t,]
      for(d in 1:d.size){
        assign(modelstate[d],Xt[d])
      }
      for(d in 1:d.size){
        diff[d,r,t] <- eval(DIFFUSION[[d]][r])
      }
    }
  }
  return(diff)
}

#extract drift term from yuima
#para: parameter of drift term (theta2)
"calc.drift" <- function(yuima,para){
  r.size <- yuima@model@noise.number
  d.size <- yuima@model@equation.number
  modelstate <- yuima@model@state.variable
  modelpara <- yuima@model@parameter@drift
  DRIFT <- yuima@model@drift
  division <- length(yuima)[1]
  drift <- matrix(0,division,d.size)
  X <- as.matrix(onezoo(yuima))
  for(i in 1:length(modelpara)){
    assign(modelpara[i],para[i])
  }
  for(t in 1:division){
    Xt <- X[t,]
    for(d in 1:d.size){
      assign(modelstate[d],Xt[d])
    }
    for(d in 1:d.size){
      drift[t,d] <- eval(DRIFT[d])
    }
  }
  return(drift)  
}

#calcurate diffusion%*%t(diffusion) matrix
calc.B <- function(diff){
  d.size <- dim(diff)[1]
  division <- dim(diff)[3]
  B <- array(0,dim=c(d.size,d.size,division))
  for(t in 1:division){
    B[,,t] <- diff[,,t]%*%t(diff[,,t])
  }
  return(B)
}

#calculate the log quasi-likelihood with parameters (theta2,theta1) and X.
##yuima : yuima object
##theta2 : parameter in drift.
##theta1 : parameter in diffusion.
##h : time width.
setGeneric("ql",
           function(yuima,theta2,theta1,h,print=FALSE)
           standardGeneric("ql")
           )
setMethod("ql", "ANY",
          function(yuima,theta2,theta1,h,print=FALSE){
  if( missing(yuima)){
      cat("\nyuima object is missing.\n")
      return(NULL)
  }
  if( missing(theta2) || missing(theta1)){
      cat("\nparameters of yuima.model is missing.\n")
      return(NULL)
  }
  if( missing(h)){
      cat("\nlength of each time is missing.\n")
      return(NULL)
  }
  if(length(yuima@model@parameter@drift)!=length(theta2)){
      cat("\nlength of drift parameter is strange.\n")
      return(NULL)
  }
  if(length(yuima@model@parameter@diffusion)!=length(theta1)){
      cat("\nlength of diffusion parameter is strange.\n")
      return(NULL)
  }
  
  
  d.size <- yuima@model@equation.number
  division <- length(yuima)[1]
  X <- as.matrix(onezoo(yuima))
  deltaX <- matrix(0,division-1,d.size)
  for(t in 1:(division-1)){
    deltaX[t,] <- X[t+1,]-X[t,]
  }
  if(is.nan(theta2) || is.nan(theta1)){
    stop("error: theta is not a namber in parameter matrix")
  }
  drift <- calc.drift(yuima,para=theta2)
  diff <- calc.diffusion(yuima,para=theta1)
  B <- calc.B(diff)
  
  QL <- 0
  pn <- numeric(division-1)
  for(t in 1:(division-1)){
    if(det(as.matrix(B[,,t]))==0){
      pn[t] <- log(1)
    }else{
      pn[t] <- log(1/((2*pi*h)^(d.size/2)*det(as.matrix(B[,,t]))^(1/2)) *
                   exp((-1/(2*h))*t(deltaX[t,]-h*drift[t,])%*%solve(as.matrix(B[,,t]))%*%(deltaX[t,]-h*drift[t,])))
      QL <- QL+pn[t]
      if(pn[t]==-Inf && FALSE){
        cat("t:",t, "\n")
        cat("B[,,t]:",B[,,t], "\n")
        cat("det(B):",det(as.matrix(B[,,t])),"\n")
        cat("deltaX[t,]", deltaX[t,], "\n")
        cat("drift[t,]", drift[t,], "\n")
      }
    }
  }
  if(QL==-Inf){
    warning("quasi likelihood is too small too calculate.")
  }
  if(print==TRUE){
    print(paste("QL:",QL,"  theta2:",theta2,"  theta1:",theta1))
  }
  return(QL)
})


#calculate the relative log quasi-likelihood with parameters (theta2,theta1) and X.
##yuima : yuima object
##theta2 : parameter in drift.
##theta1 : parameter in diffusion.
##ptheta2 : parameter in drift in the prevous mcmc step.
##ptheta1 : parameter in diffusion in the prevous mcmc step.
##h : time width.
setGeneric("rql",
           function(yuima,theta2,theta1,ptheta2,ptheta1,h,print=FALSE)
           standardGeneric("rql")
           )
setMethod("rql", "ANY" ,
          function(yuima,theta2,theta1,ptheta2,ptheta1,h,print=FALSE){
  if( missing(yuima)){
      cat("\nyuima object is missing.\n")
      return(NULL)
  }
  if( missing(theta2) || missing(theta1)){
      cat("\nparameters of yuima.model is missing.\n")
      return(NULL)
  }
  if( missing(h)){
      cat("\nlength of each time is missing.\n")
      return(NULL)
  }
  if(length(yuima@model@parameter@drift)!=length(theta2)){
      cat("\nlength of drift parameter is strange.\n")
      return(NULL)
  }
  if(length(yuima@model@parameter@diffusion)!=length(theta1)){
      cat("\nlength of diffusion parameter is strange.\n")
      return(NULL)
  }
  
  
  d.size <- yuima@model@equation.number
  division <- length(yuima)[1]
  X <- as.matrix(onezoo(yuima))
  deltaX <- matrix(0,division-1,d.size)
  for(t in 1:(division-1)){
    deltaX[t,] <- X[t+1,]-X[t,]
  }
  if(is.nan(theta2) || is.nan(theta1)){
    stop("error: theta is not a namber in parameter matrix")
  }
  drift <- calc.drift(yuima,para=theta2)
  diff <- calc.diffusion(yuima,para=theta1)
  B <- calc.B(diff)

  pdrift <- calc.drift(yuima,para=ptheta2)
  pdiff <- calc.diffusion(yuima,para=ptheta1)
  pB <- calc.B(pdiff)
  
  rQL <- 0
  pn <- numeric(division-1)
  for(t in 1:(division-1)){
    if(det(as.matrix(B[,,t]))*det(as.matrix(pB[,,t]))==0){
      pn[t] <- log(1)
    }else{
      pn[t] <- -log(det(as.matrix(B[,,t]))^(1/2))-1/(2*h)*t(deltaX[t,]-h*drift[t,])%*%solve(as.matrix(B[,,t]))%*%(deltaX[t,]-h*drift[t,])
	  pn[t] <- pn[t]-(-log(det(as.matrix(pB[,,t]))^(1/2))-1/(2*h)*t(deltaX[t,]-h*pdrift[t,])%*%solve(as.matrix(pB[,,t]))%*%(deltaX[t,]-h*pdrift[t,]))
      rQL <- rQL+pn[t]
    }
  }
  if(print==TRUE){
    print(paste("relative QL:",rQL,"  theta2:",theta2,"  theta1:",theta1))
  }
  return(rQL)
})

#estimate parameters(theta2,theta1) with a constraint ui%*%theta-ci=0  
##yuima : yuima object
##theta2 : init parameter in drift term.
##theta1 : init parameter in diffusion term.
##h : length between each observation time.
##theta1.lim, theta2.lim : limit of those parameters.
###example: 0 <= theta1 <= 1 theta1.lim = matrix(c(0,1),1,2)
###if theta1, theta2 are matrix, theta.lim can be defined matrix like rbind(c(0,1),c(0,1))
setGeneric("ml.ql",
           function(yuima,theta2,theta1,h,theta2.lim=matrix(c(0,1),1,2),theta1.lim=matrix(c(0,1),1,2),print=FALSE)
           standardGeneric("ml.ql")
           )
setMethod("ml.ql", "ANY" ,
          function(yuima,theta2,theta1,h,theta2.lim=matrix(c(0,1),1,2),theta1.lim=matrix(c(0,1),1,2),print=FALSE){
  if( missing(yuima)){
      cat("\nyuima object is missing.\n")
      return(NULL)
  }
  if( missing(theta2) || missing(theta1)){
      cat("\nparameters of yuima.model is missing.\n")
      return(NULL)
  }
  if(length(yuima@model@parameter@drift)!=length(theta2)){
      cat("\nlength of drift parameter is strange.\n")
      return(NULL)
  }
  if(length(yuima@model@parameter@diffusion)!=length(theta1)){
      cat("\nlength of diffusion parameter is strange.\n")
      return(NULL)
  }
  if( missing(h)){
      cat("\nlength of each time is missing.\n")
      return(NULL)
  }
  if(is.matrix(theta2.lim) && is.matrix(theta1.lim)){
    if(ncol(theta2.lim)!=2 || ncol(theta1.lim)!=2){
      cat("\ntheta.lim is not available.\n")
      return(NULL)    
    }
  }

  if( length(theta2)!=1 && length(theta2)!=nrow(theta2.lim)){
      cat("\nsize of theta2 and theta2.lim are different.\n")
      return(NULL)    
  }
  if( length(theta1)!=1 && length(theta1)!=nrow(theta1.lim)){
      cat("\nsize of theta1 and theta1.lim are different.\n")
      return(NULL)    
  }
  
  

  ui <- rbind(diag(length(theta1)+length(theta2)),(-1)*diag(length(theta1)+length(theta2)))
  ci <- c(rbind(theta2.lim,theta1.lim))
  ci[(length(ci)/2+1):length(ci)] <- (-1)*ci[(length(ci)/2+1):length(ci)]
  
  ql.opt <- function(theta=c(theta2,theta1)){
    return(ql(yuima,theta2=theta[1:length(theta2)],theta1=theta[(length(theta2)+1):length(theta)],h=h,print=print))
  }
  
  opt <- constrOptim(c(theta2,theta1),ql.opt,NULL,ui=ui,ci=ci,control=list(fnscale=-1),outer.iterations=500)
  if( opt$convergence != 0){
  	print("WARNING:optimization did not converge.")
   }
  return(opt)
})

