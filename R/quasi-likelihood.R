##::quasi-likelihood function

##::extract drift term from yuima
##::para: parameter of drift term (theta2)
"calc.drift" <- function(yuima, para){
  r.size <- yuima@model@noise.number
  d.size <- yuima@model@equation.number
  modelstate <- yuima@model@state.variable
  modelpara <- yuima@model@parameter@drift
  DRIFT <- yuima@model@drift
  n <- length(yuima)[1]
  drift <- matrix(0,n,d.size)
  X <- as.matrix(onezoo(yuima))
  for(i in 1:length(modelpara)){
    assign(modelpara[i],para[i])
  }
  for(t in 1:n){
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


##::extract diffusion term from yuima
##::para: parameter of diffusion term (theta1)
"calc.diffusion" <- function(yuima, para){
  r.size <- yuima@model@noise.number
  d.size <- yuima@model@equation.number
  modelstate <- yuima@model@state.variable
  modelpara <- yuima@model@parameter@diffusion
  DIFFUSION <- yuima@model@diffusion
  n <- length(yuima)[1]
  diff <- array(0, dim=c(d.size, r.size, n))
  X <- as.matrix(onezoo(yuima))
  for(k in 1:length(modelpara)){
    assign(modelpara[k], para[k])
  }
  for(r in 1:r.size){
    for(t in 1:n){
      Xt <- X[t, ]
      for(d in 1:d.size){
        assign(modelstate[d], Xt[d])
      }
      for(d in 1:d.size){
        diff[d, r, t] <- eval(DIFFUSION[[d]][r])
      }
    }
  }
  return(diff)
}

##::calcurate diffusion%*%t(diffusion) matrix
calc.B <- function(diff){
  d.size <- dim(diff)[1]
  n <- dim(diff)[3]
  B <- array(0, dim=c(d.size, d.size, n))
  for(t in 1:n){
    B[, , t] <- diff[, , t]%*%t(diff[, , t])
  }
  return(B)
}

calc.B.grad <- function(yuima, para){
  r.size <- yuima@model@noise.number
  d.size <- yuima@model@equation.number
  modelstate <- yuima@model@state.variable
  modelpara <- yuima@model@parameter@diffusion
  DIFFUSION <- yuima@model@diffusion
  n <- length(yuima)[1]
  X <- as.matrix(onezoo(yuima))
  
  for(k in 1:length(modelpara)){
    assign(modelpara[k], para[k])
  }
  ##   B <- list(NULL)
  ##   for(d in 1:d.size){
  ##     B[[d]] <- list(NULL)
  ##     for(d2 in 1:d.size){
  ##       if(d2<d){
  ##         B[[d]][[d2]] <- B[[d2]][[d]]
  ##       }else{
  ##         B[[d]][[d2]] <- expression(0)
  ##         for(r in 1:r.size){
  ##           B[[d]][[d2]] <- parse(text=paste(as.character(B[[d]][[d2]]), "+(", as.character(DIFFUSION[[d]][r]), ")*(", as.character(DIFFUSION[[d2]][r]), ")"))
  ##         }
  ##       }
  ##     }
  ##   }
  
  B.grad <- array(0, dim=c(d.size, d.size, n, length(modelpara)))
  for(k in 1:length(modelpara)){
    for(t in 1:n){
      for(d in 1:d.size){
        assign(modelstate[d], X[t, d])
      }
      for(d in 1:d.size){
        for(d2 in 1:d.size){
          if(d2<d){
            B.grad[d, d2, t, k] <- B.grad[d2, d, t, k]
          }else{
            B <- expression(0)
            for(r in 1:r.size){
              B <- parse(text=paste(as.character(B), "+(", as.character(DIFFUSION[[d]][r]), ")*(", as.character(DIFFUSION[[d2]][r]), ")"))
            }
            B.grad[d, d2, t, k] <- eval(D(B, modelpara[k]))
          }
        }
      }
    }
  }
  return(B.grad)
}

##::calculate the log quasi-likelihood with parameters (theta2, theta1) and X.
##::yuima : yuima object
##::theta2 : parameter in drift.
##::theta1 : parameter in diffusion.
##::h : time width.
setGeneric("ql",
           function(yuima, theta2, theta1, h, print=FALSE, param)
           standardGeneric("ql")
           )
setMethod("ql", "yuima",
          function(yuima, theta2, theta1, h, print=FALSE, param){
            ##QLG <- ql.grad(yuima, theta2, theta1, h, print=FALSE)
            ##print(QLG)
            if( missing(yuima)){
              cat("\nyuima object is missing.\n")
              return(NULL)
            }
            
            ## param handling
            if( missing(param) ){
              if( missing(theta2) || missing(theta1) ){
                cat("\nparameters of yuima.model is missing.\n")
                return(NULL)
              }              
            }else{
              if( missing(theta2) && missing(theta1) ){
                if( !is.list(param) ){
                  cat("\nparam must be list.\n")
                  return(NULL)
                }

                if( length(param)!=2){
                  cat("\nlength of param is strange.\n")
                  return(NULL)
                }
                
                ## get theta2 and theta1 from param
                if( is.null(names(param)) ){
                  theta2 <- as.vector(param[[1]])
                  theta1 <- as.vector(param[[2]])
                }
                else if( sum( names(param)==c("theta2", "theta1") ) == 2 ){
                  theta2 <- as.vector(param[[1]])
                  theta1 <- as.vector(param[[2]])
                }
                else if( sum( names(param)==c("theta1", "theta2") ) == 2 ){
                  theta2 <- as.vector(param[[2]])
                  theta1 <- as.vector(param[[1]])
                }
                else{
                  cat("\nnames of param are strange.\n")
                  return(NULL)
                }
              }else{
                cat("\nConflict in parameter specification method.\n")
                return(NULL)
              }
            }
            ## END param handling
                
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
            n <- length(yuima)[1]
            X <- as.matrix(onezoo(yuima))
            deltaX <- matrix(0, n-1, d.size)
			for(t in 1:(n-1))
			 deltaX[t, ] <- X[t+1, ]-X[t, ]
			
#		  
#			  print(system.time(for(t in 1:(n-1))
#              deltaX[t, ] <- X[t+1, ]-X[t, ]
#            ))
#			print(system.time(deltaX1 <- diff(X)))
#				  print(deltaX)
#				  print(deltaX1)
				  
            if(is.nan(theta2) || is.nan(theta1)){
              stop("error: theta is not a namber in parameter matrix")
            }
            drift <- calc.drift(yuima, para=theta2)
            diff <- calc.diffusion(yuima, para=theta1)
            B <- calc.B(diff)
            
            QL <- 0
            pn <- numeric(n-1)
            for(t in 1:(n-1)){
              if(det(as.matrix(B[, , t]))==0){
                pn[t] <- log(1)
              }else{
                pn[t] <- log( 1/((2*pi*h)^(d.size/2)*det(as.matrix(B[, , t]))^(1/2)) *
                             exp((-1/(2*h))*t(deltaX[t, ]-h*drift[t, ])%*%solve(as.matrix(B[, , t]))%*%(deltaX[t,]-h*drift[t, ])) )
                QL <- QL+pn[t]
                if(pn[t]==-Inf && FALSE){
                  cat("t:", t, "\n")
                  cat("B[, , t]:", B[, , t], "\n")
                  cat("det(B):",det(as.matrix(B[, , t])), "\n")
                  cat("deltaX[t, ]", deltaX[t, ], "\n")
                  cat("drift[t, ]", drift[t, ], "\n")
                }
              }
            }
            if(QL==-Inf){
              warning("quasi likelihood is too small to calculate.")
            }
            if(print==TRUE){
              print(paste("QL:", QL, "  theta2:", theta2, "  theta1:", theta1))
            }
            return(QL)
          })


ql.grad <- function(yuima, theta2, theta1, h, print=FALSE){
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
  n <- length(yuima)[1]
  X <- as.matrix(onezoo(yuima))
  deltaX <- matrix(0, n-1, d.size)
  for(t in 1:(n-1)){
    deltaX[t, ] <- X[t+1, ]-X[t, ]
  }
  if(is.nan(theta2) || is.nan(theta1)){
    stop("error: theta is not a namber in parameter matrix")
  }
  drift <- calc.drift(yuima, para=theta2)
  diff <- calc.diffusion(yuima, para=theta1)
  B <- calc.B(diff)
  B.grad <- calc.B.grad(yuima, para=theta1)
  
  QLG <- numeric(dim(B.grad)[4])
  for(k in 1:length(QLG)){
    pg1 <- 0
    pg2 <- 0
    for(t in 1:(n-1)){
      B.tmp <- as.matrix(B[, , t])
      if(det(B.tmp)!=0){
        B.grad.tmp <- as.matrix(B.grad[, , t, k])
        aa <- as.matrix(deltaX[t, ]-h*drift[t, ])
        pg1 <- pg1 + sum(diag(solve(B.tmp)%*%B.grad.tmp))
        pg2 <- pg2 + t(aa)%*%solve(B.tmp)%*%B.grad.tmp%*%solve(B.tmp)%*%aa        
      }
    }
    QLG[k] <- (-1/2)*pg1 + 1/(2*h)*pg2
    
    if(QLG[k]==-Inf){
      warning(paste("gradient[", k, "] of quasi likelihood is too small too calculate.", sep=""))
    }
  }
  QLG <- QLG/sqrt(sum(QLG^2))
  if(print){
    print(paste("QLG:", QLG, "  theta2:", theta2, "  theta1:", theta1))
  }
  return(QLG)
}

##::calculate the relative log quasi-likelihood with parameters (theta2, theta1) and X.
##::yuima : yuima object
##::theta2 : parameter in drift.
##::theta1 : parameter in diffusion.
##::ptheta2 : parameter in drift in the prevous mcmc step.
##::ptheta1 : parameter in diffusion in the prevous mcmc step.
##::h : time width.
setGeneric("rql",
           function(yuima, theta2, theta1, ptheta2, ptheta1, h, print=FALSE, param, prevparam)
           standardGeneric("rql")
           )
setMethod("rql", "yuima",
          function(yuima, theta2, theta1, ptheta2, ptheta1, h, print=FALSE, param, prevparam){
            if(missing(yuima)){
              cat("\nyuima object is missing.\n")
              return(NULL)
            }
            
            ## param handling
            if( missing(param) ){
              if( missing(theta2) || missing(theta1) ){
                cat("\nparameters of yuima.model is missing.\n")
                return(NULL)
              }              
            }else{
              if( missing(theta2) && missing(theta1) ){
                if( !is.list(param) ){
                  cat("\nparam must be list.\n")
                  return(NULL)
                }

                if( length(param)!=2){
                  cat("\nlength of param is strange.\n")
                  return(NULL)
                }
                
                ## get theta2 and theta1 from param
                if( is.null(names(param)) ){
                  theta2 <- as.vector(param[[1]])
                  theta1 <- as.vector(param[[2]])
                }
                else if( sum( names(param)==c("theta2", "theta1") ) == 2 ){
                  theta2 <- as.vector(param[[1]])
                  theta1 <- as.vector(param[[2]])
                }
                else if( sum( names(param)==c("theta1", "theta2") ) == 2 ){
                  theta2 <- as.vector(param[[2]])
                  theta1 <- as.vector(param[[1]])
                }
                else{
                  cat("\nnames of param are strange.\n")
                  return(NULL)
                }
              }else{
                cat("\nConflict in parameter specification method.\n")
                return(NULL)
              }
            }
            ## END param handling
            
            ## prevparam handling
            if( missing(prevparam) ){
              if( missing(ptheta2) || missing(ptheta1) ){
                cat("\nparameters of yuima.model is missing.\n")
                return(NULL)
              }              
            }else{
              if( missing(ptheta2) && missing(ptheta1) ){
                if( !is.list(prevparam) ){
                  cat("\nparam must be list.\n")
                  return(NULL)
                }
                if( length(prevparam)!=2){
                  cat("\nlength of param is strange.\n")
                  return(NULL)
                }
                
                ## get theta2 and theta1 from param
                if( is.null(names(prevparam)) ){
                  ptheta2 <- as.vector(prevparam[[1]])
                  ptheta1 <- as.vector(prevparam[[2]])
                }
                else if( sum( names(prevparam)==c("ptheta2", "ptheta1") ) == 2 ){
                  ptheta2 <- as.vector(prevparam[[1]])
                  ptheta1 <- as.vector(prevparam[[2]])
                }
                else if( sum( names(prevparam)==c("ptheta1", "ptheta2") ) == 2 ){
                  ptheta2 <- as.vector(prevparam[[2]])
                  ptheta1 <- as.vector(prevparam[[1]])
                }
                else{
                  cat("\nnames of prevparam are strange.\n")
                  return(NULL)
                }
                
              }else{
                cat("\nConflict in parameter specification method.")
                return(NULL)
              }
            }
            ## END prevparam handling

            if(missing(h)){
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
            n <- length(yuima)[1]
            X <- as.matrix(onezoo(yuima))
            deltaX <- matrix(0, n-1, d.size)
            for(t in 1:(n-1)){
              deltaX[t, ] <- X[t+1, ]-X[t, ]
            }
            if(is.nan(theta2) || is.nan(theta1)){
              stop("error: theta is not a namber in parameter matrix")
            }
            drift <- calc.drift(yuima, para=theta2)
            diff <- calc.diffusion(yuima, para=theta1)
            B <- calc.B(diff)

            pdrift <- calc.drift(yuima, para=ptheta2)
            pdiff <- calc.diffusion(yuima, para=ptheta1)
            pB <- calc.B(pdiff)
            
            rQL <- 0
            pn <- numeric(n-1)
            for(t in 1:(n-1)){
              if(det(as.matrix(B[, , t]))*det(as.matrix(pB[, , t]))==0){
                pn[t] <- log(1)
              }else{
                pn[t] <- -log(det(as.matrix(B[, , t]))^(1/2))-1/(2*h)*t(deltaX[t, ]-h*drift[t, ])%*%solve(as.matrix(B[, , t]))%*%(deltaX[t, ]-h*drift[t, ])
                pn[t] <- pn[t]-(-log(det(as.matrix(pB[, , t]))^(1/2))-1/(2*h)*t(deltaX[t, ]-h*pdrift[t, ])%*%solve(as.matrix(pB[, , t]))%*%(deltaX[t, ]-h*pdrift[t, ]))
                rQL <- rQL+pn[t]
              }
            }
            if(print==TRUE){
              print(paste("relative QL:", rQL, "  theta2:", theta2, "  theta1:", theta1))
            }
            return(rQL)
          })

##::estimate parameters(theta2,theta1) with a constraint ui%*%theta-ci=0  
##::yuima : yuima object
##::theta2 : init parameter in drift term.
##::theta1 : init parameter in diffusion term.
##::h : length between each observation time.
##::theta1.lim, theta2.lim : limit of those parameters.
##::example: 0 <= theta1 <= 1 theta1.lim = matrix(c(0,1),1,2)
##::if theta1, theta2 are matrix, theta.lim can be defined matrix like rbind(c(0,1),c(0,1))
setGeneric("ml.ql",
           function(yuima, theta2, theta1, h, theta2.lim=matrix(c(0, 1), 1, 2), theta1.lim=matrix(c(0, 1), 1, 2), print=FALSE, BFGS=FALSE, param, interval)
           standardGeneric("ml.ql")
           )
setMethod("ml.ql", "yuima",
          function(yuima, theta2, theta1, h, theta2.lim=matrix(c(0, 1), 1, 2), theta1.lim=matrix(c(0, 1), 1, 2), print=FALSE, BFGS=FALSE, param, interval){
            if( missing(yuima)){
              cat("\nyuima object is missing.\n")
              return(NULL)
            }

            ## param handling
            if( missing(param) ){
              if( missing(theta2) || missing(theta1) ){
                cat("\nparameters of yuima.model is missing.\n")
                return(NULL)
              }              
            }else{
              if( missing(theta2) && missing(theta1) ){
                if( !is.list(param) ){
                  cat("\nparam must be list.\n")
                  return(NULL)
                }

                if( length(param)!=2){
                  cat("\nlength of param is strange.\n")
                  return(NULL)
                }
                
                ## get theta2 and theta1 from param
                if( is.null(names(param)) ){
                  theta2 <- as.vector(param[[1]])
                  theta1 <- as.vector(param[[2]])
                }
                else if( sum( names(param)==c("theta2", "theta1") ) == 2 ){
                  theta2 <- as.vector(param[[1]])
                  theta1 <- as.vector(param[[2]])
                }
                else if( sum( names(param)==c("theta1", "theta2") ) == 2 ){
                  theta2 <- as.vector(param[[2]])
                  theta1 <- as.vector(param[[1]])
                }
                else{
                  cat("\nnames of param are strange.\n")
                  return(NULL)
                }
              }else{
                cat("\nConflict in parameter specification method.\n")
                return(NULL)
              }
            }
            ## END param handling
            

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

            ## interval handling
            if( !missing(interval) ){
              if( missing(theta2.lim) && missing(theta2.lim) ){
                if( !is.list(interval) ){
                  cat("\ninterval must be list.\n")
                  return(NULL)
                }

                theta2.len <- length(yuima@model@parameter@drift)
                theta1.len <- length(yuima@model@parameter@diffusion)
                
                if( length(interval) !=  (theta2.len+theta1.len) ){
                  cat("\nlength of interval is strange.\n")
                  return(NULL)
                }
                
                ## get theta2.lim and theta1.lim from interval
                theta2.lim <- NULL
                theta1.lim <- NULL
                for(i in 1:theta2.len){
                  theta2.lim <- rbind(theta2.lim, interval[[i]])
                }
                for(i in 1:theta1.len){
                  theta1.lim <- rbind(theta1.lim, interval[[theta2.len+i]])
                }
              }else{
                cat("\nConflict in parameter specification method.\n")
                return(NULL)
              }
            }
            ## END interval handling

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
            
            ui <- rbind(diag(length(theta1)+length(theta2)), (-1)*diag(length(theta1)+length(theta2)))
            ci <- c(rbind(theta2.lim, theta1.lim))
            ci[(length(ci)/2+1):length(ci)] <- (-1)*ci[(length(ci)/2+1):length(ci)]
            
            ql.opt <- function(theta=c(theta2, theta1)){
              return(ql(yuima, theta2=theta[1:length(theta2)], theta1=theta[(length(theta2)+1):length(theta)], h=h, print=print))
            }
            ql.grad.opt <- function(theta=c(theta2, theta1)){
              return(ql.grad(yuima, theta2=theta[1:length(theta2)], theta1=theta[(length(theta2)+1):length(theta)], h=h, print=print))
            }
            
            if(BFGS){
              opt <- constrOptim(c(theta2, theta1), ql.opt, ql.grad.opt, ui=ui, ci=ci, control=list(fnscale=-1), method="BFGS", outer.iterations=500)
            }else{
              opt <- constrOptim(c(theta2, theta1), ql.opt, NULL, ui=ui, ci=ci, control=list(fnscale=-1), outer.iterations=500)
            }
            
            if(opt$convergence != 0){
              print("WARNING:optimization did not converge.")
            }
            return(opt)
          })

