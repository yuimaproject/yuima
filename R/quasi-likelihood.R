##::quasi-likelihood function

##::extract drift term from yuima
##::para: parameter of drift term (theta2)

### TO BE FIXED: all caculations should be made on a private environment to
### avoid problems.
### I have rewritten drift.term and diff.term instead of calc.drift and
### calc.diffusion to make them independent of the specification of the 
### parameters.  S.M.I. 22/06/2010

drift.term <- function(yuima, theta){
	r.size <- yuima@model@noise.number
	d.size <- yuima@model@equation.number
	modelstate <- yuima@model@state.variable
#	modelpara <- yuima@model@parameter@drift
	DRIFT <- yuima@model@drift
	n <- length(yuima)[1]
	drift <- matrix(0,n,d.size)
	X <- as.matrix(onezoo(yuima))
#	for(i in 1:length(modelpara)){
#		assign(modelpara[i],para[i])
#	}

	for(i in 1:length(theta)){
		assign(names(theta)[i],theta[[i]])
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



diffusion.term <- function(yuima, theta){
	r.size <- yuima@model@noise.number
	d.size <- yuima@model@equation.number
	modelstate <- yuima@model@state.variable
#	modelpara <- yuima@model@parameter@diffusion
	DIFFUSION <- yuima@model@diffusion
	n <- length(yuima)[1]
	diff <- array(0, dim=c(d.size, r.size, n))
	X <- as.matrix(onezoo(yuima))
#	for(k in 1:length(modelpara)){
#		assign(modelpara[k], para[k])
#	}
	for(i in 1:length(theta)){
		assign(names(theta)[i],theta[[i]])
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
              yuima.warn("yuima object is missing.")
              return(NULL)
            }
            
            ## param handling
            if( missing(param) ){
              if( missing(theta2) || missing(theta1) ){
                yuima.warn("Parameters of yuima.model are missing.")
                return(NULL)
              }              
            }else{
              if( missing(theta2) && missing(theta1) ){
                if( !is.list(param) ){
                  yuima.warn("param must be list.")
                  return(NULL)
                }

                if( length(param)!=2){
                  yuima.warn("length of param is strange.")
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
                  yuima.warn("names of param are strange.")
                  return(NULL)
                }
              }else{
                yuima.warn("Conflict in parameter specification method.")
                return(NULL)
              }
            }
            ## END param handling
                
            if( missing(h)){
              yuima.warn("length of each time is missing.")
              return(NULL)
            }

            if(length(yuima@model@parameter@drift)!=length(theta2)){
              yuima.warn("length of drift parameter is strange.")
              return(NULL)
            }
            
            if(length(yuima@model@parameter@diffusion)!=length(theta1)){
              yuima.warn("length of diffusion parameter is strange.")
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
    yuima.warn("yuima object is missing.")
    return(NULL)
  }
  if( missing(theta2) || missing(theta1)){
    yuima.warn("parameters of yuima.model are missing.")
    return(NULL)
  }
  if( missing(h)){
    yuima.warn("length of each time is missing.")
    return(NULL)
  }
  if(length(yuima@model@parameter@drift)!=length(theta2)){
    yuima.warn("length of drift parameter is strange.")
    return(NULL)
  }
  if(length(yuima@model@parameter@diffusion)!=length(theta1)){
    yuima.warn("length of diffusion parameter is strange.")
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
              yuima.warn("yuima object is missing.")
              return(NULL)
            }
            
            ## param handling
            if( missing(param) ){
              if( missing(theta2) || missing(theta1) ){
                yuima.warn("parameters of yuima.model is missing.")
                return(NULL)
              }              
            }else{
              if( missing(theta2) && missing(theta1) ){
                if( !is.list(param) ){
                  yuima.warn("param must be list.")
                  return(NULL)
                }

                if( length(param)!=2){
                  yuima.warn("length of param is strange.")
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
                  yuima.warn("names of param are strange.")
                  return(NULL)
                }
              }else{
                yuima.warn("Conflict in parameter specification method.")
                return(NULL)
              }
            }
            ## END param handling
            
            ## prevparam handling
            if( missing(prevparam) ){
              if( missing(ptheta2) || missing(ptheta1) ){
                yuima.warn("parameters of yuima.model is missing.")
                return(NULL)
              }              
            }else{
              if( missing(ptheta2) && missing(ptheta1) ){
                if( !is.list(prevparam) ){
                  yuima.warn("param must be list.")
                  return(NULL)
                }
                if( length(prevparam)!=2){
                  yuima.warn("length of param is strange.")
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
                  yuima.warn("names of prevparam are strange.")
                  return(NULL)
                }
                
              }else{
                yuima.warn("Conflict in parameter specification method.")
                return(NULL)
              }
            }
            ## END prevparam handling

            if(missing(h)){
              yuima.warn("length of each time is missing.")
              return(NULL)
            }
            if(length(yuima@model@parameter@drift)!=length(theta2)){
              yuima.warn("length of drift parameter is strange.")
              return(NULL)
            }
            if(length(yuima@model@parameter@diffusion)!=length(theta1)){
              yuima.warn("length of diffusion parameter is strange.")
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


## ml.ql by newton method.
newton.ml.ql <- function(yuima, theta2, theta1, h, iteration=1, param.only=FALSE, verbose=FALSE, ...){
    
  ## get param
  r.size <- yuima@model@noise.number
  d.size <- yuima@model@equation.number
  modelstate <- yuima@model@state.variable
  modelpara.drift <- yuima@model@parameter@drift
  modelpara.diff <- yuima@model@parameter@diffusion

  ## get expression of functions ##
  ## a
  a <- yuima@model@drift
  
  ## b
  b <- yuima@model@diffusion

  ## B = b %*% t(b)
  if(length(b)>=2){
    ## expression matrix calculation
    B <- list()
    for( i in 1:d.size){
      for( j in i:d.size){
        tmp <- NULL   
        for(l in 1:r.size){
          B.il <- as.character(b[[i]][l])
          B.jl <- as.character(b[[j]][l])
          if(l==1){
            tmp <- paste(B.il, "*", B.jl)
          }else{
            tmp <- paste(tmp, "+", B.il, "*", B.jl)
          }
        }
        ## update B list
        B[[ (i-1)*d.size + j ]] <- parse(text=tmp)        
        if(i!=j) B[[ (j-1)*d.size + i ]] <- parse(text=tmp)
      }
    }
    dim(B) <- c(d.size, d.size)
  }else{
    b <- yuima@model@diffusion[[1]]
    B <- parse(text=paste(as.character(b),
                 " * ", as.character(b), sep=""))
    B <- list(B)
    dim(B) <- c(1,1)
  }
  
  ## some func def about B
  deriv.B <- function(B, var=character()){
    B.deriv <- B
    for(i in 1:nrow(B)){
      for( j in 1:ncol(B) ){
        B.deriv[i,j][[1]] <- D(B[i,j][[1]], var)
      }
    }
    return(B.deriv)
  }
  eval.B <- function(B, theta1, theta2, Xt.iMinus1, ...){
    ##assign variable
    for(j in 1:length(Xt.iMinus1)){
      assign(modelstate[j], Xt.iMinus1[j])
    }
    for(j in 1:length(theta1)){
      assign(modelpara.diff[j], theta1[j])
    }
    for(j in 1:length(theta2)){
      assign(modelpara.drift[j], theta2[j])
    }
    B.subs <- matrix(0, nrow(B), ncol(B))
    for( i in 1:nrow(B) ){
      for( j in 1:ncol(B) ){
        B.subs[i,j] <- eval(B[i,j][[1]])
      }
    }
    return(B.subs)
  }
  eval.inv.B <- function(B, theta1, theta2, Xt.iMinus1, ...)
    return(solve(eval.B(B, theta1, theta2, Xt.iMinus1, ...)))
  ## END
  
  ## some func about a
  deriv.a <- function(a, var=character()){
    a.deriv <- NULL
    for(i in 1:length(a)){
      a.deriv <- c(a.deriv, as.expression(D(a[i], var)))
    }
    return(a.deriv)
  }
  
  eval.a <- function(a, theta1, theta2, Xt.iMinus1, ...){
    ##assign variable
    for(j in 1:length(Xt.iMinus1)){
      assign(modelstate[j], Xt.iMinus1[j])
    }
    for(j in 1:length(theta1)){
      assign(modelpara.diff[j], theta1[j])
    }
    for(j in 1:length(theta2)){
      assign(modelpara.drift[j], theta2[j])
    }
    a.subs <- matrix(0, length(a), 1)
    for(i in 1:length(a)){
      a.subs[i,] <- eval(a[i])      
    }
    return(a.subs)
  }
  ##END
  
  ## dGi_theta1
  dGi.theta1 <- function(theta1, theta2, h, Xt.iMinus1, delta.i.X, B, a, ...){    
    ##d_theta1 log(detB)+d_theta1 1/2h*(-)T%*%B^(-1)%*%(-)
    ## 2nd term d_theta1 log(det(B))
    term.2 <- matrix(0, length(theta1), 1)
    for( j in 1:length(theta1)){
      term.2[j,1] <- sum(diag(eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
                              eval.B(deriv.B(B, modelpara.diff[j]), theta1, theta2, Xt.iMinus1)))
    }
    
    ## 3rd term d_theta1 1/2h *(-)B^(-1)(-)
    term.3 <- matrix(0, length(theta1),1)
    for( j in 1:length(theta1)){
      tmp <- delta.i.X - h * eval.a(a, theta1, theta2, Xt.iMinus1)
      term.3[j,1] <-  -1/(2*h) * t(tmp) %*% eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
          eval.B(deriv.B(B, modelpara.diff[j]), theta1, theta2, Xt.iMinus1) %*%
            eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*% tmp
    }

    ret <- 1/2*term.2+term.3
    return(ret)
  }
  
  ## dH_theta1
  dH.theta1 <- function(theta1, theta2, h, yuima, B, a, ...){
    ret <- matrix(0, length(theta1), 1)
    data <- as.matrix(onezoo(yuima))
    for( i in 1:(nrow(data)-1)){
      ##Xt.iMinus1
      Xt.iMinus1 <- data[i,]
      ##delta.i.X
      delta.i.X <- data[i+1,] - data[i,]
      ##calc
      ret <- ret + dGi.theta1(theta1, theta2, h, Xt.iMinus1, delta.i.X, B, a, ...)
    }
    ret <- -ret
    return(ret)
  }

  
  ## d2Gi_theta1
  d2Gi.theta1 <- function(theta1, theta2, h, Xt.iMinus1, delta.i.X, B, a,...){
    ##2nd term
    term.2 <- matrix(0, length(theta1), length(theta1))
    for(j in 1:length(theta1)){
      for(k in 1:length(theta1)){
        tmp <- -eval.inv.B(B,theta1, theta2, Xt.iMinus1) %*%
          eval.B(deriv.B(B,modelpara.diff[k]), theta1, theta2, Xt.iMinus1) %*%
          eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
            eval.B(deriv.B(B,modelpara.diff[j]), theta1, theta2, Xt.iMinus1) +
            eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
              eval.B( deriv.B(deriv.B(B,modelpara.diff[j]), modelpara.diff[k]),
                     theta1, theta2, Xt.iMinus1)
        term.2[j,k] <- sum(diag(tmp)) / 2
      }
    }
    ##3rd term
    term.3 <- matrix(0, length(theta1), length(theta1))
    for(j in 1:length(theta1)){
      for(k in 1:length(theta1)){
        tmp <- -2 * eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
          eval.B( deriv.B(B, modelpara.diff[k]), theta1, theta2, Xt.iMinus1 ) %*%
            eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
              eval.B( deriv.B(B, modelpara.diff[j]), theta1, theta2, Xt.iMinus1) %*%
                eval.inv.B(B, theta1, theta2, Xt.iMinus1) +
                  eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
                    eval.B( deriv.B(deriv.B(B,modelpara.diff[j]), modelpara.diff[k]),
                           theta1, theta2, Xt.iMinus1 ) %*%
                             eval.inv.B(B, theta1, theta2, Xt.iMinus1)
        tmp2 <- delta.i.X -h * eval.a(a, theta1, theta2, Xt.iMinus1)
        term.3[j,k] <- - 1 / (2*h) * t(tmp2) %*% tmp %*% tmp2 
      }
    }

    ##ret
    ret <- term.2+term.3
    return(ret)
  }
  
  ## d2H_theta1
  d2H.theta1 <-function(theta1, theta2, h, yuima, B, a, ...){
    ret <- matrix(0, length(theta1), length(theta1))
    data <- as.matrix(onezoo(yuima))
    for(i in 1:(nrow(data)-1)){
      ##Xt.iMinus1
      Xt.iMinus1 <- data[i,]
      ##delta.i.X
      delta.i.X <- data[i+1,] - data[i,]
      ##calc
      ret <- ret + d2Gi.theta1(theta1,  theta2, h, Xt.iMinus1, delta.i.X, B, a,...)
    }
    ret <- -ret
    return(ret)
  }

  ## dGi_theta2
  dGi.theta2 <- function(theta1, theta2, h, Xt.iMinus1, delta.i.X, B, a, ...){
    ##calc
    ret <- matrix(0, length(theta2), 1)
    for( j in 1:length(theta2) ){
      ret[j,1] <- -t(delta.i.X - h * eval.a(a,theta1, theta2, Xt.iMinus1)) %*%
        eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
          eval.a( deriv.a(a, modelpara.drift[j]), theta1, theta2, Xt.iMinus1)
    }

    return(ret)
  }
  
  ## dH_theta2
  dH.theta2 <-function(theta1, theta2, h, yuima, B, a, ...){
    ret <- matrix(0, length(theta2), 1)
    data <- as.matrix(onezoo(yuima))
    for( i in 1:(nrow(data)-1)){
      ##Xt.iMinus1
      Xt.iMinus1 <- data[i,]
      ##delta.i.X
      delta.i.X <- data[i+1,] - data[i,]
      ##calc
      ret <- ret + dGi.theta2(theta1, theta2, h, Xt.iMinus1, delta.i.X, B, a, ...)
      #cat("dH.theta2-tmp:", ret, "\n")
    }
    ret <- -ret
    return(ret)
  }
  
  ## d2Gi_theta2
  d2Gi.theta2 <- function(theta1, theta2, h, Xt.iMinus1, delta.i.X, B, a, ...){
    ##calc
    ret <- matrix(0, length(theta2), length(theta2))
    for(j in 1:length(theta2)){
      for(k in 1:length(theta2)){
        ret[j,k] <- - t(delta.i.X) %*% eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
          eval.a( deriv.a( deriv.a(a,modelpara.drift[j]), modelpara.drift[k]),
                 theta1, theta2, Xt.iMinus1) +
                   h * 
                     t(eval.a(deriv.a(a, modelpara.drift[j]), theta1, theta2, Xt.iMinus1)) %*%
                       eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
                         eval.a(deriv.a(a, modelpara.drift[k]), theta1, theta2, Xt.iMinus1) +
                           h *
                             t(eval.a(a, theta1, theta2, Xt.iMinus1)) %*%
                               eval.inv.B(B, theta1, theta2, Xt.iMinus1) %*%
                                 eval.a( deriv.a(deriv.a(a,modelpara.drift[j]), modelpara.drift[k]),
                                        theta1, theta2, Xt.iMinus1)
      }
    }
    return(ret)
  }
  
  ## d2H_theta2
  d2H.theta2 <- function(theta1, theta2, h, yuima, B, a, ...){
    ret <- matrix(0, length(theta2), length(theta2))
    data <- as.matrix(onezoo(yuima))
    for(i in 1:(nrow(data)-1)){
      ##Xt.iMinus1
      Xt.iMinus1 <- data[i,]
      ##delta.i.X
      delta.i.X <- data[i+1,] - data[i,]
      ##calc
      ret <- ret + d2Gi.theta2(theta1,  theta2, h, Xt.iMinus1, delta.i.X, B, a, ...)
    }
    ret <- -ret
    return(ret)
  }
  ## END function define ##

  ## newton algorithm main ##
  if(!param.only){
    if(verbose){
      #cat("theta2 init:", theta2, "\n")
      cat("theta1 init:", theta1, "\n")
      cat("theta2 init:", theta2, "\n")
    }
    
    for(ite in 1:iteration){
      
      #dHtheta2 <- dH.theta2(theta1, theta2, h, yuima, B, a, ...)
      #d2Htheta2 <- d2H.theta2(theta1, theta2, h, yuima, B, a, ...)
      #theta2 <- as.vector(theta2 - solve(d2Htheta2) %*% dHtheta2)
      
      dHtheta1 <- dH.theta1(theta1, theta2, h, yuima, B, a, ...)
      d2Htheta1 <- d2H.theta1(theta1, theta2, h, yuima, B, a, ...)
      theta1 <- as.vector(theta1 - solve(d2Htheta1) %*% dHtheta1)

      dHtheta2 <- dH.theta2(theta1, theta2, h, yuima, B, a, ...)
      d2Htheta2 <- d2H.theta2(theta1, theta2, h, yuima, B, a, ...)
      theta2 <- as.vector(theta2 - solve(d2Htheta2) %*% dHtheta2)
      
      if(verbose){
        cat("\n## Iteration", ite, "##\n")
        #cat("theta2 new:", theta2, "\n")
        cat("theta1 new:", theta1, "\n")
        cat("theta2 new:", theta2, "\n")
      }
    }
  }
  ## END newtom algorithm ##

  ##calc Sigma1, Sigma2, Hessian?##
  Sigma1 <- - solve(d2H.theta1(theta1, theta2, h, yuima, B, a, ...))
  Sigma2 <- - solve(d2H.theta2(theta1, theta2, h, yuima, B, a, ...))
  hessian <- 1 ## Not Implemented yet!
  ##END 
  
  ##temp
  return(list(theta1.new=theta1, theta2.new=theta2, Sigma1=Sigma1, Sigma2=Sigma2, hessian=hessian))

}
## END ml.ql by newton method

##::estimate parameters(theta2,theta1) with a constraint ui%*%theta-ci=0  
##::yuima : yuima object
##::theta2 : init parameter in drift term.
##::theta1 : init parameter in diffusion term.
##::h : length between each observation time.
##::theta1.lim, theta2.lim : limit of those parameters.
##::example: 0 <= theta1 <= 1 theta1.lim = matrix(c(0,1),1,2)
##::if theta1, theta2 are matrix, theta.lim can be defined matrix like rbind(c(0,1),c(0,1))
setGeneric("ml.ql",
           function(yuima, theta2, theta1, h, theta2.lim=matrix(c(0, 1), 1, 2),
                    theta1.lim=matrix(c(0, 1), 1, 2), print=FALSE,
                    method="optim",
                    param, interval)
           standardGeneric("ml.ql")
           )
setMethod("ml.ql", "yuima",
          function(yuima, theta2, theta1, h, theta2.lim=matrix(c(0, 1), 1, 2),
                   theta1.lim=matrix(c(0, 1), 1, 2), print=FALSE,
                   method="optim", #(BFGS, Newton)
                   param, interval){
            if( missing(yuima)){
              yuima.warn("yuima object is missing.")
              return(NULL)
            }

            ## param handling
            if( missing(param) ){
              if( missing(theta2) || missing(theta1) ){
                yuima.warn("parameters of yuima.model is missing.")
                return(NULL)
              }              
            }else{
              if( missing(theta2) && missing(theta1) ){
                if( !is.list(param) ){
                  yuima.warn("param must be list.")
                  return(NULL)
                }

                if( length(param)!=2){
                  yuima.warn("length of param is strange.")
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
                  yuima.warn("names of param are strange.")
                  return(NULL)
                }
              }else{
                yuima.warn("Conflict in parameter specification method.")
                return(NULL)
              }
            }
            ## END param handling
            

            if(length(yuima@model@parameter@drift)!=length(theta2)){
              yuima.warn("length of drift parameter is strange.")
              return(NULL)
            }
            if(length(yuima@model@parameter@diffusion)!=length(theta1)){
              yuima.warn("length of diffusion parameter is strange.")
              return(NULL)
            }
            if( missing(h)){
              yuima.warn("length of each time is missing.")
              return(NULL)
            }

            ## interval handling
            if( !missing(interval) ){
              if( missing(theta2.lim) && missing(theta2.lim) ){
                if( !is.list(interval) ){
                  yuima.warn("interval must be list.")
                  return(NULL)
                }

                theta2.len <- length(yuima@model@parameter@drift)
                theta1.len <- length(yuima@model@parameter@diffusion)
                
                if( length(interval) !=  (theta2.len+theta1.len) ){
                  yuima.warn("length of interval is strange.")
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
                yuima.warn("Conflict in parameter specification method.")
                return(NULL)
              }
            }
            ## END interval handling

            ql.opt <- function(theta=c(theta2, theta1)){
              return(ql(yuima, theta2=theta[1:length(theta2)],
                        theta1=theta[(length(theta2)+1):length(theta)], h=h, print=print))
            }
            ql.grad.opt <- function(theta=c(theta2, theta1)){
              return(ql.grad(yuima, theta2=theta[1:length(theta2)],
                             theta1=theta[(length(theta2)+1):length(theta)], h=h, print=print))
            }
            
            if(method=="Newton"){
              opt <- newton.ml.ql(yuima, theta2, theta1, h,
                                  iteration=10, param.only=FALSE, verbose=print)

              coef <- c(opt$theta2.new, opt$theta1.new)
              min <- ql(yuima, opt$theta2.new, opt$theta1.new, h, print, param)
              
            }else{ ## optim
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

              ql.opt.theta1 <- function(theta1){
                return(ql(yuima, theta2=theta2,
                          theta1=theta1, h=h, print=print))
              }

              ql.opt.theta2 <- function(theta2){
                return(ql(yuima, theta2=theta2,
                          theta1=theta1, h=h, print=print))
              }
              
              ##theta1 estim
              if(length(theta1) != 1){
                ui <- rbind(diag(length(theta1)), (-1)*diag(length(theta1)))
                ci <- c(theta1.lim)
                ci[(length(ci)/2+1):length(ci)] <- (-1)*ci[(length(ci)/2+1):length(ci)]
                opt1 <- constrOptim(c(theta1), ql.opt.theta1, NULL, ui=ui,
                                    ci=ci, control=list(fnscale=-1), outer.iterations=500)
                if(opt1$convergence != 0) print("WARNING:optimization did not converge.")
                theta1 <- opt1$par
              }else{
                cat("\none-diml optimization : Initial value (theta1) is ignored.\n")
                opt1 <- optimize(ql.opt.theta1, interval=theta1.lim, maximum=TRUE)
                theta1 <- opt1$maximum
              }

              #theta2 estim
              if(length(theta2) != 1){
                ui <- rbind(diag(length(theta2)), (-1)*diag(length(theta2)))
                ci <- c(theta2.lim)
                ci[(length(ci)/2+1):length(ci)] <- (-1)*ci[(length(ci)/2+1):length(ci)]
                opt2 <- constrOptim(c(theta2), ql.opt.theta2, NULL, ui=ui,
                                    ci=ci, control=list(fnscale=-1), outer.iterations=500)
                if(opt2$convergence != 0) print("WARNING:optimization did not converge.")

                theta2 <- opt2$par
                min <- opt2$value
                
              }else{
                cat("\none-diml optimization : Initial value (theta2) is ignored.\n")
                opt2 <- optimize(ql.opt.theta2, interval=theta2.lim, maximum=TRUE)
                theta2 <- opt2$maximum
                min <- opt2$objective
              }

              opt <- NULL
              opt$opt1 <- opt1
              opt$opt2 <- opt2
              coef <- c(theta2, theta1)

              #theta2 <- coef[1:length(yuima@model@parameter@drift)]
              #theta1 <- coef[(length(yuima@model@parameter@drift)+1):length(coef)]

              opt3 <- newton.ml.ql(yuima, theta2, theta1, h,
                                   iteration=0, param.only=TRUE, verbose=print)
              opt$Sigma1 <- opt3$Sigma1
              opt$Sigma2 <- opt3$Sigma2
              opt$hessian <- opt3$hessian
            }

            ##convert to mle object
            call <- match.call()
            names(coef) <- yuima@model@parameter@all
            fullcoef <- coef
            method <- method

            hessian <- opt$hessian
            vcov <- if(length(coef)) solve(hessian)
            else matrix(numeric(0L), 0L, 0L)
            
            opt <- new("mle", call=call, coef=coef, fullcoef=fullcoef,
                       vcov=vcov, min=min, details=opt, minuslogl=ql.opt,
                       method=method)
            ##END convert to mle
            
            return(opt)
          })

setMethod("confint", "yuima",
          function(object, parm, level=0.95, ...){
            yuima <- object
            #print("This is yuima confint !")
            Sigma1 <- yuima@data@mle@details$Sigma1
            Sigma2 <- yuima@data@mle@details$Sigma2
            #theta2 <- yuima@data@mle@coef[yuima@model@parameter@drift]
            #theta1 <- yuima@data@mle@coef[yuima@model@parameter@diffusion]
            coef <- yuima@data@mle@coef
            
            Calpha2 <- qnorm(level)

            ## make matrix
            ret <- matrix(0, length(yuima@model@parameter@all), 2)
            rownames(ret) <- yuima@model@parameter@all
            a <- (1 - level)/2
            a <- c(a, 1 - a)
            colnames(ret) <- paste(round(100 * a, 1), "%")
            
            ## get confint
            Sigma.diag <- c( diag(Sigma2), diag(Sigma1) )
            names(Sigma.diag) <- yuima@model@parameter@all
            for(param in yuima@model@parameter@all){
              ret[param,] <- c( (coef[param]-sqrt(Sigma.diag[param])*Calpha2),
                              (coef[param]+sqrt(Sigma.diag[param])*Calpha2))
            }

            return(ret)
            
          })

##estimate theta2 by LSE
##function name LSE
setGeneric("LSE", function(yuima, h, theta2.init, interval, ...)
           standardGeneric("LSE")
           )
setMethod("LSE", "yuima",
          function(yuima, h, theta2.init=c(), interval=c(0,1), ...){
            ##theta2.init : multi-dim param used by optim()
            ##interval : 1-dim param used by optimize()

            ##objective function
            Mn <-function(theta2){
              Mn.part <- function(deltaX, h, yuima){
                tmp <- deltaX - h %*% eval(yuima@model@drift)
                ret <- t(tmp) %*% tmp
                return(ret)
              }
              
              ##init
              sum.tmp <- 0
              X <- NULL
              X.tmp <- get.zoo.data(yuima)
              for(i in 1:length(X.tmp)){
                X <- cbind(X, as.matrix(X.tmp[[i]]))
              }
              modelpara.drift <- yuima@model@parameter@drift
              modelstate <- yuima@model@state.variable
              for(i in 1:length(theta2)){
                assign(modelpara.drift[i], theta2[i])
              }
              
              ##sum loop
              for(j in 2:dim(X)[1]){
                ##get param
                x <- X[j-1,]
                deltaX <- X[j,] - X[j-1,]
                for(k in 1:length(x)){
                  assign(modelstate[k], x[k])
                }
                ##calc
                sum.tmp <- sum.tmp + Mn.part(deltaX, h, yuima)
              }
              
              #cat("theta2 value:", theta2, "\n")
              #cat("Mn value:", sum.tmp, "\n\n")
              
              return(sum.tmp)
  
            }
            ##END objective function

            if(length(yuima@model@parameter@drift)==1){

              if(missing(interval)){
                stop("\ninterval missing.\n")
              }
              if(!missing(theta2.init)){
                yuima.warn("theta2.init is ignored.")
              }
              
              opt <- optimize(f=Mn, interval=interval, tol=1e-100, ...)
              opt <- list(par=opt$minimum, value=opt$objective)

            }else{

              if(missing(theta2.init)){
                stop("\ntheta2.init missing.\n")
              }
              if(!missing(interval)){
                yuima.warn("interval is ignored.")
              }
              
              opt <- optim(c(theta2.init), fn=Mn, gr=NULL, ...)
              opt <- list(par=opt$par, value=opt$value)
              
            }            
            return(opt)            
          })
          



### I have rewritten qmle as a version of ml.ql
### This function has an interface more similar to mle.
### ml.ql is limited in that it uses fixed names for drift and diffusion
### parameters, while yuima model allows for any names.
### also, I am using the same interface of optim to specify upper and lower bounds
### S.M.I. 22/06/2010


qmle <- function(yuima, start, method="BFGS", fixed = list(), print=FALSE, ...){

	call <- match.call()
	
	if( missing(yuima))
		yuima.stop("yuima object is missing.")
	
## param handling
	if( missing(start) )
	 yuima.warn("Starting values for the parameters are missing. Using random initial values.")

	diff.par <- yuima@model@parameter@diffusion
	drift.par <- yuima@model@parameter@drift
	jump.par <- yuima@model@parameter@jump
	measure.par <- yuima@model@parameter@measure
	common.par <- yuima@model@parameter@common
	
	if(length(common.par)>0)
		yuima.stop("Drift and diffusion parameters must be different.")
	
	if(length(jump.par)+length(measure.par)>0)
		yuima.stop("Cannot estimate the jump models, yet")
	
	if(!is.list(start))
		yuima.stop("Argument 'start' must be of list type.")

	fullcoef <- c(diff.par, drift.par)
	npar <- length(fullcoef)
	
	
	fixed.par <- names(fixed)

	if (any(!(fixed.par %in% fullcoef))) 
	 yuima.stop("Some named arguments in 'fixed' are not arguments to the supplied yuima model")
	
	nm <- names(start)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
		yuima.stop("some named arguments in 'start' are not arguments to the supplied yuima model")
    start <- start[order(oo)]
    nm <- names(start)

	idx.diff <- match(diff.par, nm)
	idx.drift <- match(drift.par, nm)

	f <- function(p) {
        mycoef <- as.list(p)
        names(mycoef) <- nm
        mycoef[fixed.par] <- fixed
        negquasilik(yuima=yuima, param=mycoef, print=print)
    }
		
	oout <- if(length(start)){ 
		optim(start, f, method = method, hessian = TRUE, ...)
    } else {
		list(par = numeric(0L), value = f(start))
	}
	
	coef <- oout$par
    vcov <- if (length(coef)) 
	solve(oout$hessian)
    else matrix(numeric(0L), 0L, 0L)
    min <- oout$value
	
  	mycoef <- as.list(coef)
	names(mycoef) <- nm
	mycoef[fixed.par] <- fixed
	
    new("mle", call = call, coef = coef, fullcoef = unlist(mycoef), 
       vcov = vcov, min = min, details = oout, minuslogl = negquasilik, 
       method = method)
}

negquasilik <- function(yuima, param, print=FALSE){

	diff.par <- yuima@model@parameter@diffusion
	drift.par <- yuima@model@parameter@drift
	fullcoef <- c(diff.par, drift.par)
	npar <- length(fullcoef)
	
	nm <- names(param)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
		yuima.stop("some named arguments in 'param' are not arguments to the supplied yuima model")
    param <- param[order(oo)]
    nm <- names(param)
	
	idx.diff <- match(diff.par, nm)
	idx.drift <- match(drift.par, nm)

### FIXME	
	h <- deltat(yuima@data@zoo.data[[1]])

    theta1 <- unlist(param[idx.diff])
    theta2 <- unlist(param[idx.drift])
	n.theta1 <- length(theta1)
	n.theta2 <- length(theta2)
	n.theta <- n.theta1+n.theta2
	
	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	X <- as.matrix(onezoo(yuima))
	deltaX <- matrix(0, n-1, d.size)
	for(t in 1:(n-1))
	deltaX[t, ] <- X[t+1, ]-X[t, ]
	
	
	drift <- drift.term(yuima, param)
	diff <- diffusion.term(yuima, param)

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
		yuima.warn("quasi likelihood is too small to calculate.")
	}
	if(print==TRUE){

		yuima.warn(sprintf("NEG-QL: %f, %s", -QL, paste(names(param),param,sep="=",collapse=", ")))
	}

	return(-QL)

}



