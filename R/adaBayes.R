##::quasi-likelihood with prior
##
##

## rmnorm
rmnorm <- function (n, mean, cov=diag(length(mean))) 
{
	if(length(mean)==1){	
		return(rnorm(n,mean,sqrt(cov)))
	}else{
		return(rmvnorm(n,mean,cov))
	}
}

## dmnorm
dmnorm <- function (x, mean, cov=diag(length(mean))) 
{
	if(length(mean)==1){	
		return(dnorm(x,mean,sqrt(cov)))
	}else{
		return(dmvnorm(x,mean,cov))
	}
}

## ml.ql by newton method.
newton.ml.qlb <- function(yuima, theta2, theta1, h, iteration=1, param.only=FALSE, verbose=FALSE, ...){
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
	eval.inv.B <- function(B, theta1, theta2, Xt.iMinus1, ...){
    return(solve(eval.B(B, theta1, theta2, Xt.iMinus1, ...)))
	}
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
	
## Temporary dH2_theta1
	d2Htemp.theta1 <- function(theta1, theta2, h, yuima, B, a, ...){
		ret <- matrix(0, length(theta1), length(theta1))
		data <- as.matrix(onezoo(yuima))
		for( i in 1:(nrow(data)-1)){
##Xt.iMinus1
			Xt.iMinus1 <- data[i,]
##delta.i.X
			delta.i.X <- data[i+1,] - data[i,]
##calc		
			temp <- dGi.theta1(theta1, theta2, h, Xt.iMinus1, delta.i.X, B, a, ...)
			ret <- ret + temp%*%t(temp)
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
	
## Temp dH_theta2
	d2Htemp.theta2 <-function(theta1, theta2, h, yuima, B, a, ...){
		ret <- matrix(0, length(theta2), length(theta2))
		data <- as.matrix(onezoo(yuima))
		for( i in 1:(nrow(data)-1)){
##Xt.iMinus1
			Xt.iMinus1 <- data[i,]
##delta.i.X
			delta.i.X <- data[i+1,] - data[i,]
##calc
			temp <- dGi.theta2(theta1, theta2, h, Xt.iMinus1, delta.i.X, B, a, ...)
			ret <- ret + temp%*%t(temp)
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
			cat("theta2 init:", theta2, "\n")
			cat("theta1 init:", theta1, "\n")
		}
		
		for(ite in 1:iteration){
			
			dHtheta2 <- dH.theta2(theta1, theta2, h, yuima, B, a, ...)
			d2Htheta2 <- d2Htemp.theta2(theta1, theta2, h, yuima, B, a, ...)
			theta2 <- as.vector(theta2 - solve(d2Htheta2) %*% dHtheta2)
			dHtheta1 <- dH.theta1(theta1, theta2, h, yuima, B, a, ...)
			d2Htheta1 <- d2Htemp.theta1(theta1, theta2, h, yuima, B, a, ...)
			theta1 <- as.vector(theta1 - solve(d2Htheta1) %*% dHtheta1)
			if(verbose){
				cat("\n## Iteration", ite, "##\n")
				cat("theta2 new:", theta2, "\n")
				cat("theta1 new:", theta1, "\n")
			}
			
		}
	}
## END newtom algorithm ##
	
##calc Sigma1, Sigma2, Hessian?##
	Sigma1 <- - solve(d2H.theta1(theta1, theta2, h, yuima, B, a, ...))
	Sigma1temp <- - solve(d2Htemp.theta1(theta1, theta2, h, yuima, B, a, ...))
	Sigma2 <- - solve(d2H.theta2(theta1, theta2, h, yuima, B, a, ...))
    Sigma2temp <- - solve(d2Htemp.theta2(theta1, theta2, h, yuima, B, a, ...))
	hessian <- 1 ## Not Implemented yet!
##END 
##temp
	return(list(theta1.new=theta1, theta2.new=theta2, Sigma1=Sigma1, Sigma1temp=Sigma1temp, Sigma2=Sigma2,Sigma2temp=Sigma2temp, hessian=hessian))
	
}

##::calculate the log quasi-likelihood with parameters (theta2, theta1) and X.
##::yuima : yuima object
##::theta2 : parameter in drift.
##::theta1 : parameter in diffusion.
##::h : time width.

##::quasi-bayes function

##::estimate parameters(theta2,theta1) with a constraint ui%*%theta-ci=0  
##::yuima : yuima object
##::theta2 : init parameter in drift term.
##::theta1 : init parameter in diffusion term.
##::h : length between each observation time.
##::theta1.lim, theta2.lim : limit of those parameters.
##::example: 0 <= theta1 <= 1 theta1.lim = matrix(c(0,1),1,2)
##::if theta1, theta2 are matrix, theta.lim can be defined matrix like rbind(c(0,1),c(0,1))
setGeneric("adaBayes",
           function(yuima, print=FALSE, start, prior,propose,n.iter=100,lower,upper,n.burnin,method="nomcmc",mhtype="independent")
           standardGeneric("adaBayes")
           )
setMethod("adaBayes", "yuima",
          function(yuima, print=FALSE,  start, prior,propose,n.iter=100,lower,upper,n.burnin,method="nomcmc",mhtype="independent"){
            if( missing(yuima)){
              cat("\nyuima object is missing.\n")
              return(NULL)
            }
			            init <- start
			  if(length(match(yuima@model@parameter@drift ,names(init),nomatch=0))!=length(yuima@model@parameter@drift)){
				  cat("\ndrift parameters in yuima model and init do not match.\n")
				  return(NULL)
			  }
			  if(length(match(yuima@model@parameter@diffusion ,names(init),nomatch=0))!=length(yuima@model@parameter@diffusion)){
				  cat("\ndiffusion parameters in yuima model and init do not match.\n")
				  return(NULL)
			  }
			              
			## BEGIN burnin handling
			if(missing(n.burnin)){
              n.burnin <- n.iter
			}
			## END burnin handling
			  
			h <- deltat(yuima@data@zoo.data[[1]])
			
			  ## BEGIN Prior construction
			  liprior<- function(prior,term){
				  mvec <- numeric(0); mdom <- numeric(0)
				  for(i in 1:length(slot(yuima@model@parameter,term))){
					  if(prior[slot(yuima@model@parameter,term)[i]][[1]]$measure.type=="density"){
						  mvec <- append(mvec,prior[slot(yuima@model@parameter,term)[i]][[1]]$density)
						  if(is.null(prior[slot(yuima@model@parameter,term)[i]][[1]]$domain)){
							  mdom <- cbind(mdom,c(-Inf,Inf))
						  }else{
							  mdom <- cbind(mdom,prior[slot(yuima@model@parameter,term)[i]][[1]]$domain)
						  }
					  }
				  }
				  mdensity <- function(param){
					  res <- 1
					  for(i in 1:length(slot(yuima@model@parameter,term))){
						  res <- res*mvec[i][[1]](param[slot(yuima@model@parameter,term)[i]])
					  }
					  return(res)
					  
				  }
				  return(list(mdom=mdom,mdensity=mdensity))
			  }
			  ## END Prior construction
			  
			n <- length(yuima)[1]
			  
			  ilparam <- function(liparam){
				  return(list("theta2"=unlist(liparam[slot(yuima@model@parameter,"drift")]),"theta1"=unlist(liparam[slot(yuima@model@parameter,"diffusion")])))
			  }
			  liparam <- function(ilparam){
				  return(relist(unlist(ilparam),init))
			  }
			
			if(method=="nomcmc"){
			  require(adapt)
              ## BEGIN numerical integration function
              nintegral <- function(term,param=init,prior=prior,lower=lower,upper=upper,print=FALSE){
				lip <- liprior(prior,term); mdom <- lip$mdom; mdensity <- lip$mdensity
				  
				## BEGIN denominator calculation
				  nparam <- param
				  denominator <- function(subparam){
					  if(is.null(dim(subparam))){
						  matparam <- matrix(subparam,ncol=length(unlist(param[slot(yuima@model@parameter,term)])))
					  }else{
						  matparam <- subparam
					  }
					  
					  nparam <- param

					  qvec <- numeric(dim(matparam)[1])
					  for(k in 1:dim(matparam)[1]){
						  for(j in 1:(length(slot(yuima@model@parameter,term)))){
							  nparam[slot(yuima@model@parameter,term)][j] <- matparam[k,][grep(slot(yuima@model@parameter,term)[j],names(unlist(param[slot(yuima@model@parameter,term)])))]
						  }
						  qvec[k] <- exp(quasilogl(yuima,param=nparam,print=print)+log(mdensity(nparam))-quasilogl(yuima,param=param,print=print)-log(mdensity(param)))
					  }
					  return(qvec)
				  }
				  if(length(slot(yuima@model@parameter,term))==1){
					  denomvalue <- integrate(denominator,mdom[1,1],mdom[2,1])$value
				  }else{
					  denomvalue <- adapt(length(unlist(param[slot(yuima@model@parameter,term)])),lower=mdom[1,],upper=mdom[2,],functn=denominator)$value
				  }
				## END denominator calculation
				  
                ## BEGIN numerator calculation
				unlistparam <- numeric(length(param[slot(yuima@model@parameter,term)]))
				for(i in 1:length(param[slot(yuima@model@parameter,term)])){
					numevalue <- numeric(1)
					nparam <- param
					numerator <- function(subparam){
						if(is.null(dim(subparam))){
							matparam <- matrix(subparam,ncol=length(unlist(param[slot(yuima@model@parameter,term)])))
						}else{
							matparam <- subparam
						}
						qvec <- numeric(dim(matparam)[1])
						for(k in 1:(dim(matparam)[1])){
							nparam[slot(yuima@model@parameter,term)] <- matparam[k,]
							qvec[k] <- matparam[k,i]*exp(quasilogl(yuima,param=nparam,print=print)+log(mdensity(nparam))-quasilogl(yuima,param=param,print=print)-log(mdensity(param)))
						}
						return(qvec)
					}
					if(length(slot(yuima@model@parameter,term))==1){
						numevalue <- integrate(numerator,mdom[1,1],mdom[2,1])$value
					}else{
						numevalue <- adapt(length(param[slot(yuima@model@parameter,term)]),lower=mdom[1,],upper=mdom[2,],functn=numerator)$value
					}
					unlistparam[i] <- numevalue/denomvalue
				}
                ## END numerator calculation
				newparam <- param
				newparam[slot(yuima@model@parameter,term)] <- relist(unlistparam,param[slot(yuima@model@parameter,term)])
                
                return(newparam)
              }
              ## END numerical integration function
              
              ## BEGIN numerical integration procedure
              param <- nintegral("diffusion",param=init,prior=prior,lower=lower,upper=upper)
              param <- nintegral("drift",param=param,prior=prior,lower=lower,upper=upper)
              ## END numerical integration procedure
              
			}else if(method=="mcmc"){
              if(mhtype=="RW"){
                ## BEGIN propose construction
                if(missing(propose)){
                  rpropose1 <- function(n,mean,varcov){
                    sd <- sqrt(varcov)
                    d <- length(init[slot(yuima@model@parameter,"diffusion")])		
                    return(matrix(data=rnorm(d*n,mean,sd),nrow=n))
                  }
                  rpropose2 <- function(n,mean,varcov){
                    sd <- sqrt(varcov)
                    d <- length(init[slot(yuima@model@parameter,"drift")])
                    return(matrix(data=rnorm(d*n,mean,sd),nrow=n))
                  }

                  propose.param1 <- list(mean=0,varcov=1/n)
                  propose.param2 <- list(mean=0,varcov=1/(n*h))
                  propose1 <- list(rpropose=rpropose1,propose.param=propose.param1)
                  propose2 <- list(rpropose=rpropose2,propose.param=propose.param2)
                  propose <- list(propose2 = propose2,propose1 = propose1)
                }
                ## END propose construction
                
                ## BEGIN mh RW-algorithm function
                mhalgorithm <- function(j,length=n.iter,param=init,prior=prior,print=FALSE){
				if(prior[[j]]$measure.type=="density"){
					mdensity <- prior[[j]]$density
					if(is.null(prior[[j]]$domain)){
						mdom <- c(-Inf,Inf)
					}else{
						mdom <- prior[[j]]$domain
					}
				}

                  rpropose <- propose[[j]][[1]]
                  propose.param <- propose[[j]][[2]]
                  d <- length(param[[j]])
                  
                  dh <- matrix(t(rpropose(n.iter,propose.param[[1]],propose.param[[2]])),nrow=d)
                  
                  u <- runif(n.iter)
                  
                  pparam <- param
                  pql <- quasilogl(yuima,param=pparam,print=print)+log(mdensity(pparam[[j]]))
                  nparam <- param
                  mhestimator <- 0
                  
				for(i in 1:n.iter){
					nparam[[j]] <- pparam[[j]] + dh[,i]									
					nql <- quasilogl(yuima,param=nparam,print=print)+log(mdensity(nparam[[j]]))
					alpha <- exp(nql-pql)
					if(is.na(alpha)){alpha <- -Inf}
					pparam[[j]] <- (alpha>u[i])*(nparam[[j]]-pparam[[j]])+pparam[[j]]
					if(alpha>u[i]){ pql <- nql}
					mhestimator <- mhestimator + pparam[[j]]
					}
                  
                  return(mhestimator/n.iter)
                }
                ## END mh RW-algorithm function
              }
              
              if(mhtype=="independent"){
				require(mvtnorm)
                ## propose distribution construction
                if(missing(propose)){
                  rpropose1 <- function(n,mean,varcov){
                    d <- length(init[slot(yuima@model@parameter,"diffusion")])			
                    
                    if(d==1) varcov <- as.numeric(varcov)
                    
                    return( t( as.matrix(rmnorm(n,mean,varcov)) ) )
                  }
                  rpropose2 <- function(n,mean,varcov){
                    d <- length(init[slot(yuima@model@parameter,"drift")])		
                    
                    if(d==1) varcov <- as.numeric(varcov)
                    
                    return( t( as.matrix(rmnorm(n,mean,varcov)) ) )
                  }
                  dpropose <- function(x,mean,varcov){
                    n <- dim(as.matrix(x))[2]
                    d <- dim(as.matrix(x))[1]
                    
                    if(n==1 & d==1) x <- matrix(x)
                    
                    if(!is.matrix(x)){
                      x <- as.matrix(x)
                    }
                    dvector <- numeric(n)
					
                    for(i in 1:n){
                      dvector[i] <- dmnorm(x[,i],mean,varcov)
                    }
                    return(dvector)
                  }
                }
                
                ## BEGIN mh Independent type algorithm function
                mhalgorithm <- function(j,length=n.iter,param=init,prior=prior,print=FALSE){
				if(j==2){term <- "diffusion"}else{term <- "drift"}
					lip <- liprior(prior,term); mdom <- lip$mdom; mdensity <- lip$mdensity
				if(sum(c("theta2","theta1") %in% names(param))==2){
					unlistparam <- param
				}else{
					unlistparam <- ilparam(param)
				}
                  newton <- newton.ml.qlb(yuima, theta2=unlistparam[[1]],theta1=unlistparam[[2]], h, iteration=1)
                  if(print){ cat("newton result\n"); print(newton) }
					
                  if(j==1){##estimate theta2
                    newton.tmp <- newton.ml.qlb(yuima, theta2=unlist(newton$theta2.new), theta1=unlistparam[[2]], h, param.only=TRUE)
                    propose.param1 <- list(mean=newton.tmp$theta1.new,varcov=newton.tmp$Sigma1temp)
                    propose.param2 <- list(mean=newton.tmp$theta2.new,varcov=newton.tmp$Sigma2temp)
                    if(length(prior)!=1) prior <- prior$prior.theta2
                    
				  }else if(j==2){ ##estimate theta1
                    newton.tmp <- newton.ml.qlb(yuima, theta2=unlistparam[[1]], theta1=newton$theta1.new, h, param.only=TRUE)
                    propose.param1 <- list(mean=newton.tmp$theta1.new,varcov=newton.tmp$Sigma1temp)
                    propose.param2 <- list(mean=newton.tmp$theta2.new,varcov=newton.tmp$Sigma2temp)
                    if(length(prior)!=1) prior <- prior$prior.theta1
                  }
                  propose1 <- list(rpropose=rpropose1, dpropose=dpropose, propose.param=propose.param1)
                  propose2 <- list(rpropose=rpropose2, dpropose=dpropose, propose.param=propose.param2)
                  propose <- list(propose2=propose2, propose1=propose1)
                  param <- list(theta2=newton$theta2.new, theta1=newton$theta1.new)
					
                  rpropose <- propose[[j]][[1]]
                  dpropose <- propose[[j]][[2]]
                  propose.param <- propose[[j]][[3]]
                  d <- length(param[[j]])

                  state <- rpropose(n.iter,propose.param[[1]],propose.param[[2]])
                  tmp.param <- param
                  weight <- numeric(n.iter)

                  for(i in 1:n.iter){
                    tmp.param[[j]] <- state[,i] 
                    weight[i] <- exp(quasilogl(yuima,param=liparam(tmp.param),print=print))*mdensity(liparam(tmp.param))/
                      dpropose(state[,i],propose.param[[1]],propose.param[[2]])
                  }

                  cstate <- param[[j]]
                  cweight <- exp(quasilogl(yuima,param=liparam(param),print=print))*mdensity(liparam(tmp.param))/
                    dpropose(param[[j]],propose.param[[1]],as.matrix(propose.param[[2]]))

                  u <- runif(n.iter)
                  
                  mhestimator <- 0

                  for(i in 1:n.iter){
                    alpha <- weight[i]/cweight
                    if(is.na(alpha)){alpha <- -Inf}
                    if(alpha>u[i]){ cweight <- weight[i]; cstate <- state[,i]}
                    mhestimator <- mhestimator + state[,i]
                  }
                  return(mhestimator/n.iter)
                }
              }
              
              ## BEGIN mh procedure
              ##theta1 estim.
              param <- liparam(list(init[slot(yuima@model@parameter,"drift")],mhalgorithm(2,length=n.burnin,param=init,prior=prior, print=print)))
              param <- liparam(list(param[slot(yuima@model@parameter,"drift")],mhalgorithm(2,length=n.iter,param=param,prior=prior, print=print)))
              ##theta2 estim.
              param <- liparam(list(list(mhalgorithm(1,length=n.burnin,param=param,prior=prior, print=print),param[slot(yuima@model@parameter,"diffusion")])))
              param <- liparam(list(mhalgorithm(1,length=n.iter,param=param,prior=prior, print=print),param[slot(yuima@model@parameter,"diffusion")]))
              ## END mh procedure
			}
            return(liparam(param))
          })
