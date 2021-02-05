##::quasi-bayes function


setGeneric("adaBayes",
           function(yuima, start,prior,lower,upper, method="mcmc",iteration=NULL,mcmc,rate=1.0,
                    rcpp=TRUE,algorithm="randomwalk",center=NULL,sd=NULL,rho=NULL,path=FALSE
                  )
             standardGeneric("adaBayes")
)
# 
# adabayes<- setClass("adabayes",contains = "mle",
#                     slots = c(mcmc="list", accept_rate="list",coef = "numeirc",
#                                          call="call",vcov="matrix",fullcoef="numeric",model="yuima.model"))
adabayes<- setClass("adabayes",
                    #contains = "mle",
                    slots = c(mcmc="list", accept_rate="list",coef = "numeric",call="call",vcov="matrix",fullcoef="numeric"))

setMethod("adaBayes", 
          "yuima",
          function(yuima, start,prior,lower,upper, method="mcmc",iteration=NULL,mcmc,rate=1.0,rcpp=TRUE,
                   algorithm="randomwalk",center=NULL,sd=NULL,rho=NULL,path=FALSE
                   )
          {
            
            
            
            if(length(iteration)>0){mcmc=iteration}
            mean=unlist(center)
            joint <- FALSE
            fixed <- numeric(0)
            print <- FALSE
            
            call <- match.call()
            
            if( missing(yuima))
              yuima.stop("yuima object is missing.")
            
            ## param handling
            
            ## FIXME: maybe we should choose initial values at random within lower/upper
            ##        at present, qmle stops	
            
            if(missing(lower) || missing(upper)){
              if(missing(prior)){
                lower = numeric(0)
                upper = numeric(0)
                pdlist <- numeric(length(yuima@model@parameter@all))
                names(pdlist) <- yuima@model@parameter@all
                for(i in 1: length(pdlist)){
                  lower = append(lower,-Inf)
                  upper = append(upper,Inf)
                }
              }
              # yuima.stop("lower or upper is missing.")
            }
            
            diff.par <- yuima@model@parameter@diffusion
            drift.par <- yuima@model@parameter@drift
            jump.par <- yuima@model@parameter@jump
            measure.par <- yuima@model@parameter@measure
            common.par <- yuima@model@parameter@common

            
            ## BEGIN Prior construction
            if(!missing(prior)){
              priorLower = numeric(0)
              priorUpper = numeric(0)
              pdlist <- numeric(length(yuima@model@parameter@all))
              names(pdlist) <- yuima@model@parameter@all
              for(i in 1: length(pdlist)){
                if(prior[[names(pdlist)[i]]]$measure.type=="code"){
                  expr <- prior[[names(pdlist)[i]]]$df
                  code <- suppressWarnings(sub("^(.+?)\\(.+", "\\1", expr, perl=TRUE))
                  args <- unlist(strsplit(suppressWarnings(sub("^.+?\\((.+)\\)", "\\1", expr, perl=TRUE)), ","))
                  pdlist[i] <- switch(code,
                                      dunif=paste("function(z){return(dunif(z, ", args[2], ", ", args[3],"))}"),
                                      dnorm=paste("function(z){return(dnorm(z,", args[2], ", ", args[3], "))}"),
                                      dbeta=paste("function(z){return(dbeta(z, ", args[2], ", ", args[3], "))}"),
                                      dgamma=paste("function(z){return(dgamma(z, ", args[2], ", ", args[3], "))}"),
                                      dexp=paste("function(z){return(dexp(z, ", args[2], "))}")
                  )
                  qf <- switch(code,
                               dunif=paste("function(z){return(qunif(z, ", args[2], ", ", args[3],"))}"),
                               dnorm=paste("function(z){return(qnorm(z,", args[2], ", ", args[3], "))}"),
                               dbeta=paste("function(z){return(qbeta(z, ", args[2], ", ", args[3], "))}"),
                               dgamma=paste("function(z){return(qgamma(z, ", args[2], ", ", args[3], "))}"),
                               dexp=paste("function(z){return(qexp(z, ", args[2], "))}")
                  )
                  priorLower = append(priorLower,eval(parse("text"=qf))(0.00))
                  priorUpper = append(priorUpper,eval(parse("text"=qf))(1.00))
                  
                  
                }
                
              }
              
              #20200518kaino
              names(priorLower)<-names(pdlist)
              names(priorUpper)<-names(pdlist)
              
              if(missing(lower) || missing(upper)){
                lower <- priorLower
                upper <- priorUpper
              }
              else{
                #20200518kaino
                #if(sum(unlist(lower)<priorLower) + sum(unlist(upper)>priorUpper) > 0){
                if(sum(unlist(lower)<priorLower[names(lower)]) + sum(unlist(upper)>priorUpper[names(upper)]) > 0){
                  yuima.stop("lower&upper of prior are out of parameter space.")
                }
              }
              
              #20200518
              #names(lower) <- names(pdlist)
              #names(upper) <- names(pdlist)
              
              
              
              pd <- function(param){
                value <- 1
                for(i in 1:length(pdlist)){
                  value <- value*eval(parse(text=pdlist[[i]]))(param[[i]])
                }
                return(value)
              }
            }else{
              pd <- function(param) return(1)
            }
            ## END Prior construction
            
            if(!is.list(lower)){
              lower <- as.list(lower)
            }
            if(!is.list(upper)){
              upper <- as.list(upper)
            }
            
            JointOptim <- joint
            if(length(common.par)>0){
              JointOptim <- TRUE
              yuima.warn("Drift and diffusion parameters must be different. Doing
                         joint estimation, asymptotic theory may not hold true.")
            }


            if(length(jump.par)+length(measure.par)>0)
              yuima.stop("Cannot estimate the jump models, yet")
            
            
            fulcoef <- NULL
            
            if(length(diff.par)>0)
              fulcoef <- diff.par
            
            if(length(drift.par)>0)
              fulcoef <- c(fulcoef, drift.par)
            
            if(length(yuima@model@parameter@common)>0){
              fulcoef<-unique(fulcoef)
            }
            
            
            fixed.par <- names(fixed)
            
            if (any(!(fixed.par %in% fulcoef))) 
              yuima.stop("Some named arguments in 'fixed' are not arguments to the supplied yuima model")
            
            nm <- names(start)
            npar <- length(nm)

            oo <- match(nm, fulcoef)
            if(any(is.na(oo))) 
              yuima.stop("some named arguments in 'start' are not arguments to the supplied yuima model")
            start <- start[order(oo)]
            if(!missing(prior)){
              #20200525kaino
              #pdlist <- pdlist[order(oo)]
              pdlist <- pdlist[names(start)]
            }
            nm <- names(start)
            
            idx.diff <- match(diff.par, nm)
            idx.drift <- match(drift.par,nm)
            if(length(common.par)>0){
              idx.common=match(common.par,nm)
              idx.drift=setdiff(idx.drift,idx.common)
            }
            idx.fixed <- match(fixed.par, nm)
            tmplower <- as.list( rep( -Inf, length(nm)))
            names(tmplower) <- nm	
            if(!missing(lower)){
              idx <- match(names(lower), names(tmplower))
              if(any(is.na(idx)))
                yuima.stop("names in 'lower' do not match names fo parameters")
              tmplower[ idx ] <- lower	
            }
            lower <- tmplower
            
            tmpupper <- as.list( rep( Inf, length(nm)))
            names(tmpupper) <- nm	
            if(!missing(upper)){
              idx <- match(names(upper), names(tmpupper))
              if(any(is.na(idx)))
                yuima.stop("names in 'lower' do not match names fo parameters")
              tmpupper[ idx ] <- upper	
            }
            upper <- tmpupper
            
            
            
            
            d.size <- yuima@model@equation.number
            #20200601kaino
            if (is.CARMA(yuima)){
              # 24/12
              d.size <-1
            }
            n <- length(yuima)[1]
            
            G <- rate
            if(G<=0 || G>1){
              yuima.stop("rate G should be 0 < G <= 1")
            }
            n_0 <- floor(n^G)
            if(n_0 < 2) n_0 <- 2
            
            #######data is reduced to n_0 before qmle(16/11/2016)
            env <- new.env()
            #assign("X",  yuima@data@original.data[1:n_0,], envir=env)
            assign("X",  as.matrix(onezoo(yuima)[1:n_0,]), envir=env)
            assign("deltaX",  matrix(0, n_0 - 1, d.size), envir=env)
            assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env)
            #20200601kaino
            if (is.CARMA(yuima)){
              #24/12 If we consider a carma model,
              # the observations are only the first column of env$X
              #     assign("X",  as.matrix(onezoo(yuima)), envir=env)
              #     env$X<-as.matrix(env$X[,1])
              #     assign("deltaX",  matrix(0, n-1, d.size)[,1], envir=env)
              env$X<-as.matrix(env$X[,1])
              env$deltaX<-as.matrix(env$deltaX[,1])
              assign("time.obs",length(env$X),envir=env)
              assign("p", yuima@model@info@p, envir=env)
              assign("q", yuima@model@info@q, envir=env)
              assign("V_inf0", matrix(diag(rep(1,env$p)),env$p,env$p), envir=env)
            }
            
            assign("Cn.r", rep(1,n_0-1), envir=env)
            
            for(t in 1:(n_0-1))
              env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
            
            assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)
            
            pp<-0 
            if(env$h >= 1){
              qq <- 2/G
            }else{
              while(1){
                if(n*env$h^pp < 0.1) break
                pp <- pp + 1
              }
              qq <- max(pp,2/G) 
            }
            
            C.temper.diff <- n_0^(2/(qq*G)-1) #this is used in pg.
            C.temper.drift <- (n_0*env$h)^(2/(qq*G)-1) #this is used in pg.
            
            mle <- qmle(yuima, "start"=start, "lower"=lower,"upper"=upper, "method"="L-BFGS-B",rcpp=rcpp)
            start <- as.list(mle@coef)

            sigma.diff=NULL
            sigma.drift=NULL
            
            if(length(sd)<npar){
              sigma=mle@vcov;
            }
              # if(length(diff.par)>0){
              #   sigma.diff=diag(1,length(diff.par))
              #   for(i in 1:length(diff.par)){
              #     sigma.diff[i,i]=sd[[diff.par[i]]]
              #   }
              # }
              # 
              # if(length(drift.par)>0){
              #   sigma.drift=diag(1,length(drift.par))
              #   for(i in 1:length(drift.par)){
              #     sigma.drift[i,i]=sd[[drift.par[i]]]
              #   }
              # }
            #}
            
            
            # mu.diff=NULL
            # mu.drift=NULL
            # 
            # 
            # if(length(mean)<npar){
            #   mu.diff=NULL
            #   mu.drift=NULL}
            #else{
            #   if(length(diff.par)>0){
            #   for(i in 1:length(diff.par)){
            #     mu.diff[i]=mean[[diff.par[i]]]
            #   }
            #   }
            #   
            #   if(length(drift.par)>0){
            #   for(i in 1:length(drift.par)){
            #     mu.drift[i]=mean[[drift.par[i]]]
            #   }
            #   }
            # }
            
            ####mpcn make proposal
            sqn<-function(x,LL){
              vv=x%*%LL
              zz=sum(vv*vv)
              if(zz<0.0000001){
                zz=0.0000001
              }
              return(zz)
            }
            
            # mpro<-function(mu,sample,low,up,rho,LL,L){
            #   d=length(mu);
            #   tmp=mu+sqrt(rho)*(sample-mu);
            #   tmp = mu+sqrt(rho)*(sample-mu);
            #   dt=(sample-mu)%*%LL%*%(sample-mu)
            #   tmp2 = 2.0/dt;
            #   while(1){
            #     prop = tmp+rnorm(d)*L*sqrt((1.0-rho)/rgamma(1,0.5*d,tmp2));
            #     if((sum(low>prop)+sum(up<prop))==0) break;
            #   }
            #   return(prop)
            # }

            make<-function(mu,sample,low,up,rho,LL,L){
              
              d=length(mu);
              tmp=mu+sqrt(rho)*(sample-mu);
              dt=sqn(sample-mu,LL);
              tmp2 = 2.0/dt;
    
              prop = tmp+rnorm(d)%*%L*sqrt((1.0-rho)/rgamma(1,0.5*d,scale=tmp2));
              
              # for(i in 1:100){
              #   prop = tmp+rnorm(d)*sqrt((1.0-rho)/rgamma(1,0.5*d,tmp2));
              #   if((sum(low>prop)+sum(up<prop))==0){break;}
              #   flg=flg+1;
              # }
              # if(flg>100){print("error")}else{
              #   return(prop)
              # }
              return(prop)
            }
            
            
            integ <- function(idx.fixed=NULL,f=f,start=start,par=NULL,hessian=FALSE,upper,lower){
              if(length(idx.fixed)==0){
                intf <- hcubature(f,lowerLimit=unlist(lower),upperLimit=unlist(upper),fDim=(length(upper)+1))$integral
              }else{
                intf <- hcubature(f,lowerLimit=unlist(lower[-idx.fixed]),upperLimit=unlist(upper[-idx.fixed]),fDim=(length(upper[-idx.fixed])+1))$integral
              }
              return(intf[-1]/intf[1])
            }
            mcinteg <- function(idx.fixed=NULL,f=f,p,start=start,par=NULL,hessian=FALSE,upper,lower,mean,vcov,mcmc){
              #vcov=vcov[nm,nm];#Song add
              if(setequal(unlist(drift.par),unlist(diff.par))&(length(mean)>length(diff.par))){mean=mean[1:length(diff.par)]}
              if(length(idx.fixed)==0){#only have drift or diffusion part
                intf <- mcintegrate(f,p,lowerLimit=lower,upperLimit=upper,mean,vcov,mcmc)
              }else{
                intf <- mcintegrate(f,p,lowerLimit=lower[-idx.fixed],upperLimit=upper[-idx.fixed],mean[-idx.fixed],vcov[-idx.fixed,-idx.fixed],mcmc)
              }
              return(intf)
            }
            
            
            mcintegrate <- function(f,p, lowerLimit, upperLimit,mean,vcov,mcmc){
              
              acc=0
              if(algorithm=="randomwalk"){
                x_c <- mean
                if(path){
                  X <- matrix(0,mcmc,length(mean))
                  acc=0
                }
              
                nm=length(mean)
                if(length(vcov)!=(nm^2)){
                  vcov=diag(1,nm,nm)*vcov[1:nm,1:nm]
                }
                
                sigma=NULL
                if(length(sd)>0&length(sd)==npar){
                    sigma=diag(sd,npar,npar);
                    if(length(idx.fixed)>0){
                      vcov=sigma[-idx.fixed,-idx.fixed]
                    }else{
                      vcov=sigma;
                    }
                  }
                
                p_c <- p(x_c)
                val <- f(x_c)
                if(path) X[1,] <- x_c
                # if(length(mean)>1){
                #   x <- rmvnorm(mcmc-1,mean,vcov)
                #   q <- dmvnorm(x,mean,vcov)
                #   q_c <- dmvnorm(mean,mean,vcov) 
                # }else{
                #   x <- rnorm(mcmc-1,mean,sqrt(vcov))
                #   q <- dnorm(x,mean,sqrt(vcov))
                #   q_c <- dnorm(mean,mean,sqrt(vcov)) 
                # }
                #
                for(i in 1:(mcmc-1)){
                  if(length(mean)>1){
                    
                      x_n <- rmvnorm(1,x_c,vcov)
                      #if(sum(x_n<lowerLimit)==0 & sum(x_n>upperLimit)==0) break
                  }else{
                      x_n <- rnorm(1,x_c,sqrt(vcov))
                      #if(sum(x_n<lowerLimit)==0 & sum(x_n>upperLimit)==0) break
                  }
                  
                  if(sum(x_n<lowerLimit)==0 & sum(x_n>upperLimit)==0){
                    p_n <- p(x_n)
                    u <- log(runif(1))
                    a <- p_n-p_c
                    if(u<a){
                      p_c <- p_n
                      # q_c <- q_n
                      x_c <- x_n
                      acc=acc+1
                    }
                  }
                  
                  if(path) X[i+1,] <- x_c
                  #}
                  val <- val+f(x_c)
                }
                if(path){
                  return(list(val=unlist(val/mcmc),X=X,acc=acc/mcmc))
                }else{
                  return(list(val=unlist(val/mcmc)))
                } 
              }
              else if(tolower(algorithm)=="mpcn"){ #MpCN
                if(path) X <- matrix(0,mcmc,length(mean))
                if((sum(idx.diff==idx.fixed)==0)){
                  if(length(center)>0){
                    if(length(idx.fixed)>0){
                      mean =unlist(center[-idx.fixed])
                    }else{
                      mean =unlist(center) ##drift or diffusion parameter only
                    }
                  }
                }
                
                # if((sum(idx.diff==idx.fixed)>0)){
                #   if(length(center)>0){
                #     mean =unlist(center[-idx.fixed])
                #   }
                # }
                x_n <- mean
                val <- mean
                
                L=diag(1,length(mean));LL=L;
                nm=length(mean)
                LL=vcov;
                if(length(vcov)!=(nm^2)){
                  LL=diag(1,nm,nm)
                }
                if(length(sd)>0&length(sd)==npar){
                  sigma=diag(sd,npar,npar);
                  if(length(idx.fixed)>0){
                    LL=sigma[-idx.fixed,-idx.fixed]
                  }else{
                    LL=sigma;
                  }
                }
                
                
                # if(length(pre_conditionnal)==1){
                #   LL=t(chol(solve(vcov)));L=chol(vcov);
                # }
                if(length(rho)==0){rho=0.8}
                
                logLik_old <- p(x_n)+0.5*length(mean)*log(sqn(x_n-mean,LL))
                if(path) X[1,] <- x_n
                
                
                for(i in 1:(mcmc-1)){
                  #browser()
                  prop <- make(mean,x_n,unlist(lowerLimit),unlist(upperLimit),rho,LL,L)
                  if((sum(prop<lowerLimit)+sum(prop>upperLimit))==0){
                  logLik_new <- p(prop)+0.5*length(mean)*log(sqn(prop-mean,LL))
                  u <- log(runif(1))
                  if( logLik_new-logLik_old > u){
                    x_n <- prop
                    logLik_old <- logLik_new
                    acc=acc+1
                  }
                  }
                  if(path) X[i+1,] <- x_n
                  val <- val+f(x_n)
                }
                #return(unlist(val/mcmc))
                if(path){
                  return(list(val=unlist(val/mcmc),X=X,acc=acc/mcmc))
                }else{
                  return(list(val=unlist(val/mcmc)))
                }
              }
            }
            
            #print(mle@coef)
            flagNotPosDif <- 0
            
              for(i in 1:npar){
                if(mle@vcov[i,i] <= 0) flagNotPosDif <- 1 #Check mle@vcov is positive difinite matrix
              }
              if(flagNotPosDif == 1){
                mle@vcov <- diag(c(rep(1 / n_0,length(diff.par)),rep(1 / (n_0 * env$h),length(npar)-length(diff.par)))) # Redifine mle@vcov
              }
            
              # cov.diff=mle@vcov[1:length(diff.par),1:length(diff.par)]
              # cov.drift=mle@vcov[(length(diff.par)+1):(length(diff.par)+length(drift.par)),(length(diff.par)+1):(length(diff.par)+length(drift.par))]
              # if(missing(cov){
              # }else{
              # cov.diff=cov[[1]]
              # cov.drift=cov[[2]]
            
            
            
            tmp <- minusquasilogl(yuima=yuima, param=mle@coef, print=print, env,rcpp=rcpp)
            
            g <- function(p,fixed,idx.fixed){
              mycoef <- mle@coef
              #mycoef <- as.list(p)
              if(length(idx.fixed)>0){
                mycoef[-idx.fixed] <- p
                mycoef[idx.fixed] <- fixed
              }else{
                mycoef <- as.list(p)
                names(mycoef) <- nm
              }
              return(c(1,p)*exp(-minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmp)*pd(param=mycoef))
            }
            
            pg <- function(p,fixed,idx.fixed){
              mycoef <- start
              #mycoef <- as.list(p)
              if(length(idx.fixed)>0){
                mycoef[-idx.fixed] <- p
                mycoef[idx.fixed] <- fixed
              }else{
                mycoef <- as.list(p)
                names(mycoef) <- nm
              }
              
              # mycoef <- as.list(p)
              # names(mycoef) <- nm
              #return(exp(-minusquasilogl(yuima=yuima, param=mycoef, print=print, env)+tmp)*pd(param=mycoef))
              #return(-minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmp+log(pd(param=mycoef)))#log
              if(sum(idx.diff==idx.fixed)>0){
                #return(C.temper.diff*(-minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmp+log(pd(param=mycoef))))#log
                return(C.temper.drift*(-minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmp+log(pd(param=mycoef))))#log
              }else{
                #return(C.temper.drift*(-minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmp+log(pd(param=mycoef))))#log
                return(C.temper.diff*(-minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)+tmp+log(pd(param=mycoef))))#log
              }
            }
            
            idf <- function(p){return(p)}
            
            #	 fj <- function(p) {
            #		 mycoef <- as.list(p)
            #		 names(mycoef) <- nm
            #		 mycoef[fixed.par] <- fixed
            #		 minusquasilogl(yuima=yuima, param=mycoef, print=print, env)
            #	 }
            
            X.diff <- NULL
            X.drift <- NULL
            
            oout <- NULL
            HESS <- matrix(0, length(nm), length(nm))
            colnames(HESS) <- nm
            rownames(HESS) <- nm
            HaveDriftHess <- FALSE
            HaveDiffHess <- FALSE
            if(length(start)){
              #		if(JointOptim){ ### joint optimization
              #			if(length(start)>1){ #multidimensional optim
              #				oout <- optim(start, fj, method = method, hessian = TRUE, lower=lower, upper=upper)
              #				HESS <- oout$hessian
              #				HaveDriftHess <- TRUE
              #				HaveDiffHess <- TRUE
              #			} else { ### one dimensional optim
              #				opt1 <- optimize(f, ...) ## an interval should be provided
              #				opt1 <- list(par=integ(f=f,upper=upper,lower=lower,fDim=length(lower)+1),objective=0)
              #               oout <- list(par = opt1$minimum, value = opt1$objective)
              #			} ### endif( length(start)>1 )
              #		} else {  ### first diffusion, then drift
              theta1 <- NULL
              acc.diff<-NULL
              
              old.fixed <- fixed 
              old.start <- start
              
              if(length(idx.diff)>0){
                ## DIFFUSION ESTIMATIOn first
                old.fixed <- fixed
                old.start <- start
                new.start <- start[idx.diff] # considering only initial guess for diffusion
                new.fixed <- fixed
                if(length(idx.drift)>0)	
                  new.fixed[nm[idx.drift]] <- start[idx.drift]
                fixed <- new.fixed
                fixed.par <- names(fixed)
                idx.fixed <- match(fixed.par, nm)
                names(new.start) <- nm[idx.diff]
                
                f <- function(p){return(g(p,fixed,idx.fixed))}
                pf <- function(p){return(pg(p,fixed,idx.fixed))}
                if(length(unlist(new.start))>1){
                  #			 oout <- do.call(optim, args=mydots)
                  if(method=="mcmc"){
                    #oout <- list(par=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=mle@vcov,mcmc=mcmc))
                    mci=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=mle@vcov,mcmc=mcmc)
                    X.diff=mci$X
                    acc.diff=mci$acc
                    oout <- list(par=mci$val)
                  }else{
                    oout <- list(par=integ(idx.fixed=idx.fixed,f=f,upper=upper,lower=lower,start=start))
                  }
                } else {
                  #			 opt1 <- do.call(optimize, args=mydots)
                  if(method=="mcmc"){
                    #opt1 <- list(minimum=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=mle@vcov,mcmc=mcmc))
                    mci=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=mle@vcov,mcmc=mcmc)
                    if(path) X.diff=mci$X
                    acc.diff=mci$acc
                    opt1 <- list(minimum=mci$val)
                  }else{
                    opt1 <- list(minimum=integ(idx.fixed=idx.fixed,f=f,upper=upper,lower=lower))
                  }
                  theta1 <- opt1$minimum
                  #names(theta1) <- diff.par
                  #			 oout <- list(par = theta1, value = opt1$objective) 
                  oout <- list(par=theta1,value=0)
                }
                theta1 <- oout$par
                #names(theta1) <- nm[idx.diff]
                names(theta1) <- diff.par
              } ## endif(length(idx.diff)>0)
              
              theta2 <- NULL
              acc.drift<-NULL
              
              if(length(idx.drift)>0 & setequal(unlist(drift.par),unlist(diff.par))==FALSE){
                ## DRIFT estimation with first state diffusion estimates
                fixed <- old.fixed
                start <- old.start
                new.start <- start[idx.drift] # considering only initial guess for drift
                new.fixed <- fixed
                new.fixed[names(theta1)] <- theta1
                fixed <- new.fixed
                fixed.par <- names(fixed)
                idx.fixed <- match(fixed.par, nm)
                names(new.start) <- nm[idx.drift]
                
                f <- function(p){return(g(p,fixed,idx.fixed))}
                pf <- function(p){return(pg(p,fixed,idx.fixed))}
                
                if(length(unlist(new.start))>1){
                  #			  oout1 <- do.call(optim, args=mydots)
                  if(method=="mcmc"){
                    #oout1 <- list(par=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=mle@vcov,mcmc=mcmc))
                    mci=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=mle@vcov,mcmc=mcmc)
                    if(path) X.drift=mci$X
                    acc.drift=mci$acc
                    #20200601kaino
                    #oout1 <- list(minimum=mci$val)
                    oout1 <- list(par=mci$val)
                  }else{
                    oout1 <- list(par=integ(idx.fixed=idx.fixed,f=f,upper=upper,lower=lower))
                  }
                } else {
                  #				opt1 <- do.call(optimize, args=mydots)
                  if(method=="mcmc"){
                    #opt1 <- list(minimum=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=mle@vcov,mcmc=mcmc))
                    mci=mcinteg(idx.fixed=idx.fixed,f=idf,p=pf,upper=upper,lower=lower,mean=mle@coef,vcov=mle@vcov,mcmc=mcmc)
                    X.drift=mci$X
                    acc.drift=mci$acc
                    opt1 <- list(minimum=mci$val)
                  }else{
                    opt1 <- list(minimum=integ(idx.fixed=idx.fixed,f=f,upper=upper,lower=lower))
                  }
                  theta2 <- opt1$minimum
                  names(theta2) <-nm[-idx.fixed]
                  oout1 <- list(par = theta2, value = as.numeric(opt1$objective)) 	
                }
                theta2 <- oout1$par
              } ## endif(length(idx.drift)>0)
              oout1 <- list(par=c(theta1, theta2))
              names(oout1$par) <-nm # c(diff.par,drift.par)
              oout <- oout1
              
              #		} ### endif JointOptim
            } else {
              list(par = numeric(0L), value = f(start))
            }
            
            
            
            # if(path){
            #   return(list(coef=oout$par,X.diff=X.diff,X.drift=X.drift,accept_rate=accept_rate))
            # }else{
            #   return(list(coef=oout$par))
            # }
            fDrift <- function(p) {
              mycoef <- as.list(p)
              names(mycoef) <- drift.par
              mycoef[diff.par] <- coef[diff.par]
              minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)
            }
            
            fDiff <- function(p) {
              mycoef <- as.list(p)
              names(mycoef) <- diff.par
              mycoef[drift.par] <- coef[drift.par]
              minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)
            }
            
            coef <- oout$par
            control=list()
            par <- coef
            names(par) <- nm
            #names(par) <- c(diff.par, drift.par)
            #nm <- c(diff.par, drift.par)
            
            
            
            #	 print(par)
            #	 print(coef)
            conDrift <- list(trace = 5, fnscale = 1, 
                             parscale = rep.int(5, length(drift.par)), 
                             ndeps = rep.int(0.001, length(drift.par)), maxit = 100L, 
                             abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
                             beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
                             factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
            conDiff <- list(trace = 5, fnscale = 1, 
                            parscale = rep.int(5, length(diff.par)), 
                            ndeps = rep.int(0.001, length(diff.par)), maxit = 100L, 
                            abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
                            beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
                            factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
            
            #	 nmsC <- names(con)
            #	 if (method == "Nelder-Mead") 
            #	 con$maxit <- 500
            #	 if (method == "SANN") {
            #		 con$maxit <- 10000
            #		 con$REPORT <- 100
            #	 }
            #	 con[(namc <- names(control))] <- control
            #	 if (length(noNms <- namc[!namc %in% nmsC])) 
            #	 warning("unknown names in control: ", paste(noNms, collapse = ", "))
            #	 if (con$trace < 0) 
            #	 warning("read the documentation for 'trace' more carefully")
            #	 else if (method == "SANN" && con$trace && as.integer(con$REPORT) == 
            #			  0) 
            #	 stop("'trace != 0' needs 'REPORT >= 1'")
            #	 if (method == "L-BFGS-B" && any(!is.na(match(c("reltol", 
            #													"abstol"), namc)))) 
            #	 warning("method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'")
            #	 npar <- length(par)
            #	 if (npar == 1 && method == "Nelder-Mead") 
            #	 warning("one-diml optimization by Nelder-Mead is unreliable: use optimize")
            #	 
            if(!HaveDriftHess & (length(drift.par)>0)){
              #hess2 <- .Internal(optimhess(coef[drift.par], fDrift, NULL, conDrift))
              hess2 <- optimHess(coef[drift.par], fDrift, NULL, control=conDrift)
              HESS[drift.par,drift.par] <- hess2	 
            }
            
            if(!HaveDiffHess  & (length(diff.par)>0)){
              #hess1 <- .Internal(optimhess(coef[diff.par], fDiff, NULL, conDiff))
              hess1 <- optimHess(coef[diff.par], fDiff, NULL, control=conDiff)
              HESS[diff.par,diff.par] <- hess1	 
            }
            
            oout$hessian <- HESS
            
            # vcov <- if (length(coef)) 
            #   solve(oout$hessian)
            # else matrix(numeric(0L), 0L, 0L)
            
            mycoef <- as.list(coef)
            names(mycoef) <- nm
            mycoef[fixed.par] <- fixed
            
            #min <- minusquasilogl(yuima=yuima, param=mycoef, print=print, env,rcpp=rcpp)
            # 
            
            # new("mle", call = call, coef = coef, fullcoef = unlist(mycoef), 
            #     
            #     #       vcov = vcov, min = min, details = oout, minuslogl = minusquasilogl, 
            #     vcov = vcov,  details = oout, 
            #     method = method)
            
            mcmc_sample<-cbind(X.diff,X.drift)
            
            mcmc<-list()
            
            lf=length(diff.par)
            vcov=matrix(0,npar,npar)
            if(path){
              if(lf>1){
                vcov[1:lf,1:lf]=cov(X.diff)
              }else if(lf==1){
                vcov[1:lf,1:lf]=var(X.diff)
              }
              
              if((npar-lf)>1){
                vcov[(lf+1):npar,(lf+1):npar]=cov(X.drift)
              }else if((npar-lf)==1){
                vcov[(lf+1):npar,(lf+1):npar]=var(X.drift)
              }
            }
            
            
            
            accept_rate<-list()
            if(path){
              for(i in 1:npar){
                mcmc[[i]]=as.mcmc(mcmc_sample[,i])
              }
              #para_drift=yuima@model@parameter@drift[yuima@model@parameter@drift!=yuima@model@parameter@common]
              names(mcmc) <-nm;# c(yuima@model@parameter@diffusion,para_drift)
              accept_rate=list(acc.drift,acc.diff)
              names(accept_rate)<-c("accepte.rate.drift","accepte.rate.diffusion")
            }else{
              mcmc=list("NULL")
              accept_rate=list("NULL")
            }
            
            fulcoef = unlist(mycoef);
   
           new("adabayes",mcmc=mcmc, accept_rate=accept_rate,coef=fulcoef,call=call,vcov=vcov,fullcoef=fulcoef)
          }
)

setGeneric("coef")
setMethod("coef", "adabayes", function(object) object@fullcoef )

