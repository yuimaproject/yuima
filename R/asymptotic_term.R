
yuima.warn <- function(x){
	cat(sprintf("\nYUIMA: %s\n", x))
}

# in this source we note formulae like latex


setGeneric("asymptotic_term",
           function(yuima,block=100, rho, g, expand.var="e")
           standardGeneric("asymptotic_term")
           )

setMethod("asymptotic_term",signature(yuima="yuima"), function(yuima,block=100, rho, g, expand.var="e"){

  if(is.null(yuima@model)) stop("model object is missing!")
  if(is.null(yuima@sampling)) stop("sampling object is missing!")
  if(is.null(yuima@functional)) stop("functional object is missing!")

  f <- getf(yuima@functional)
  F <- getF(yuima@functional)

  ##:: fix bug 07/23
  #e <- gete(yuima@functional)
  assign(expand.var, gete(yuima@functional))
  
##  Terminal <- yuima@sampling@Terminal
  Terminal <- yuima@sampling@Terminal[1]
##  division <- yuima@sampling@n
  division <- yuima@sampling@n[1]
  xinit <- getxinit(yuima@functional)
  state <- yuima@model@state.variable
  V0 <- yuima@model@drift
  V <- yuima@model@diffusion
  r.size <- yuima@model@noise.number
  d.size <- yuima@model@equation.number
  k.size <- length(F)

  print("compute X.t0")
  X.t0 <- Get.X.t0(yuima, expand.var=expand.var)
  delta <- deltat(X.t0)

  ##:: fix bug 07/23
  pars <- expand.var #yuima@model@parameter@all[1]  #epsilon
	
	# fix bug 20110628
	if(k.size==1){
		G <- function(x){
			n <- length(x)
			res <- numeric(0)
			for(i in 1:n){
				res <- append(res,g(x[i]))
			}
			return(res)
		}
	}else{
		G <- function(x) return(g(x))
	}


  # function to return symbolic derivatives of myfunc by mystate(multi-state)
  Derivation.vector <- function(myfunc,mystate,dim1,dim2){
    tmp <- vector(dim1*dim2,mode="expression")
    for(i in 1:dim1){
      for(j in 1:dim2){
        tmp[(i-1)*dim2+j] <- parse(text=deparse(D(myfunc[i],mystate[j])))
      }
    }
    return(tmp)
  }
  
  # function to return symbolic derivatives of myfunc by mystate(single state)
  Derivation.scalar <- function(myfunc,mystate,dim){
    tmp <- vector(dim,mode="expression")
    for(i in 1:dim){
      tmp[i] <- parse(text=deparse(D(myfunc[i],mystate)))
    }
    return(tmp)
  }

  # function to solve Y_{t} (between (13.9) and (13.10)) using runge kutta method. Y_{t} is GL(d) valued (matrices)
  Get.Y <- function(){
    ## init
    dt <- Terminal/division
    assign(pars,0)  ## epsilon=0
    Yinit <- as.vector(diag(d.size))
    Yt <- Yinit
    Y <- Yinit
    k <- numeric(d.size*d.size)
    k1 <- numeric(d.size*d.size)
    k2 <- numeric(d.size*d.size)
    k3 <- numeric(d.size*d.size)
    k4 <- numeric(d.size*d.size)
    Ystate <- paste("y",1:(d.size*d.size),sep="")

    F <- NULL
    F.n <- vector(d.size,mode="expression")
    for(n in 1:d.size){
      for(i in 1:d.size){
        F.tmp <- dx.drift[((i-1)*d.size+1):(i*d.size)]

##C³
#        F.n[i] <- parse(text=paste(paste(F.tmp,"*",Ystate[((n-1)*d.size+1):(n*d.size)],sep=""),collapse="+"))
        F.n[i] <- parse(text=paste(paste("(",F.tmp,")","*",Ystate[((n-1)*d.size+1):(n*d.size)],sep=""),collapse="+"))
##

      }
      F <- append(F,F.n)
    }
    ## runge kutta
    for(t in 1:division){
      ## Xt
      for(i in 1:d.size){
        assign(state[i],X.t0[t,i])   ## state[i] is x_i, for example.
      }
      ## k1
      for(i in 1:(d.size*d.size)){
        assign(Ystate[i],Yt[i])
      }

      for(i in 1:(d.size*d.size)){
        k1[i] <- dt*eval(F[i])
      }
      ## k2
      for(i in 1:(d.size*d.size)){
        assign(Ystate[i],Yt[i]+k1[i]/2)
      }
      for(i in 1:(d.size*d.size)){
        k2[i] <- dt*eval(F[i])
      }
      ## k3
      for(i in 1:(d.size*d.size)){
        assign(Ystate[i],Yt[i]+k2[i]/2)
      }
      for(i in 1:(d.size*d.size)){
        k3[i] <- dt*eval(F[i])
      }
      ## k4

##C³
      ## Xt+1
      for(i in 1:d.size){
        assign(state[i],X.t0[t+1,i])   
      }
##

      for(i in 1:(d.size*d.size)){
        assign(Ystate[i],Yt[i]+k3[i])
      }
      for(i in 1:(d.size*d.size)){
        k4[i] <- dt*eval(F[i])
      }
      ## F(Y(t+dt))
      k <- (k1+k2+k2+k3+k3+k4)/6
      Yt <- Yt+k
      Y <- rbind(Y,Yt)
    }
    ## return matrix : (division+1)*(d.size*d.size)
    rownames(Y) <- NULL
    colnames(Y) <- Ystate
    return(ts(Y,deltat=dt,start=0))  
  }

  # function to calculate Y_{t}Y_{s}^{-1}
  Get.YtYis <- function(t,s,range){

##C³
#    yt <- matrix(Y[(range[t]-1)*delta/deltat(Y)+1,],d.size,d.size)
    yt <- matrix(Y[range[t],],d.size,d.size)
#    yis <- solve(matrix(Y[(range[s]-1)*delta/deltat(Y)+1,],d.size,d.size))
    yis <- solve(matrix(Y[range[s],],d.size,d.size))
##

    return(yt%*%yis)
  }
  
  # function to calculate lambda_{t,s}
  ## require: de.diffusion, ytyis
  Get.lambda.ts <- function(t,s,range){
    tmp <- matrix(0,d.size,r.size)
    assign(pars,0)  ## epsilon=0
    #ytyis <- Get.YtYis((range[t]-1)*delta,(range[s]-1)*delta)
    for(i in 1:d.size){

##C³
#      assign(state[i],X.t0[(range[s]-1)*delta/deltat(X.t0)+1,i])
      assign(state[i],X.t0[range[s],i])
##
    }
    for(i in 1:d.size){
      for(j in 1:r.size){
        tmp[i,j] <- eval(de.diffusion[[i]][j]) # dV/de
      }
    }
    return(ytyis[t,s,,]%*%tmp)
  }
  
  # function to calculate lambda_{t,s,0}
  ## require: de.drift, ytyis
  Get.lambda.ts0 <- function(t,s,range){
    tmp <- seq(0,0,length=d.size)
    assign(pars,0)  ## epsilon=0
    #ytyis <- Get.YtYis((range[t]-1)*delta,(range[s]-1)*delta)
    for(i in 1:d.size){

##C³
#      assign(state[i],X.t0[(range[s]-1)*delta/deltat(X.t0)+1,i])
      assign(state[i],X.t0[range[s],i])
##

    }
    for(i in 1:d.size){
      tmp[i] <- eval(de.drift[i]) # dV0/de
    }
    return(ytyis[t,s,,]%*%tmp)
  }

  # function to calculate mu_{i,t,s}
  ## require: ytyis
  Get.mu.its <- function(i.state,t,s,range){
    tmp <- matrix(0,d.size,r.size)
    assign(pars,0)  ## epsilon=0
    #ytyis <- Get.YtYis((range[t]-1)*delta,(range[s]-1)*delta)
    for(i in 1:d.size){

##C³
#      assign(state[i],X.t0[(range[s]-1)*delta/deltat(X.t0)+1,i])
      assign(state[i],X.t0[range[s],i])
##
    }
	# make expression of dV/di
    diV <- as.list(NULL)
    for(i in 1:d.size){
      diV[i] <- list(Derivation.scalar(V[[i]],state[i.state],r.size))
    }
	# make expression of d(dV/di)/de
    dideV <- as.list(NULL)
    for(i in 1:d.size){
      dideV[i] <- list(Derivation.scalar(diV[[i]],pars,r.size))
    }
	# evaluate expression
    for(i in 1:d.size){
      for(j in 1:r.size){
        tmp[i,j] <- eval(dideV[[i]][j])
      }
    }
    return(ytyis[t,s,,]%*%tmp)
  }
  
  # function to calculate mu_{i,t,s,0}
  ## require: ytyis
  Get.mu.its0 <- function(i.state,t,s,range){
    tmp <- seq(0,0,length=d.size)
    assign(pars,0)  ## epsilon=0
    #ytyis <- Get.YtYis((range[t]-1)*delta,(range[s]-1)*delta)
    for(i in 1:d.size){

##C³
#      assign(state[i],X.t0[(range[t]-1)*delta/deltat(X.t0)+1,i])
#      assign(state[i],X.t0[(range[s]-1)*delta/deltat(X.t0)+1,i])
      assign(state[i],X.t0[range[s],i])
##

    }
    diV0 <- Derivation.scalar(V0,state[i.state],d.size) # dV0/di

##C³
#    dideV0 <- Derivation.scalar(V0,pars,d.size) #d(dV0/di)/de
    dideV0 <- Derivation.scalar(diV0,pars,d.size) #d(dV0/di)/de
##
    for(i in 1:d.size){
      tmp[i] <- eval(dideV0[i])
    }
    return(ytyis[t,s,,]%*%tmp)
  }
  
  # function to calculate nu_{i,j,t,s}
  ## require: ytyis
  Get.nu.ijts <- function(i.state,j.state,t,s,range){
    tmp <- seq(0,0,length=d.size)
    assign(pars,0)  ## epsilon=0
    #ytyis <- Get.YtYis((range[t]-1)*delta,(range[s]-1)*delta)
    for(i in 1:d.size){

##C³
#      assign(state[i],X.t0[(range[s]-1)*delta/deltat(X.t0)+1,i])
      assign(state[i],X.t0[range[s],i])
##
    }
    diV0 <- Derivation.scalar(V0,state[i.state],d.size)     #dV0/di
    didjV0 <- Derivation.scalar(diV0,state[j.state],d.size) #d(dV0/di)/dj
    for(i in 1:d.size){
      tmp[i] <- eval(didjV0[i])
    }
    return(ytyis[t,s,,]%*%tmp)
  }
  
  # function to calculate nu_{t,s}
  ## require: dede.diffusion, ytyis
  Get.nu.ts <- function(t,s,range){
    tmp <- matrix(0,d.size,r.size)
    assign(pars,0)  ## epsilon=0
    #ytyis <- Get.YtYis((range[t]-1)*delta,(range[s]-1)*delta)
    for(i in 1:d.size){

##C³
#      assign(state[i],X.t0[(range[s]-1)*delta/deltat(X.t0)+1,i])
      assign(state[i],X.t0[range[s],i])
##
    }
    for(i in 1:d.size){
      for(j in 1:r.size){
        tmp[i,j] <- eval(dede.diffusion[[i]][j])
      }
    }
    return(ytyis[t,s,,]%*%tmp)
  }
  
  # function to calculate nu_{t,s,0}
  ## require: dede.drift, ytyis
  Get.nu.ts0 <- function(t,s,range){
    tmp <- seq(0,0,length=d.size)
    assign(pars,0)  ## epsilon=0
    #ytyis <- Get.YtYis((range[t]-1)*delta,(range[s]-1)*delta)
    for(i in 1:d.size){

##C³
#      assign(state[i],X.t0[(range[s]-1)*delta/deltat(X.t0)+1,i])
      assign(state[i],X.t0[range[s],i])
##
    }
    for(i in 1:d.size){
      tmp[i] <- eval(dede.drift[i])
    }
    return(ytyis[t,s,,]%*%tmp)
  }
  
  # function to calculate mu in thesis p5
  funcmu <- function(e=0){

    division <- nrow(X.t0)
    XT <- X.t0[division,] #data X0 observated last. size:vector[d.size]
    
    ## calculate derived F by e with XT and e=0. deF(XT,0)
    deF0 <- c() #size:vector[k.size]
    for(d in 1:d.size){
      assign(state[d],XT[d]) #input XT in state to use eval function
    }
    for(k in 1:k.size){
      deF <- deriv(F[k],"e") #expression of derived F by e
      deF0[k] <- attr(eval(deF),"gradient") #calculate deF(derived F by e) with XT
    }
    
    ## calculate derived f0 by e with Xt and e=0. def0(Xt,0)
    def0 <- matrix(0,k.size,division) #size:matrix[k.size,division]
    for(k in 1:k.size){
      def <- deriv(f[[1]][k],"e") #expression of derived f0 by e
      for(t in 1:division){
        X0t <- X.t0[t,]  #data X0 observated on time t
        for(d in 1:d.size){ 
          assign(state[d],X0t[d]) #input X0t in state
        }
        def0[k,t] <- attr(eval(def),"gradient")  #calculate def(derived f0 by e) with X0t
      }
    }
    
	# integrate def0 (just sum it)
    def0 <- apply(def0,1,sum) #sum of def0. size:vector[k.size]

##C³
#    def0 <- def0*(1/(division-1))
    def0 <- def0*(Terminal/(division-1))
##

    mu <- def0+deF0  #size:vector[k.size]

##C³
    dxF <- c()
    dxf <- c()

    for(k in 1:k.size){
      dxF[k] <- deriv(F[k],state)
      dxf[k] <- deriv(f[[1]][k],state)
    }

    tmp1 <- matrix(0,d.size,division-1)
    for(t in 1:(division-1)){
	for(d in 1:d.size){
	  assign(state[d],X.t0[t,])
	}
#	browser()
		deV0 <- numeric(d.size)
		for(i in 1:d.size) deV0[i]	<- eval(de.drift[i])

##C³
	tmp1[,t] <- deV0 %*% solve(matrix(Y[t,],d.size,d.size))
##
    }

    for(d in 1:d.size){
	assign(state[d],X.t0[division,d])
    }

    dxFT <- matrix(0,k.size,d.size)
    for(k in 1:k.size){
      dxFT[k,] <- attr(eval(dxF[k]),"gradient")
    }

    tmp2 <- as.matrix(rep(1,division-1))
    third <- dxFT %*% matrix(Y[division,],d.size,d.size) %*% (tmp1 %*% tmp2) * (Terminal/(division-1))

    dxf0 <- array(0,dim=c(k.size,d.size,division))
    for(t in 1:division){
	for(d in 1:d.size){
	  assign(state[d],X.t0[t,d])
	}
	for(k in 1:k.size){
	  dxf0[k,,t] <- attr(eval(dxf[k]),"gradient")
	}
    }

    tmp3 <- matrix(0,k.size,division-2)
    for(t in 2:(division-1)){
	tmp4 <- tmp1[,1:(t-1)]
	tmp5 <- as.matrix(rep(1,t-1))
##C³
	tmp3[,t-1] <- matrix(dxf0[,,t],k.size,d.size) %*% matrix(Y[t,],d.size,d.size) %*% (tmp4 %*% tmp5) * (Terminal/(division-1))
##
    }

    tmp6 <- as.matrix(rep(1,division-2))
    fourth <- tmp3 %*% tmp6 * (Terminal/(division-1))

    mu <- mu + third + fourth
##

    return(mu)
  }

  # function to calculate a_{s}^{alpha} in bookchapter p5
  funca <- function(e=0){ 
    #init
    division <- nrow(X.t0)
    XT <- X.t0[division,] #data X0 observated last. size:vector[d.size]
    defa <- array(0,dim=c(k.size,r.size,division)) #size:array[k.size,r.size,division]
    deva <- array(0,dim=c(d.size,r.size,division)) #size:array[d.size,r.size,division]
    dxF0 <- matrix(0,k.size,d.size) #size:matrix[k.size,d.size]
    dxf0 <- array(0,dim=c(k.size,d.size,division)) #size:array[k.size,d.size,division]
    dxf <- c()
    dxF <- c()
    def <- list()
    dev <- list()
	
	# prepare expression of derivatives
    for(k in 1:k.size){
      dxf[k] <- deriv(f[[1]][k],state) #expression of d f0/dx
      dxF[k] <- deriv(F[k],state) #expression of d F/dx
      def[[k]] <- list()
      for(r in 2:(r.size+1)){
        def[[k]][[r-1]] <- deriv(f[[r]][k],"e") #expression of derived fa by e
      }
    }
    for(r in 1:r.size){
      dev[[r]] <- list()
      for(d in 1:d.size){
        dev[[r]][[d]] <- deriv(V[[d]][r],"e") #expression of derived Vr by e
      }
    }

	# evaluate derivative expressions
    for(t in 1:division){
      X0t <- X.t0[t,]
      for(d in 1:d.size){
        assign(state[d],X0t[d]) #input X0t in state to use eval function to V
      }
      for(k in 1:k.size){
        ##calculate derived f0 by x with Xt and e=0. dxf0(Xt,0)
        dxf0[k,,t] <- attr(eval(dxf[k]),"gradient") #calculate dxf(derived f0 by x) with X0t
        ##calculate derived F by x with XT and e=0. dxF(XT,0)
        dxF0[k,] <- attr(eval(dxF[k]),"gradient") #calculate dxF(derived F by x) with XT
        for(r in 2:(r.size+1)){
          ##calculate derived fa by e with Xt and e=0. defa(Xt,0)
          defa[k,r-1,t] <- attr(eval(def[[k]][[r-1]]),"gradient")  #calculate def(derived fa by e) with X0t   
        }
      }
      for(r in 1:r.size){
        for(d in 1:d.size){
          ##calculate derived Va by e with Xt and e=0. deVa(Xt,0)
          deva[d,r,t] <- attr(eval(dev[[r]][[d]]),"gradient") #calculate dev(derived Vr by e) with X0t
        }
      }
    }  

	# prepare Y and Y^{-1}
    arrayY <- array(0,dim=c(d.size,d.size,division))
    invY <- array(0,dim=c(d.size,d.size,division))
    for(t in 1:division){
      arrayY[,,t] <- matrix(Y[t,],d.size,d.size)
      invY[,,t] <- solve(arrayY[,,t])
    }

	# calculate dxF*Y^{T}*Y^{-1}*deV_{a}
	second <- array(0,dim=c(k.size,r.size,division))
	temp <- dxF0 %*% arrayY[,,division]
	for(t in 1:division) {
		second[,,t] <- temp %*% invY[,,t] %*% deva[,,t]
	}

	#calculate integral
	fIntegral <- array(0,dim=c(k.size,r.size,division))
	third <- array(0,dim=c(k.size,r.size,division))

##C³
#	dt <- Terminal / division
	dt <- Terminal / (division-1)
#	third[,,division] <- dxf0[,,division] %*% arrayY[,,division] %*% invY[,,division] %*% deva[,,division] * dt
	third[,,division] <- 0
##

	for(s in (division-1):1){

##C³
#		third[,,s] <- third[,,s+1] + dxf0[,,s] %*% arrayY[,,s] %*% invY[,,s] %*% deva[,,s] * dt
		third[,,s] <- third[,,s+1] + matrix(dxf0[,,s],k.size,d.size) %*% arrayY[,,s] %*% invY[,,s] %*% matrix(deva[,,s],d.size,r.size) * dt
##
	}

#	for(s in 1:(division-1)){
#	  tmp <- array(0,dim=c(k.size,r.size,division-s))
#	  for(t in s:(division-1)){
#print(c(s,t))
#	    tmp[,,t-s+1] <- dxf0[,,t] %*% arrayY[,,t] %*% invY[,,s] %*% deva[,,s] * dt

#	    third[,,s] <- third[,,s] + tmp[,,t-s+1]
#	  }
#	}

 #   defa <- aperm(defa,c(1,3,2))*1.0
 #   deva <- aperm(deva,c(1,3,2))*1.0
 #   dxF0 <- dxF0*1.0
 #   dxf0 <- dxf0*1.0
    
    ##use C source
 #   dyn.load("yuima.so")
 #   a <- .Call("get_a",defa,dxF0,arrayY,invY,deva,dxf0,
 #              1.0*dim(defa),1.0*dim(arrayY),1.0*dim(invY),1.0*dim(deva),1.0*dim(dxf0),
 #              1.0*a,1.0*dim(a))

    return(defa + second + third) #size:array[k.size,r.size,division]
  }
  
  # function to calculate sigma in thesis p5
  # require: aMat
  funcsigma <- function(e=0){
    division <- nrow(X.t0)
    sigma <- matrix(0,k.size,k.size) #size:matrix[k.size,k.size]

##C³
#    for(t in 1:division){
    for(t in 1:(division-1)){

#      sigma <- sigma+(aMat[,,t]%*%t(aMat[,,t])) /(division-1) #calculate sigma
#      sigma <- sigma+(aMat[,,t]%*%t(aMat[,,t])) * Terminal/(division-1) #calculate sigma
	a.t <- matrix(aMat[,,t],k.size,r.size)
	sigma <- sigma+(a.t%*%t(a.t)) * Terminal/(division-1)
##
    }
    if(any(eigen(sigma)$value<=0.0001)){
    # Singularity check
      yuima.warn("Eigen value of covariance matrix in very small.")
    }
    return(sigma)
  }
  
  ## integrate start:1 end:t number to integrate:block
  # because integration at all 0:T takes too much time,
  # we sample only 'block' points and use trapezium rule to integrate  
  make.range.for.trapezium.fomula <- function(t,block){
    if(t/block <= 1){ # block >= t : just use all points 
      range <- c(1:t)
    }else{ # make array which includes points to use
      range <- as.integer( (c(0:block) * (t/block))+1)
      if( range[block+1]  < t){
        range[block+2] <- t
      }else if( range[block+1] > t){
        range[block+1] <- t
      }
    }
    return(range)
  }

  # function to return expressions of df0/dxi
  deriv.f0.for.state<- function(f0){
    tmp_deriv_f0 <- function(i,f){
      d_xi_f <- c()
      for(j in 1:k.size){
        d_xi_f[j] <- deriv(f[j],state[i])
      }
      return(d_xi_f)
    }
    tmp <- apply(as.matrix(1:length(state)),1,tmp_deriv_f0,f0)
    ##return list of (df0/dxi)
    return(tmp)
  }
  
  ## This function return value of expr(X0[t])
  ## expr:list of derived for state x1,x2,...
  ## ex) expr[[1]] : f0 derived for x1
  ## t: time index (t=1,2,...,division+1)
  input.deriv <- function(t,expr,l=1,X0){
    df <- c()
    ##input x1,x2,...,
    for(i in  1:d.size){
      df[i] <- expr[[i]][l]
      assign(state[i],X0[t,i])
    }
    ##epsilon = 0
    assign(pars[1],0)
    tmp <- c()
    for(i in 1:d.size){
      tmp[i] <- attr(eval(df[i]),"gradient")
    }
    return(tmp)
  }

  ## get hessian for(state,"e")
  hessian.f0.di.de<- function(f0){
    tmp.hessian.f0 <- function(i,f){
      d.xi.de.f <- c()
      for(j in 1:k.size){
        d.xi.de.f[j] <- deriv(f[j],c(state[i],"e"),hessian=T)
      }
      return(d.xi.de.f)
    }
    tmp <- apply(as.matrix(1:length(state)),1,tmp.hessian.f0,f0)
    ##return list of (d^2 f0/dxi de)
    return(tmp)
  }
  
  hessian.f.dxi.de<- function(f){
    list_l_dxi_de <- list(NULL)
    list_dxi_de <- list(NULL)
    de <- c()
    for(k in 1:k.size){
      for(i in 1:d.size){
        for(r in 1:r.size){
          de[r] <- deriv(f[[r+1]][k],c(state[i],"e"),hessian=T)
        }
        list_dxi_de[[i]] <- de
      }
      list_l_dxi_de[[k]] <- list_dxi_de
    }
    ##return list of (d^2 f0/dxi de)
    return(list_l_dxi_de)
  }

  ## This function return list deriv f0 for (xi,xj)
  ## list_k_dxi_dxj[[ k ]][[ i ]][ j ] is expression f0[k] derived for (state[i],state[j])
  ## f0:1,...,k.size expression
  hessian.f0.di.dj<- function(f0){
    list_k_dxi_dxj <- as.list(NULL)
    list_dxi_dxj<- as.list(NULL)
    dxj <- c()
    for(k in 1:k.size){
      for(i in 1:d.size){
        for(j in 1:d.size){
          dxj[j] <- deriv(f0[k],c(state[i],state[j]),hessian=T)
        }
        list_dxi_dxj[[i]] <- dxj
      }
      list_k_dxi_dxj[[k]] <- list_dxi_dxj
    }
    return(list_k_dxi_dxj)
  }

  # following funcs (named 'input.hessian~') solve expression of hessian
  ## This function return (d x 1 vector)
  input.hessian <- function(t,expr,n=1,m=2,l=1,X0){
    h_f <- c()
    ##input x1,x2,...,
    for(i in  1:d.size){
      h_f[i] <- expr[[i]][l]
      assign(state[i],X0[t,i])
    }
    ##epsilon = 0
    assign(pars[1],0)
    tmp <- c()
    for(i in 1:d.size){
      tmp[i] <- attr(eval(h_f[i]),"hessian")[1,n,m]
    }
    return(tmp)
  }

  ## This function returns (d x d matrix)
  ## t: time index expr:derived expression 
  input.hessian.dxi.dxj <- function(t,expr,n=1,m=2,l=1,X0){
    h_f <- list(NULL)
    ##input x1,x2,...,
    for(i in  1:d.size){
      h_f[[i]] <- expr[[l]][[i]]
      assign(state[i],X0[t,i])
    }
    ##epsilon = 0
    assign(pars[1],0)
    
    tmp <- matrix(0,d.size,d.size)
    for(i in 1:d.size){
      for(j in 1:d.size){
        if(i == j){
          tmp[i,j] <- attr(eval(h_f[[i]][j]),"hessian")[1,1,1]
        }else{
          tmp[i,j] <- attr(eval(h_f[[i]][j]),"hessian")[1,i,j]
        }
      }
    }
    dim(tmp) <- c()
    return(tmp)
  }
  
  ## This function returns (d x r matrix)
  ## t:index , expr:hessian for f, l=1,...,k.size
  input.hessian.dxi.de.f<- function(t,expr,n=1,m=2,l=1,X0){
    h.f <- list(NULL)
    h.f.r <- c()
    for(i in 1:d.size){
      for(r in 1:r.size){
        h.f.r[r] <- expr[[l]][[i]][r]
      }
      assign(state[i],X0[t,i])
      h.f[[i]] <- h.f.r
    }
    assign(pars[1],0)
    tmp <- matrix(0,d.size,r.size)
    for(d in 1:d.size){
      for(r in 1:r.size){
        tmp[d,r] <- attr(eval(h.f[[d]][r]),"hessian")[1,n,m]
      }
    }
    return(tmp)
  }
  
  ## multi dimension gausian distribusion
  normal <- function(x,mu,Sigma){
    if(length(x)!=length(mu)){
      print("Error:length of x != one of mu")
    }
    dimension <- length(x)
    tmp <- 1/((sqrt(2*pi))^dimension * sqrt(det(Sigma))) * exp(-1/2 * t((x-mu)) %*% solve(Sigma) %*% (x-mu) )
    return(tmp)
  }

  ## get d0
  ## required library(adapt)
  get.d0.term <- function(){
    lambda <- eigen(Sigma)$values
    matA <- eigen(Sigma)$vector
    ## get g(z)*pi0(z)

##C³
#    gz_pi0 <- function(z){
#      return( G(z) * H0 *normal(z,mu=mu,Sigma=Sigma))      
#    }

    gz_pi0 <- function(z){
	tmpz <- matA %*% z
      return( G(tmpz) * H0 *normal(tmpz,mu=mu,Sigma=Sigma))   #det(matA) = 1
    }
##

    gz_pi02 <- function(z){
      return( G(z) * H0 *dnorm(z,mean=mu,sd=sqrt(Sigma)))
    }

    ## integrate
    if( k.size ==1){ # use 'integrate' if k=1

##’ˆÓF’†S‚©‚ç‚¸‚ê‚·‚¬‚é‚ÆÏ•ª‚ğŒvZ‚µ‚È‚¢	#”ÍˆÍ‚ª•Ï‚í‚é‚Æ’l‚ª•Ï‚í‚Á‚Ä‚µ‚Ü‚¤
      tmp <- integrate(gz_pi02,mu-5*sqrt(Sigma),mu+5*sqrt(Sigma))$value
    }else if( 2 <= k.size || k.size <= 20 ){ # use library 'adapt' to solve multi-dimentional integration

##C³
#	  max <- 6 * sqrt(lambda.max)
#      min <- -6 * sqrt(lambda.max)
#      L <- (max - min)

	lower <- c()
	upper <- c()

	for(k in 1:k.size){
	  max <- 7 * sqrt(lambda[k])
	  lower[k] <- - max
	  upper[k] <- max
	}
##

#	if(require(adapt)){
	  if(require(cubature)){
#tmp <- adapt(ndim=k.size,lower=rep(min,k.size),upper=rep(max,k.size),functn=gz_pi0)$value

##C³
#        tmp <- adaptIntegrate(gz_pi0, lower=rep(min,k.size),upper=rep(max,k.size))$integral
	tmp <- adaptIntegrate(gz_pi0, lower=lower,upper=upper)$integral
##
		} else {
	   tmp <- NA		
	  }
    }else{
      stop("length k is too big.")
    }
    return(tmp)
  }
  
  ###############################################################################
  # following funcs are part of d1 term
  # because they are finally integrated at 'get.d1.term()',
  # these funcs are called over and over again.
  # so, we use trapezium rule for integration to save time and memory.
  
  # these funcs almost calculate each formulas including trapezium integration.
  # see each formulas in thesis to know what these funcs do.
  
  # some funcs do alternative calculation at k=1.
  # it depends on 'integrate()' function
  ###############################################################################
  
  # p.9 Lemma2 (a)
  
  ## This function returns First term of di.bar (d.size x block martix) and
  ## part of second term(Second.tmp)
  ## Second.tmp is ((d x block) x k) matrix
  ## aMat.tmp is k x (r x block) matrix

  get.di.bar.init <- function(){
    # trapezium rule
    tmp.mat <-  rep(1,block)
    tmp.mat2 <- rep(Diff,d.size)
    dim(tmp.mat2) <- c((block*block),d.size)
    tmp.mat2 <- t(tmp.mat2)
    dim(tmp.mat2) <- c((d.size*block),block)
    First <- (lambda.ts0 * tmp.mat2) %*% tmp.mat * delta
    dim(First) <- c(d.size,block)
	
    tmp.mat3 <- tmp.mat2
    for(i in 1:r.size){
      if(i != 1){

##C³?
#        tmp.mat3 <- rbind(tmp.mat3,tmp.mat2)
        tmp.mat3 <- cbind(tmp.mat3,tmp.mat2)
##

      }
    }
    dim(tmp.mat3) <- c(d.size*block,r.size*block)
    Second.tmp <- (lambda.ts * tmp.mat3) %*% t(aMat.tmp) * delta
    tmp <- list(First=First,
                Second.tmp=Second.tmp)
    return(tmp)
  }

  ## dependency:dat.di.bar
  get.di.bar <- function(x){
    if(k.size ==1){
      First <- dat.di.bar$First
      
      tmp.di.bar <- dat.di.bar$Second.tmp
      dim(tmp.di.bar) <- c(d.size,block)
      Second.tmp <- tmp.di.bar
      for(i in 1:(length(x))){
        if(i!=1){
          First <- rbind(First,dat.di.bar$First)  
        }
      }
      for(i in 1:length(x)){
        if(i!=1){
          Second.tmp <- rbind(Second.tmp,tmp.di.bar)
        }
      }
      tmp.x <- x
      for(i in 1:d.size){
        if(i !=1){
          tmp.x <- rbind(tmp.x,x)
        }
      }
      dim(tmp.x) <- c()
      tmp.x <- rep(tmp.x,block)
      dim(tmp.x) <- c((d.size*length(x)),block)
      Second <-  tmp.x * Second.tmp / as.double(Sigma)
      tmp <- First + Second
    }else{
      Second <- dat.di.bar$Second.tmp %*% solve(Sigma) %*% x
      dim(Second) <- c(d.size,block)
      tmp <- dat.di.bar$First  + Second
    }
    return(tmp)
  }

  ## h is (d x block matrix)
  ## x is k dimension vector
  get.Di.bar <- function(h,x){
    if(k.size==1){
      tmp.Diff <- Diff[block,]
      for(i in 1:d.size){
        if( i != 1 ){
          tmp.Diff <- rbind( tmp.Diff,Diff[block,] )
        }
      }
      tmp.h <- h * tmp.Diff
      for(i in 1:length(x)){
        if(i !=1){
          tmp.h <- rbind( tmp.h , (h*tmp.Diff) )
        }
      }
      tmp <- tmp.h * get.di.bar(x) *delta
      tmp <- tmp %*% rep(1,block)
      dim(tmp) <- c(d.size,length(x))
    }else{
      tmp.Diff <- Diff[block,]
      for(i in 1:d.size){
        if( i != 1 ){
          tmp.Diff <- rbind( tmp.Diff,Diff[block,] )
        }
      }
      tmp <- h * get.di.bar(x) * tmp.Diff * delta
      tmp <- as.vector(tmp %*% rep(1,block))
    }
    return(tmp)
  }
  
  # p.10 (b)
  
  ## this function returns first term and part of second term of Di
  ## h: d x (r x block)  matrix
  ## x: k dimension vector
  ## dependency:dat.di.bar
  get.Di.init <- function(h){
    ## First term
    tmp1 <- dat.di.bar$First
    for(i in 1:r.size){
      if(i !=1){
        tmp1 <- rbind( tmp1 , dat.di.bar$First )
      }
    }
    dim(tmp1) <- c( d.size , r.size * block)
    tmp.Diff <- Diff[block,]
    for(i in 1:(k.size * r.size)){
      if(i != 1) tmp.Diff <- rbind(tmp.Diff,Diff[block,])
    }
    dim(tmp.Diff) <- c(k.size,(r.size*block))
    First <- (tmp1 * h) %*% t(aMat.tmp * tmp.Diff) %*% solve(Sigma) * delta
    ## End of First term
    
    tmp.mat1 <- t(dat.di.bar$Second.tmp)
    dim(tmp.mat1) <- c((k.size * d.size),block)
    Second.tmp <- array(0,dim=c(k.size,k.size,d.size))
    for( i in 1:d.size){
      tmp.mat2 <- tmp.mat1[((i-1)*k.size + 1):(i*k.size),]
      dim(tmp.mat2) <- c()
      tmp.mat3 <- tmp.mat2
      tmp2 <-  h[i,]
      dim(tmp2) <- c(r.size,block)
      tmp.mat4 <- tmp2
      for(j in 1:r.size){
        if(j != 1) tmp.mat3 <- rbind( tmp.mat3 , tmp.mat2 )
      }
      for(j in 1:k.size){
        if(j !=1) tmp.mat4 <- rbind( tmp.mat4 , tmp2 )
      }
      dim(tmp.mat4) <- c(r.size,(k.size * block))
      tmp3 <- (tmp.mat4 * tmp.mat3)
      tmp4 <- matrix(0,(r.size*block),k.size)
      for(k in 1:k.size){
        tmp4[,k] <- tmp3[,(c(1:block)*k.size + ( k - k.size))]
      }
      tmp.Diff <- Diff[block,]

##’ˆÓF‘äŒ`Œö®‚ª‚æ‚­‚È‚¢‚ª‚ ‚éB

      for( j in 1:(k.size*r.size)){
        if(j !=1){
          tmp.Diff <- rbind(tmp.Diff,Diff[block,])
        }
      }
      dim(tmp.Diff)<-c(k.size,(r.size*block))
      Second.tmp[,,i] <- (aMat.tmp * tmp.Diff) %*% tmp4 * delta
    }
    tmp <- list(First=First,
                Second.tmp=Second.tmp
                )
  }

  ## dependency: dat.Di.init include First and Second.tmp
  get.Di <- function(dat.Di.init,x){
    if(k.size==1){
      tmp <- numeric(d.size)
      for(i in 1:d.size){
        tmp[i] <- dat.Di.init$Second.tmp[,,i]
      }
      ##dat.Di.init$First is d x k matrix

##C³
#      First <- dat.Di.init$First %*% x / as.double(Sigma)
      First <- dat.Di.init$First %*% x 	#get.Di.init‚Åsolve(Sigma)‚ğ‚©‚¯‚Ä‚¢‚é
      Second <- (tmp/(as.double(Sigma)^2) ) %*% t(x^2 - as.double(Sigma))
    }else{
#      First <- dat.Di.init$First %*% solve(Sigma) %*% x
      First <- dat.Di.init$First %*% x	#get.Di.init‚Åsolve(Sigma)‚ğ‚©‚¯‚Ä‚¢‚é
##
      Second <- numeric(d.size)
      for(i in 1:d.size){
        tmp <- dat.Di.init$Second.tmp[,,i]
        Second[i] <-  sum( diag( solve(Sigma) %*% tmp %*% solve(Sigma) %*% (x%*%t(x) - Sigma) ) )
      }
    }
    Di <- First + Second
    return(Di)
  }

  # p.10 (c)
  
  ## preparation to calculate d^{ij}(x)_{t}
  ## dependency: dat.di.bar
  get.dij.init <- function(){
    ## First term
    tmp1.1<- dat.di.bar$First
    tmp1.2 <- rep(dat.di.bar$First,d.size)
    dim(tmp1.2) <- c((d.size*block),d.size)
    tmp1.2 <- t(tmp1.2)
    dim(tmp1.2) <- c((d.size^2),block)
    for(i in 1:d.size){
      if(i !=1){
        tmp1.1 <- rbind( tmp1.1 , dat.di.bar$First )
      }
    }
    First <- tmp1.1 * tmp1.2 ## dim(First) = c( d.size^2 , block)
    
    ## Second term
    tmp2.1 <- tmp1.1
    dim(tmp2.1) <- c()
    tmp2.1 <- rep(tmp2.1,k.size)
    dim(tmp2.1) <- c((d.size^2)*block,k.size)
    tmp2.2 <- rep(dat.di.bar$Second.tmp,d.size)
    dim(tmp2.2) <- c(d.size*k.size*block,d.size)
    tmp2.2 <- t(tmp2.2)
    dim(tmp2.2) <- c(d.size^2*block,k.size)
    Second.tmp <- (tmp2.1 * tmp2.2) %*% solve(Sigma) ## dim(Second.tmp) = c(d.size^2 * block ,k.size)

    ## Third term
    tmp3.1 <- rep(dat.di.bar$First,d.size)
    dim(tmp3.1) <- c(d.size*block,d.size)
    tmp3.1 <- as.vector(t(tmp3.1))
    tmp3.1 <- rep(tmp3.1,k.size)
    dim(tmp3.1) <- c((d.size^2)*block,k.size)
    
    tmp3 <- dat.di.bar$Second.tmp
    dim(tmp3) <- c(d.size,k.size*block)
    tmp3.2 <- tmp3
    for( i in 1:d.size){
      if(i != 1){
        tmp3.2 <- rbind(tmp3.2,tmp3)
      }
    }
    dim(tmp3.2) <- c((d.size^2)*block,k.size)
    Third.tmp <- (tmp3.1 * tmp3.2) %*% solve(Sigma) ## dim(Third.t,p) = c(d.size^2 * block ,k.size)
    
    ##Fourth term
    tmp4 <- t(dat.di.bar$Second.tmp)
    dim(tmp4) <- c(k.size*d.size,block)
    tmp4.1 <- tmp4
    for(i in 1:d.size){
      if(i != 1) tmp4.1 <- rbind(tmp4.1,tmp4)
    }
    # rm(tmp4)
    dim(tmp4.1) <- c(k.size,(d.size^2)*block)
    tmp4.2 <- tmp4.1
    for(k in 1:k.size){
      if(k != 1) tmp4.2 <- rbind(tmp4.2,tmp4.1)
    }
    # rm(tmp4.1)
    dim(tmp4.2) <- c(k.size,k.size*(d.size^2)*block)
    
    tmp4.3 <- t(dat.di.bar$Second.tmp)
    for(i in 1:d.size){
      if(i != 1) tmp4.3 <- rbind(tmp4.3,t(dat.di.bar$Second.tmp))
    }
    dim(tmp4.3) <- c()
    tmp4.4 <- tmp4.3
    for( k in 1:k.size){
      if( k != 1) tmp4.4 <- rbind(tmp4.4,tmp4.3)
    }
    # rm(tmp4.3)
    dim(tmp4.4) <- c(k.size,k.size*(d.size^2)*block)
    tmp4.5 <- tmp4.2 * tmp4.4
    tmp4.6 <- matrix(0,k.size*(d.size^2)*block,k.size)
    for(k in 1:k.size){
      tmp4.6[,k] <- tmp4.5[,(1:(d.size^2 * block)-1)*k.size + k ]
    }
    tmp4.6 <- t(tmp4.6)
    tmp4.7 <- solve(Sigma) %*%  (tmp4.5 + tmp4.6)
    # rm(tmp4.5)
    # rm(tmp4.6)
    tmp4.8 <- matrix(0,k.size*(d.size^2)*block,k.size)
    for(k in 1:k.size){
      tmp4.8[,k] <- tmp4.7[,(1:(d.size^2 * block)-1)*k.size + k ]
    }
    Fourth.tmp <- tmp4.8 %*% solve(Sigma) ## dim(Fourth.tmp) =c( k.size * d.size^2 *  block , k.size)
    
    ## Fifth term
    tmp.Diff <- rep(Diff,d.size)
    dim(tmp.Diff) <- c(block^2,d.size)
    tmp.Diff <- t(tmp.Diff)
    dim(tmp.Diff) <- c(d.size*block,block)
    tmp.Diff2 <- tmp.Diff
    for(r in 1:r.size){
      if(r != 1) tmp.Diff <- rbind(tmp.Diff,tmp.Diff2)
    }
    dim(tmp.Diff) <- c(d.size*block,r.size*block)
    # rm(tmp.Diff2)  
    Fifth <- matrix(0,d.size^2,block)
    
    for(t in 1:block){
      start <- (t-1) * d.size + 1
      end <- (t-1) * d.size +d.size
      if( d.size == 1 ){
        Fifth[,t] <- (lambda.ts[start:end,] * tmp.Diff[start:end,]) %*% lambda.ts[start:end,] *delta
      }else{
        Fifth[,t] <- (lambda.ts[start:end,] * tmp.Diff[start:end,]) %*% t(lambda.ts[start:end,]) *delta
      }
    }
    ## dim(Fifth) = c( d.size^2 ,block )
    tmp <- list(First=First,
                Second.tmp=Second.tmp,
                Third.tmp=Third.tmp,
                Fourth.tmp=Fourth.tmp,
                Fifth=Fifth
                )
    return(tmp)    
  }

  ## dependency: dat.dij
  get.dij <- function(x){
    if( k.size == 1 ){
      First <- dat.dij$First
      First.tmp <- First
      
      Fifth <- dat.dij$Fifth
      Fifth.tmp <- Fifth
      
      Second <- dat.dij$Second.tmp
      dim(Second) <- c(d.size^2,block)
      Second.tmp <- Second
      for( i in 1:length(x)){
        if(i != 1){
          First <- rbind(First,First.tmp)
          Second <- rbind(Second,Second.tmp)
          Fifth <- rbind(Fifth,Fifth.tmp)
        }
      }
      # rm(First.tmp)
      # rm(Second.tmp)
      # rm(Fifth.tmp)
      
      tmp.x <- rep(x,d.size^2 * block)
      dim(tmp.x) <- c(length(x),d.size^2 * block)
      tmp.x <- t(tmp.x)
      dim(tmp.x) <- c(block,d.size^2*length(x))
      tmp.x <- t(tmp.x)
      Second <- Second * tmp.x
      
      Third <- dat.dij$Third.tmp
      dim(Third) <- c(d.size^2,block)
      Third.tmp <- Third
      for( i in 1:length(x)){
        if(i != 1) Third <- rbind(Third,Third.tmp)
      }
      # rm(Third.tmp)
      Third <- Third * tmp.x
      
      Fourth <- dat.dij$Fourth.tmp
      dim(Fourth) <- c(d.size^2,block)
      Fourth.tmp <- Fourth
      for( i in 1:length(x)){
        if(i != 1) Fourth <- rbind(Fourth,Fourth.tmp)
      }
      # rm(Fourth.tmp)
      tmp.x <- rep((x^2 - Sigma),d.size^2 * block)
      dim(tmp.x) <- c(length(x),d.size^2 * block)
      tmp.x <- t(tmp.x)
      dim(tmp.x) <- c(block,d.size^2*length(x))
      tmp.x <- t(tmp.x)
      Fourth <- Fourth * tmp.x /2
      
      tmp <- First + Second + Third + Fourth + Fifth
    }else{
      Second <- dat.dij$Second.tmp %*% x
      dim(Second) <- c(d.size^2,block)
      
      Third <- dat.dij$Third.tmp %*% x
      dim(Third) <- c(d.size^2,block)
      
      tmp.x <- rep((x %*% t(x)) - Sigma, d.size^2 * block )
      dim(tmp.x) <- c(k.size , k.size * d.size^2 * block )
      
      Fourth <- rep(1,k.size) %*% t( dat.dij$Fourth.tmp * t(tmp.x) )
      dim(Fourth) <- c(k.size ,d.size^2 * block) 
      Fourth <- ( rep(1,k.size) %*% Fourth ) /2
      dim(Fourth) <- c(d.size^2,block)
      tmp <- dat.dij$First + Second + Third + Fourth + dat.dij$Fifth
    }    
    return( tmp )
  }

  get.Dij <- function(h,x){
    if(k.size == 1){
      tmp.Diff <- Diff[block,]
      for(i in 1:(d.size^2)){
        if( i != 1 ){
          tmp.Diff <- rbind( tmp.Diff,Diff[block,] )
        }
      }
      tmp.h <- h * tmp.Diff
      for(i in 1:length(x)){
        if(i !=1){
          tmp.h <- rbind( tmp.h , (h*tmp.Diff) )
        }
      }
      tmp <- tmp.h * get.dij(x) *delta
      tmp <- tmp %*% rep(1,block)
      dim(tmp) <- c(d.size^2,length(x))
    }else{
      tmp.Diff <- Diff[block,]
      for(i in 1:(d.size^2)){
        if( i != 1 ){
          tmp.Diff <- rbind( tmp.Diff,Diff[block,] )
        }
      }
      tmp <- h * get.dij(x) * tmp.Diff * delta
      tmp <- as.vector(tmp %*% rep(1,block))
    }
    return(tmp)
  }
  
  # p.11 (d)
  
  get.el.init <- function(){
    tmp.mat <-  rep(1,block)
    tmp.mat2 <- rep(Diff,d.size)
    dim(tmp.mat2) <- c((block*block),d.size)
    tmp.mat2 <- t(tmp.mat2)
    dim(tmp.mat2) <- c((d.size*block),block)
    First <- (nu.ts0 * tmp.mat2) %*% tmp.mat * delta
    
    dim(First) <- c(d.size,block)
    tmp.mat3 <- tmp.mat2
    for(i in 1:r.size){
      if(i != 1){
        tmp.mat3 <- rbind(tmp.mat3,tmp.mat2)
      }
    }
    dim(tmp.mat3) <- c(d.size*block,r.size*block)
    Second.tmp <- (nu.ts * tmp.mat3) %*% t(aMat.tmp) * delta
    tmp <- list(Third=First,
                Fifth.tmp=Second.tmp)
    return(tmp)
  }

  
  ## x : k.size dimension vector
  get.el <- function(x){
    if(k.size==1){
      First <- matrix(0,d.size*length(x),block)
      Second <- matrix(0,d.size*length(x),block)
      Third <- matrix(0,d.size*length(x),block)
      Fourth <- matrix(0,d.size*length(x),block)
      Fifth <- matrix(0,d.size*length(x),block)
      mu.its.l <- matrix(0,d.size,(block*r.size))
      mu.its0.l <- matrix(0,d.size,block)
      nu.ijts.l <- matrix(0,d.size^2,block)
      for(l in 1:d.size){
        for(t in 1:block){
          nu.ijts.l[1:(d.size)^2,] <- nu.ijts[,,(t-1)*d.size + l,]
          mu.its.l[1:d.size,] <- mu.its[,( (t-1)*d.size + l ),]
          mu.its0.l[1:d.size,] <- mu.its0[,( (t-1)*d.size + l ),]        
          
          idx <-((1:length(x))-1) * d.size + l
          ## First term
          First[idx,t] <- rep(1,d.size^2) %*% get.Dij(nu.ijts.l,x)
          ## Second term
          Second[idx,t] <- 2 * rep(1,d.size) %*%  get.Di.bar(mu.its0.l,x) 
          
          ## Fourth term
          Di.init <- get.Di.init( mu.its.l )
          Fourth[idx,t] <- 2 * rep(1,d.size) %*% get.Di( Di.init ,x)
        }
      }
      ## Third term
      Third <- dat.el$Third
      for(i in 1:length(x) ){
        if(i != 1) Third <- rbind(Third,dat.el$Third)
      }
      ## Fifth term
      tmp.el <- dat.el$Fifth.tmp
      dim(tmp.el) <- c(d.size,block)
      Fifth.tmp <- tmp.el
      for(i in 1:length(x)){
        if(i!=1) Fifth.tmp <- rbind(Fifth.tmp,tmp.el)
      }
      tmp.x <- x
      for(i in 1:d.size){
        if(i !=1) tmp.x <- rbind(tmp.x,x)
      }
      dim(tmp.x) <- c()
      tmp.x <- rep(tmp.x,block)
      dim(tmp.x) <- c((d.size*length(x)),block)
      
      Fifth <- tmp.x * Fifth.tmp / as.double(Sigma)

##C³      
#      tmp <- First + Second * Third + Fourth + Fifth
      tmp <- First + Second + Third + Fourth + Fifth
##

    }else{
      First <- matrix(0,d.size,block)
      Second <- matrix(0,d.size,block)
      Third <- matrix(0,d.size,block)
      Fourth <- matrix(0,d.size,block)
      Fifth <- matrix(0,d.size,block)
      mu.its.l <- matrix(0,d.size,(block*r.size))
      mu.its0.l <- matrix(0,d.size,block)
      nu.ijts.l <- matrix(0,d.size^2,block)
      for(l in 1:d.size){
        for(t in 1:block){
          nu.ijts.l[1:(d.size)^2,] <- nu.ijts[,,(t-1)*d.size + l,]
          mu.its.l[1:d.size,] <- mu.its[,( (t-1)*d.size + l ),]
          mu.its0.l[1:d.size,] <- mu.its0[,( (t-1)*d.size + l ),]        
          
          First[l,t] <- sum( get.Dij(nu.ijts.l,x) )        
          Second[l,t] <- 2 * sum( get.Di.bar(mu.its0.l,x) )
          
          Di.init <- get.Di.init( mu.its.l )
          Fourth[l,t] <- 2 * sum( get.Di( Di.init ,x) )
        }
      }
      Third <- dat.el$Third
      
      Fifth <- dat.el$Fifth.tmp %*% solve(Sigma) %*% x
      dim(Fifth) <- c(d.size,block)
      tmp <- First + Second + Third + Fourth + Fifth
    }
    
    return(tmp)
  }
  
  # function to calculate E^{l}(h;x)_{T}(or 'El') is included at 'get.Pl1()' as follows.
  # because El is used at both term 1 and 7 of Pl1,
  # El is better calculated IN 'get.Pl1()'
  
	
  ## This function returns Pl1(p.12)
  ## l=1,...,d.size
  get.Pl1 <- function(l,z){
  
  
    ## h: ( d x block ) matrix
    ## x: k dimension vector
    ## "get.Pl1.get.el" needed!
    get.El <- function(h,x){
      if(k.size==1){
        tmp.Diff <- Diff[block,]
        for(i in 1:d.size){
          if( i != 1 ){
            tmp.Diff <- rbind( tmp.Diff,Diff[block,] )
          }
        }
        tmp.h <- h * tmp.Diff
        for(i in 1:length(x)){
          if(i !=1){
            tmp.h <- rbind( tmp.h , (h*tmp.Diff) )
          }
        }
        tmp <- tmp.h * get.Pl1.get.el *delta
        tmp <- tmp %*% rep(1,block)
        dim(tmp) <- c(d.size,length(x))
      }else{
        tmp.Diff <- Diff[block,]
        for(i in 1:d.size){
          if( i != 1 ){
            tmp.Diff <- rbind( tmp.Diff,Diff[block,] )
          }
        }
        tmp <- h * get.Pl1.get.el * tmp.Diff * delta
        tmp <- as.vector(tmp %*% rep(1,block))
      }
      return(tmp)
    }
	
    #prepare get.el for First and Seventh
	get.Pl1.get.el <- get.el(z)
  
    if(k.size == 1){
      First <- (t(rep(1,d.size)) %*% get.El(di.f0.l[l,,],z) ) /2
      Second <- (t(rep(1,d.size^2)) %*% get.Dij(di.dj.f0.l[l,,],z)  ) /2
      Third <- (t(rep(1,d.size)) %*% get.Di.bar(di.de.f0.l[l,,],z) )
      Fourth <- sum( de.de.f0.l[l,] * Diff[block,] *delta) /2
      Fourth <- rep(Fourth,length(z))

		if(d.size==1){ 
			Di.init <- get.Di.init( t(as.matrix(di.de.f.l[l,,] )))
		}else{		
			Di.init <- get.Di.init( as.matrix(di.de.f.l[l,,]))
		}      
		Fifth <- (t(rep(1,d.size)) %*% get.Di( Di.init , z) )

      Sixth <- ( (de.de.f.l[l,] * Diff[block,]) %*% t(aMat.tmp) *delta) / as.double(Sigma) * z /2
      
      tmp.el <- get.Pl1.get.el[,block]
      dim(tmp.el) <- c(d.size,length(z))
      Seventh <- t(di.F.l[l,]) %*% tmp.el /2
      
      tmp.dij <- get.dij(z)[,block]
      dim(tmp.dij) <- c(d.size^2,length(z))
      Eighth <- di.dj.F.l[l,] %*% tmp.dij /2
      
      tmp.di.bar <-  get.di.bar(z)[,block]
      dim(tmp.di.bar) <- c(d.size,length(z))
      Ninth <- di.de.F.l[l,] %*% tmp.di.bar /2
      
      Tenth <- rep(de.de.F.l[l] /2 , length(z))
    }else{
      First <- sum( get.El(di.f0.l[l,,],z) ) /2
      Second <- sum( get.Dij(di.dj.f0.l[l,,],z) ) /2
      Third <- sum( get.Di.bar(di.de.f0.l[l,,],z) )
      Fourth <- sum( de.de.f0.l[l,] * Diff[block,] *delta) /2

      Di.init <- get.Di.init( di.de.f.l[l,,] )
      Fifth <- sum( get.Di( Di.init , z) )
      Sixth <- ((de.de.f.l[l,] * Diff[block,]) %*% t(aMat.tmp) *delta) %*% solve(Sigma) %*% z /2
      Seventh <-sum( di.F.l[l,] * get.el(z)[,block] )/2
      Eighth <- sum(di.dj.F.l[l,] * get.dij(z)[,block]) /2
      Ninth <- sum( di.de.F.l[l,] * get.di.bar(z)[,block] ) /2
      Tenth <- de.de.F.l[l] /2
    }
    Pl1 <- First + Second + Third + Fourth + Fifth + Sixth + Seventh + Eighth + Ninth + Tenth
    return(Pl1)
  }

  ## This function returns Pl2
  get.P2 <- function(z){
    if(k.size==1){
      First <- t(rep(1,d.size)) %*% get.Di.bar(di.rho,z)
      Second <- rep(de.rho %*% Diff[block,] *delta ,length(z))
    }else{
      First <- sum(get.Di.bar(di.rho,z))
      Second <- de.rho %*% Diff[block,] *delta
    }
    tmp <- First + Second
  }
  
  # now calculate pi1 using funcs above
  
  get.pi1 <- function(z){
    First <- 0
    z.tilda <- z - mu
    if(k.size ==1){
      First <- (get.Pl1(1,(z.tilda+delta)) * dnorm((z.tilda+delta),0,sqrt(Sigma)) -
                get.Pl1(1,z.tilda) * dnorm(z.tilda,0,sqrt(Sigma)))/delta
      Second <- get.P2(z.tilda) * dnorm(z.tilda,0,sqrt(Sigma))
    }else{
      tmp <- numeric(k.size)
      for(k in 1:k.size){
        dif <- numeric(k.size)
        dif[k] <- dif[k] + delta
        z.delta <- z.tilda + dif
        tmp[k] <- (get.Pl1(k,(z.delta)) * normal((z.delta),numeric(k.size),Sigma) -
                  get.Pl1(k,z.tilda) * normal(z.tilda,numeric(k.size),Sigma))/delta
      }
      First <- sum(tmp)
      Second <- get.P2(z.tilda) * normal(z.tilda,numeric(k.size),Sigma)
    }
    pi1 <- -H0*( First + Second )
    return(pi1)
  }

  # calculate d1
  ## required library(adapt)
  get.d1.term<- function(){

##C³
#    lambda.max <- max(eigen(Sigma)$values)
	lambda <- eigen(Sigma)$values
##

    ## get g(z)*pi1(z)
    
    gz_pi1 <- function(z){
      tmp <- G(z) * get.pi1(z)
      return( tmp  )
    }
    ## integrate
    if( k.size ==1){ # use 'integrate()'
      tmp <- integrate(gz_pi1,mu-5*sqrt(Sigma),mu+5*sqrt(Sigma))$value
    }else if(2 <= k.size || k.size <= 20 ){ # use sampling approximation for solving multi-dim integration. 2010/11/13
      set.seed(123)

##C³
#      max <- 10*lambda.max
#      min <- -10*lambda.max
#      tmp.x <- seq(min,max,length=100)
#      my.x <- NULL
#      for(k in 1:k.size){
#        my.x <- rbind(my.x,tmp.x)
#      }

#      est.ind <- seq(1,ncol(my.x),by=10)
#      est.points <- my.x[,est.ind]
#      width <- diff(tmp.x)
      
#      tmp <- 0
#      for(i in 1:(ncol(est.points)-1)){
#        tmp <- tmp + gz_pi1(est.points[,i])*width[i]
#      }


#      max <- 7*sqrt(lambda.max)
#      min <- -7*sqrt(lambda.max)
#      tmp.x <- seq(min,max,length=1000)
#      my.x <- NULL
#      for(k in 1:k.size){
#	  samp <- sample(1000,20)
#        my.x <- rbind(my.x,tmp.x[samp])
#      }

#      est.points <- my.x
#      width <- diff(tmp.x)
      
#      tmp <- 0
#      for(i in 1:(ncol(est.points)-1)){
#        tmp <- tmp + gz_pi1(est.points[,i])*(width[i]^k.size)
#      }

##C³
      matA <- eigen(Sigma)$vector

      gz_pi1 <- function(z){
	  tmpz <- t(matA) %*% z
        tmp <- G(tmpz) * get.pi1(tmpz)	#det(matA) = 1
        return( tmp  )
      }
 
      my.x <- NULL

      for(k in 1:k.size){
        max <- 7 * sqrt(lambda[k])
        min <- -7 * sqrt(lambda[k])
        tmp.x <- seq(min,max,length=1000)
	  samp <- sample(1000,20)
        my.x <- rbind(my.x,tmp.x[samp])
      }

      est.points <- my.x
      
      tmp <- 0
      for(i in 1:20){
        tmp <- tmp + gz_pi1(est.points[,i])
      }

	tmp <- tmp/20


##    }else if( 2 <= k.size || k.size <= 20 ){ # use 'cubatur()' to solve multi-dim integration.
    }else if( (2 <= k.size || k.size <= 20) && FALSE ){ # use 'cubatur()' to solve multi-dim integration. 
      max <- 10*lambda.max
##      max <- 10*lambda.max/k.size
      min <- -10*lambda.max
	  if(require(cubature)){
#		  tmp <- adapt(ndim=k.size,lower=rep(min,k.size),upper=rep(max,k.size),functn=gz_pi1)$value
        print("Messages bellow are for debug...")
        print(date())
        tmp <- adaptIntegrate(gz_pi1,lower=rep(min,k.size),upper=rep(max,k.size),tol=1e-5,maxEval=500/k.size)
        print(date())
        print("tolerance")
        print(tmp$error)
        print("function evaluated times")
        print(tmp$functionEvaluations)
        print("return code: if it is not 0, error occured in integration (may be does not converge)")
        print(tmp$returnCode)
        tmp <- tmp$integral
	  } else {
	    tmp <- NA	  
	  }
    }else{
      stop("length k is too long.")
    }
    return(tmp)
  }

  # 'rho' is a given function of (X,e) (p.7 formula(13.19))
  
  # d(rho)/di
  get.di.rho <- function(){
    di.rho <- numeric(d.size)
    tmp <- matrix(0,d.size,block)
    for(i in 1:d.size){
      di.rho[i] <- deriv(rho,state[i])
    }
    for(t in 1:block){
      for(i in 1:d.size){
        assign(state[i],X.t0[t,i])
      }
      for(j in 1:d.size){
        tmp[j,t]<- attr(eval(di.rho[j]),"grad")
      }
    }
    return(tmp)
  }
  
  # d(rho)/de
  get.de.rho <- function(){
    tmp <- matrix(0,1,block)
    de.rho <- deriv(rho,pars[1])
    for(t in 1:block){
      for(i in 1:d.size){
        assign(state[i],X.t0[t,i])
      }
      tmp[1,t]<- attr(eval(de.rho),"grad")
    }
    return(tmp)
  }
  
  # H_{t}^{e} at t=0 (see formula (13.19))
  # use trapezium integration
  get.H0 <- function(){

    assign(pars[1],0)

    tmp <- matrix(0,1,block)
    for(t in 1:block){
      for(i in 1:d.size){
        assign(state[i],X.t0[range[t],i])
      }
      tmp[1,t] <- eval(rho)
    }
    H0 <- exp( -(tmp %*% Diff[block,] *delta) )
    return(as.double(H0) )
  }


  #################################################
  # Here is an entry point of 'asymptotic_term()' #
  #################################################

  ## initialization part

  # make expressions of derivation of V0
  dx.drift <- Derivation.vector(V0,state,d.size,d.size)
  de.drift <- Derivation.scalar(V0,pars,d.size)
  dede.drift <- Derivation.scalar(de.drift,pars,d.size)  ##

  # make expressions of derivation of V
  dx.diffusion <- as.list(NULL)
  for(i in 1:d.size){
    dx.diffusion[i] <- list(Derivation.vector(V[[i]],state,d.size,r.size))
  }
  de.diffusion <- as.list(NULL)
  for(i in 1:d.size){
    de.diffusion[i] <- list(Derivation.scalar(V[[i]],pars,r.size))
  }
  dede.diffusion <- as.list(NULL)
  for(i in 1:d.size){
    dede.diffusion[i] <- list(Derivation.scalar(de.diffusion[[i]],pars,r.size))
  }

  yuima.warn("Get variables ...")
  Y <- Get.Y()
  mu <- funcmu()
  aMat <- funca()   ## ¤³¤³¤Ç¥¨¥é¡¼ 2010/11/24, TBC
  Sigma <- funcsigma()

  # calculate each variables shown in p.9 and
  # prepare for trapezium integration
  
  range <- make.range.for.trapezium.fomula(division+1,(block-1))
  lambda.ts <- matrix(0,(d.size*block),(r.size*block))
  lambda.ts0 <- matrix(0,(block*d.size),block)
  nu.ijts <- array(0,dim=c(d.size,d.size,(d.size*block),block))
  nu.ts <- matrix(0,(d.size*block),(r.size*block))
  nu.ts0 <- matrix(0,(block*d.size),block)
  mu.its0 <- array(0,dim=c(d.size,(d.size*block),block))
  mu.its <- array(0,dim=c(d.size,(d.size*block),(r.size*block)))
  ytyis <- array(0,dim=c(length(range),length(range),d.size,d.size))
  for(t in 1:length(range)) {
    for(s in 1:length(range)) {
	  ytyis[t,s,,] <- Get.YtYis(t,s,range)
	}
  }
  # calculate width of trapezoid
  diff <- numeric(block-1)

##C³
#  diff[1:(block-2)] <- range[2] - range[1]
#  diff[block-1] <- range[block]-range[block-1]
  diff <- range[2:block] - range[1:(block-1)]	#diff‚Ì·‚Íˆê’è‚Å‚Í‚È‚¢
##

  Diff <- matrix(0,block,block)
  if(k.size == 1){
    aMat.tmp <- t(as.matrix(aMat[,,1]))
  }else{
    aMat.tmp <- cbind(aMat[,,1])
  }    
  for(t in 2:length(range)){
    if(k.size==1){
      tmp.matrix <- t(as.matrix(aMat[,,range[t]]))
      aMat.tmp <- cbind(aMat.tmp, tmp.matrix)
    }else{  
      aMat.tmp <- cbind(aMat.tmp, aMat[,,range[t]])
    }
  
    dbase <- (t-1)*d.size
    for(s in 1:t){
      rbase <- (s-1)*r.size
      if( s==1 || t==s ){
        if(s==1){
          Diff[t,s] <- diff[s]
        }else if(t==s){
          Diff[t,s] <- diff[s-1]
        }
        lambda.ts0[(dbase+1):(dbase+d.size),s] <- Get.lambda.ts0(t,s,range)
        nu.ts0[(dbase+1):(dbase+d.size),s] <- Get.nu.ts0(t,s,range)
        lambda.ts[(dbase+1):(dbase+d.size),(rbase+1):(rbase+r.size)] <- Get.lambda.ts(t,s,range)
        nu.ts[(dbase+1):(dbase+d.size),(rbase+1):(rbase+r.size)] <- Get.nu.ts(t,s,range)
      }else{
        Diff[t,s] <- diff[s-1]+diff[s]
        lambda.ts0[(dbase+1):(dbase+d.size),s] <- Get.lambda.ts0(t,s,range)
        nu.ts0[(dbase+1):(dbase+d.size),s] <- Get.nu.ts0(t,s,range)        
        lambda.ts[(dbase+1):(dbase+d.size),(rbase+1):(rbase+r.size)] <- Get.lambda.ts(t,s,range)
        nu.ts[(dbase+1):(dbase+d.size),(rbase+1):(rbase+r.size)] <- Get.nu.ts(t,s,range)
      }
    }
  }
  
  Diff<- Diff / 2 # trapezium integrations need "/2"
  for( t in 1:length(range)){
    dbase <- (t-1)*d.size
    for( s in 1:t){
      rbase <- (s-1)*r.size
      for( i in 1:d.size){
        mu.its[i,( dbase+1 ):( dbase+d.size ),( (rbase+1):(rbase+r.size) )] <- Get.mu.its(i,t,s,range)
		mu.its0[i,(dbase+1):(dbase+d.size),s] <- Get.mu.its0(i,t,s,range)
        for(j in 1:d.size){
          nu.ijts[i,j,(dbase+1):(dbase+d.size) , s] <- Get.nu.ijts(i,j,t,s,range)
        }
      }
    }
  }
  ## end of Get variables

  ## derivations
  func1 <- deriv.f0.for.state(f[[1]])
  func2 <- hessian.f0.di.dj(f[[1]])
  func3 <- hessian.f0.di.de(f[[1]])
  func5 <- hessian.f.dxi.de(f)
  func7 <- deriv.f0.for.state(F)
  func8 <- hessian.f0.di.dj(F)
  func9 <- hessian.f0.di.de(F)
  di.f0.l <- array(0,dim=c(k.size,d.size,block))
  di.dj.f0.l <- array(0,dim=c(k.size,d.size^2,block))
  di.de.f0.l <- array(0,dim=c(k.size,d.size,block))
  de.de.f0.l <- matrix(0,k.size,block)
  di.de.f.l <- array(0,dim=c(k.size,d.size,r.size*block))
  de.de.f.l <- matrix(0,k.size,r.size*block)
  di.F.l <- matrix(0,k.size,d.size)
  di.dj.F.l <- matrix(0,k.size,d.size^2)
  di.de.F.l <- matrix(0,k.size,d.size)
  de.de.F.l <- numeric(k.size)
  for(k in 1:k.size){
    for(t in 1:block){
      di.f0.l[k,,t] <- input.deriv(range[t],expr=func1,l=k,X0=X.t0)
      di.dj.f0.l[k,,t] <- input.hessian.dxi.dxj(range[t],expr=func2,l=k,X0=X.t0)
      di.de.f0.l[k,,t] <- input.hessian(range[t],expr=func3,n=1,m=2,l=k,X0=X.t0)
      de.de.f0.l[k,t] <- input.hessian(t,expr=func3,n=2,m=2,l=k,X0=X.t0)[1]
      di.de.f.l[k,,((t-1)*r.size+1):((t-1)*r.size+r.size)] <- input.hessian.dxi.de.f(range[t],expr=func5,l=k,X0=X.t0)
      de.de.f.l[k,((t-1)*r.size+1):((t-1)*r.size+r.size)] <- input.hessian.dxi.de.f(range[t],expr=func5,n=2,m=2,l=k,X0=X.t0)[1,]
    }
    di.F.l[k,] <- input.deriv(block,expr=func7,l=k,X0=X.t0)
    di.dj.F.l[k,] <- input.hessian.dxi.dxj(block,expr=func8,l=k,X0=X.t0)
    di.de.F.l[k,] <- input.hessian(block,expr=func9,n=1,m=2,l=k,X0=X.t0)
    de.de.F.l[k] <- input.hessian(block,expr=func9,n=2,m=2,l=k,X0=X.t0)[1]
  }
  di.rho <- get.di.rho()
  de.rho <- get.de.rho()
  H0 <- get.H0()

  yuima.warn("Done.")
  yuima.warn("Initializing ...")
  dat.di.bar <- get.di.bar.init()
  dat.dij <- get.dij.init()
  dat.el <- get.el.init()
  yuima.warn("Done.")
  
  ## calculation
  yuima.warn("Calculating d0 ...")
  d0 <- get.d0.term()

  yuima.warn("Done\n")
  yuima.warn("Calculating d1 term ...")
  d1 <- get.d1.term()
  yuima.warn("Done\n")
  terms <- list(d0=d0, d1=d1)
  return(terms)
 })
