euler<-function(xinit,yuima,dW,env){

	
	sdeModel<-yuima@model
	
	modelstate <- sdeModel@solve.variable
	modeltime <- sdeModel@time.variable
	V0 <- sdeModel@drift
	V <- sdeModel@diffusion
	r.size <- sdeModel@noise.number
	d.size <- sdeModel@equation.number
	Terminal <- yuima@sampling@Terminal[1]
	division <- yuima@sampling@division[1]
	
	dX <- xinit
	
##:: set time step	
	delta <- Terminal/division 
	
##:: check if DRIFT and/or DIFFUSION has values
	has.drift <- sum(as.character(sdeModel@drift) != "(0)")
	var.in.diff <- is.logical(any(match(unlist(lapply(sdeModel@diffusion, all.vars)), sdeModel@state.variable)))
	
  ##:: function to calculate coefficients of dW(including drift term)
  ##:: common used in Wiener and CP
  p.b <- function(t, X=numeric(d.size)){
    ##:: assign names of variables
    for(i in 1:length(modelstate)){
      assign(modelstate[i], X[i], env)
    }
    assign(modeltime, t, env)
    ##:: solve diffusion term
    if(has.drift){
      tmp <- matrix(0, d.size, r.size+1)
      for(i in 1:d.size){
        tmp[i,1] <- eval(V0[i], env)
        for(j in 1:r.size){
          tmp[i,j+1] <- eval(V[[i]][j],env)
        }
      }
    }else{  ##:: no drift term (faster)
      tmp <- matrix(0, d.size, r.size)
      for(i in 1:d.size){
        for(j in 1:r.size){
          tmp[i,j] <- eval(V[[i]][j],env)
        }
      }
    }
    return(tmp)
  }
  
  X_mat <- matrix(0, d.size, division+1)
  X_mat[,1] <- dX
  
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
          pbdata[(i-1)*(r.size+1)+1, ] <- eval(V0[i], env)
          for(j in 1:r.size){
            pbdata[(i-1)*(r.size+1)+j+1, ] <- eval(V[[i]][j], env)
          }
        }
        dim(pbdata)<-(c(r.size+1, d.size*t.size))
      }else{
        pbdata <- matrix(0, d.size*r.size, t.size)
        for(i in 1:d.size){
          for(j in 1:r.size){
            pbdata[(i-1)*r.size+j, ] <- eval(V[[i]][j], env)
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
        assign(modelstate[i], X[i], env)
      }
      assign(modeltime, t, env)
      tmp <- numeric(d.size)
      for(i in 1:d.size){
        tmp[i] <-  eval(JP[i], env)
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
      if(grep("^[dexp|dnorm|dgamma]", sdeModel@measure$df$expr)){
        ##:: e.g. dnorm(z,1,1) -> rnorm(mu.size*N_sharp,1,1)
        F <- suppressWarnings(parse(text=gsub("^d(.+?)\\(.+?,", "r\\1(mu.size*N_sharp,", sdeModel@measure$df$expr, perl=TRUE)))
      }else{
        stop("Sorry. CP only supports dexp, dnorm and dgamma yet.")
      }
      randJ <- eval(F)  ## this expression is evaluated locally not in the yuimaEnv
      j <- 1
      for(i in 1:division){
        if(JAMP==FALSE || sum(i==ij)==0){
          Pi <- 0
        }else{
          J <- eta0*randJ[j]/lambda
          j <- j+1
          ##cat(paste(J,"\n"))
          ##Pi <- zeta(dX, J)
          assign(sdeModel@jump.variable, J, env)
          
          if(sdeModel@J.flag){
            J <- 1
          }
          
          Pi <- p.b.j(t=i*delta,X=dX) * J
          ##Pi <- p.b.j(t=i*delta, X=dX)
        }
        dX <- dX + p.b(t=i*delta, X=dX) %*% dW[, i] + Pi
        X_mat[, i+1] <- dX
      }
      tsX <- ts(data=t(X_mat), deltat=delta, start=0)
      ##::end CP
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
        assign(sdeModel@jump.variable, dZ[i], env)
        
        if(sdeModel@J.flag){
          dZ[i] <- 1
        }
          
        dX <- dX + p.b(t=i*delta, X=dX) %*% dW[, i] +p.b.j(t=i*delta, X=dX) * dZ[i]
        X_mat[, i+1] <- dX
      }
      tsX <- ts(data=t(X_mat), deltat=delta, start=0)
      ##::end code
    }else{
      cat(paste("Type \"", sdeModel@measure.type, "\" not supported yet.\n", sep=""))
      return(NULL)
    }
  }##::end levy
  yuimaData <- setData(original.data=tsX)
  return(yuimaData)


}
