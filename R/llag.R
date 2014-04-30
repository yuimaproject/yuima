#lead-lag estimation

#x:data

#setGeneric( "llag", function(x,verbose=FALSE) standardGeneric("llag") )
#setMethod( "llag", "yuima", function(x,verbose=FALSE) llag(x@data,verbose=verbose ))
#setMethod( "llag", "yuima.data", function(x,verbose=FALSE) {

#setGeneric( "llag",
#		function(x,from=FALSE,to=FALSE,division=FALSE,verbose=FALSE)
#		standardGeneric("llag") )
#setMethod( "llag", "yuima",
#		function(x,from=FALSE,to=FALSE,division=FALSE,verbose=FALSE)
#		llag(x@data,from=FALSE,to=FALSE,division=FALSE,verbose=verbose ))
#setMethod( "llag", "yuima.data", function(x,from=FALSE,to=FALSE,division=FALSE,verbose=FALSE) {

setGeneric( "llag",
		function(x,from=-Inf,to=Inf,division=FALSE,verbose=FALSE,grid)
		standardGeneric("llag") )
setMethod( "llag", "yuima",
		function(x,from=-Inf,to=Inf,division=FALSE,verbose=FALSE,grid)
		llag(x@data,from=from,to=to,division=division,verbose=verbose,grid=grid))
setMethod( "llag", "yuima.data", function(x,from=-Inf,to=Inf,division=FALSE,
                                          verbose=FALSE,grid) {
  
  if((is(x)=="yuima")||(is(x)=="yuima.data")){
    zdata <- get.zoo.data(x)
  }else{
    print("llag:invalid argument")
    return(NULL)
  }
  
  d <- length(zdata)
  
  # allocate memory
  ser.times <- vector(d, mode="list") # time index in 'x'
  ser.diffX <- vector(d, mode="list") # difference of data
  vol <- double(d)
  
  for(i in 1:d){
    
    # NA data must be skipped
    idt <- which(is.na(zdata[[i]]))
    if(length(idt>0)){
      zdata[[i]] <- zdata[[i]][-idt]
    }
    if(length(zdata[[i]])<2) {
      stop("length of data (w/o NA) must be more than 1")
    }
    
    # set data and time index
    ser.times[[i]] <- as.numeric(time(zdata[[i]]))
    # set difference of the data 
    ser.diffX[[i]] <- diff( as.numeric(zdata[[i]]) )
    vol[i] <- sum(ser.diffX[[i]]^2)
  }
  
  theta <- matrix(0,d,d)
  covmat <- diag(vol)
  
  d.size <- d*(d-1)/2
  crosscov <- vector(d.size,mode="list")
  
  if(missing(grid)){
    
    if(length(from) != d.size){
      from <- c(from,rep(-Inf,d.size - length(from)))
    }
    
    if(length(to) != d.size){
      to <- c(to,rep(Inf,d.size - length(to)))
    }
    
    if(length(division) == 1){
      division <- rep(division,d.size)
    }
    
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        
        time1 <- ser.times[[i]]
        time2 <- ser.times[[j]] 
        
        # calculate the maximum of correlation by substituting theta to lagcce
        
        #n:=2*length(data)
        
        num <- d*(i-1) - (i-1)*i/2 + (j-i)
        
        if(division[num]==FALSE){
          n <- round(2*max(length(time1),length(time2)))+1
        }else{
          n <- division[num]
        }
        
        # maximum value of lagcce
        
        tmptheta <- time2[1]-time1[1] # time lag (sec)
        
        num1 <- time1[length(time1)]-time1[1] # elapsed time for series 1
        num2 <- time2[length(time2)]-time2[1] # elapsed time for series 2
        
        # modified
        
        if(is.numeric(from[num])==TRUE && is.numeric(to[num])==TRUE){
          num2 <- min(-from[num],num2+tmptheta)
          num1 <- min(to[num],num1-tmptheta)
          tmptheta <- 0
          
          if(-num2 >= num1){
            print("llag:invalid range")
            return(NULL)
          }
        }else if(is.numeric(from[num])==TRUE){
          num2 <- min(-from[num],num2+tmptheta)
          num1 <- num1-tmptheta
          tmptheta <- 0
          
          if(-num2 >= num1){
            print("llag:invalid range")
            return(NULL)
          }
        }else if(is.numeric(to[num])==TRUE){
          num2 <- num2+tmptheta
          num1 <- min(to[num],num1-tmptheta)
          tmptheta <- 0
          
          if(-num2 >= num1){
            print("llag:invalid range")
            return(NULL)
          }
        }
        
        y <- seq(-num2-tmptheta,num1-tmptheta,length=n)[2:(n-1)]
        
        tmp <- .C("HYcrosscov",
                  as.integer(n-2),
                  as.integer(length(time1)),
                  as.integer(length(time2)),
                  as.double(y),
                  as.double(time1),
                  as.double(time2),
                  double(length(time2)),
                  as.double(ser.diffX[[i]]),
                  as.double(ser.diffX[[j]]),
                  value=double(n-2),PACKAGE="yuima")$value
        
        idx <- which.max(abs(tmp))
        mlag <- -y[idx] # make the first timing of max or min
        cov <- tmp[idx]
        
        theta[i,j] <- mlag
        covmat[i,j] <- cov
        theta[j,i] <- -mlag
        covmat[j,i] <- covmat[i,j]
        
        crosscov[[num]] <- zoo(tmp,y)
      }
    }
    
  }else{
    
    if(!is.list(grid)){
      if(is.numeric(grid)){
        grid <- data.frame(matrix(grid,length(grid),d.size))
      }else{
        print("llag:invalid grid")
        return(NULL)
      }
    }
    
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        
        time1 <- ser.times[[i]]
        time2 <- ser.times[[j]] 
        
        num <- d*(i-1) - (i-1)*i/2 + (j-i)
        
        tmp <- .C("HYcrosscov",
                  as.integer(length(grid[[num]])),
                  as.integer(length(time1)),
                  as.integer(length(time2)),
                  as.double(grid[[num]]),
                  as.double(time1),
                  as.double(time2),
                  double(length(time2)),
                  as.double(ser.diffX[[i]]),
                  as.double(ser.diffX[[j]]),
                  value=double(length(grid[[num]])),PACKAGE="yuima")$value
        
        idx <- which.max(abs(tmp))
        mlag <- -grid[[num]][idx] # make the first timing of max or min
        cov <- tmp[idx]
        
        theta[i,j] <- mlag
        covmat[i,j] <- cov
        theta[j,i] <- -mlag
        covmat[j,i] <- covmat[i,j]
        
        crosscov[[num]] <- zoo(tmp,grid[[num]])
      }
    }
  }
  
  cormat <- diag(1/sqrt(diag(covmat)))%*%covmat%*%diag(1/sqrt(diag(covmat)))
  colnames(theta) <- names(zdata)
  rownames(theta) <- names(zdata)

  if(verbose==TRUE){
    colnames(covmat) <- names(zdata)
    rownames(covmat) <- names(zdata)
    colnames(cormat) <- names(zdata)
    rownames(cormat) <- names(zdata)
    return(list(lagcce=theta,covmat=covmat,cormat=cormat,crosscov=crosscov))
  }else{
    return(theta)
  }
})
  

## Old version
if(0){
setMethod( "llag", "yuima.data", function(x,from=-Inf,to=Inf,division=FALSE,verbose=FALSE) {
  
  
  if(!is(x)=="yuima.data"){
    if(is(x)=="yuima"){
      dat <- x@data
    }else{
      print("llag:invalid argument")
      return(NULL)
    }
  }else{
    dat <- x
  }
  
  d <- length(dat@zoo.data)
  
  lagccep <- function(datp,theta){
    time(datp[[2]]) <- time(datp[[2]])+theta
    return(cce(setData(datp))$covmat[1,2])
  }
  
  lagcce <- function(datzoo,theta){
    d <- dim(theta)[1]
    lcmat <- cce(setData(datzoo))$covmat
    
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        datp <- datzoo[c(i,j)]
        lcmat[i,j] <- lagccep(datp,theta[i,j])
        lcmat[j,i] <- lcmat[i,j]
      }
    }		
    return(lcmat)	
  }
  
  d.size <- d*(d-1)/2
  
  if(length(from) != d.size){
    from <- c(from,rep(-Inf,d.size - length(from)))
  }
  
  if(length(to) != d.size){
    to <- c(to,rep(Inf,d.size - length(to)))
  }
  
  if(length(division) != d.size){
    division <- c(division,rep(FALSE,d.size - length(division)))
  }
  
  find_lag <- function(i,j){
    datp <- dat@zoo.data[c(i,j)]
    time1 <- time(datp[[1]])
    time2 <- time(datp[[2]])  
    
    # calculate the maximum of correlation by substituting theta to lagcce
    
    #n:=2*length(data)
    
    num <- d*(i-1) - (i-1)*i/2 + (j-i)
    
    if(division[num]==FALSE){
      n <- round(2*max(length(datp[[1]]),length(datp[[2]])))+1
    }else{
      n <- division[num]
    }
    
    # maximum value of lagcce
    
    tmptheta <- as.numeric(time2[1])-as.numeric(time1[1]) # time lag (sec)
    
    num1 <- as.numeric(time1[length(time1)])-as.numeric(time1[1]) # elapsed time for series 1
    num2 <- as.numeric(time2[length(time2)])-as.numeric(time2[1]) # elapsed time for series 2
    
    # modified
    
    if(is.numeric(from[num])==TRUE && is.numeric(to[num])==TRUE){
      num2 <- min(-from[num],num2+tmptheta)
      num1 <- min(to[num],num1-tmptheta)
      tmptheta <- 0
      
      if(-num2 >= num1){
        print("llag:invalid range")
        return(NULL)
      }
    }else if(is.numeric(from[num])==TRUE){
      num2 <- min(-from[num],num2+tmptheta)
      num1 <- num1-tmptheta
      tmptheta <- 0
      
      if(-num2 >= num1){
        print("llag:invalid range")
        return(NULL)
      }
    }else if(is.numeric(to[num])==TRUE){
      num2 <- num2+tmptheta
      num1 <- min(to[num],num1-tmptheta)
      tmptheta <- 0
      
      if(-num2 >= num1){
        print("llag:invalid range")
        return(NULL)
      }
    }
    
    y <- seq(-num2-tmptheta,num1-tmptheta,length=n)
    tmp <- double(n)
    
    for(i.tmp in 2:(n-1)){
      tmp[i.tmp] <- lagccep(datp,y[i.tmp])
    }
    
    mat <- cbind(y[2:(n-1)],tmp[2:(n-1)])
    
    #idx <- abs(mat[,2])==max(abs(mat[,2]))
    #mlag <- mat[,1][idx][1] # make the first timing of max or min
    #cov <- mat[,2][idx][1]
    idx <- which.max(abs(mat[,2]))
    mlag <- mat[,1][idx] # make the first timing of max or min
    cov <- mat[,2][idx]
    return(list(lag=-mlag,cov=cov))
  }
  
  theta <- matrix(numeric(d^2),ncol=d)
  covmat <- cce(dat)$covmat
  
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      fl <- find_lag(i,j)
      theta[i,j] <- fl$lag
      covmat[i,j] <- fl$cov
      theta[j,i] <- -theta[i,j]
      covmat[j,i] <- covmat[i,j]
    }
  }
  
  
  #  covmat <- lagcce(dat@zoo.data,theta)
  cormat <- diag(1/sqrt(diag(covmat)))%*%covmat%*%diag(1/sqrt(diag(covmat)))
  if(verbose==TRUE){
    return(list(lagcce=theta,covmat=covmat,cormat=cormat))
  }else{
    return(theta)
  }
})
}