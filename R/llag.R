#lead-lag estimation

#x:data

#setGeneric( "llag", function(x,verbose=FALSE) standardGeneric("llag") )
#setMethod( "llag", "yuima", function(x,verbose=FALSE) llag(x@data,verbose=verbose ))
#setMethod( "llag", "yuima.data", function(x,verbose=FALSE) {

setGeneric( "llag",
		function(x,from=FALSE,to=FALSE,division=FALSE,verbose=FALSE)
		standardGeneric("llag") )
setMethod( "llag", "yuima",
		function(x,from=FALSE,to=FALSE,division=FALSE,verbose=FALSE)
		llag(x@data,from=FALSE,to=FALSE,division=FALSE,verbose=verbose ))
setMethod( "llag", "yuima.data", function(x,from=FALSE,to=FALSE,division=FALSE,verbose=FALSE) {


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

	if(length(division) == 1){
		division <- rep(division,d.size)
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
		tmp <- real(n)

		for(i in 2:(n-1)){
			tmp[i] <- lagccep(datp,y[i])
		}
  
		mat <- cbind(y[2:(n-1)],tmp[2:(n-1)])
		
		idx <- abs(mat[,2])==max(abs(mat[,2]))
		mlag <- mat[,1][idx][1] # make the first timing of max or min
		cov <- mat[,2][idx][1]
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
		return(list(lagcce=theta,covmat=covmat,cormat=cormat))
	}
})
  
