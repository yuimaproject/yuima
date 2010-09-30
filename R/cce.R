# cce() Cumulative Covariance Estimator

#
# CumulativeCovarianceEstimator
#

# returns a matrix of var-cov estimates

setGeneric( "cce", function(x) standardGeneric("cce") )

setMethod("cce", "yuima", function(x) cce(x@data) )

setMethod("cce", "yuima.data", function(x) {  
	
	data <- get.zoo.data(x)
	n.series <- length(data)
#if(n.series <2)
# stop("Please provide at least 2-dimensional time series")
	
# allocate memory
	ser.X <- vector(n.series, mode="list")     # data in 'x'
	ser.times <- vector(n.series, mode="list") # time index in 'x'
	ser.diffX <- vector(n.series, mode="list") # difference of data
	
	for(i in 1:n.series){
# set data and time index
		ser.X[[i]] <- as.numeric(data[[i]]) # we need to transform data into numeric to avoid problem with duplicated indexes below
		ser.times[[i]] <- as.numeric(time(data[[i]]))
		
# NA data must be skipped
		idt <- which(is.na(ser.X[[i]]))
		if(length(idt>0)){
			ser.X[[i]] <- (ser.X[[i]])[-idt]
			ser.times[[i]] <- (ser.times[[i]])[-idt]
		}
		if(length(ser.X[[i]])<2) {
			stop("length of data (w/o NA) must be more than 1")
		}
# set difference of the data
		ser.diffX[[i]] <- diff( ser.X[[i]] )
	}
	
	
# core part of cce
	
	cmat <- matrix(0, n.series, n.series)  # cov
	for(i in 1:n.series){
		for(j in i:n.series){ 
			I <- rep(1,n.series)
#Checking Starting Point
			repeat{
				if(ser.times[[i]][I[i]] >= ser.times[[j]][I[j]+1]){
					I[j] <- I[j]+1	 
				}else if(ser.times[[i]][I[i]+1] <= ser.times[[j]][I[j]]){
					I[i] <- I[i]+1	 
				}else{
					break
				}
			}
			
			
#Main Component
			if(i!=j){
				while((I[i]<length(ser.times[[i]])) && (I[j]<length(ser.times[[j]]))) {
					cmat[j,i] <- cmat[j,i] + (ser.diffX[[i]])[I[i]]*(ser.diffX[[j]])[I[j]]
					if(ser.times[[i]][I[i]+1]>ser.times[[j]][I[j]+1]){
						I[j] <- I[j] + 1
					}else if(ser.times[[i]][I[i]+1]<ser.times[[j]][I[j]+1]){
						I[i] <- I[i] +1
					}else{
						I[i] <- I[i]+1
						I[j] <- I[j]+1
					}
				}
			}else{
				cmat[i,j] <- sum(ser.diffX[[i]]^2)
			}
			cmat[i,j] <- cmat[j,i]  
		}
		
	}
	return( cmat )            
})

##:: add 2010/09/19 for list data handling
setMethod("cce", "list", function(x){
  ##:: init
  cce.result <- list()
  
  ##:: error check and  run cce
  for( i in 1:length(x) ){
    if(class(x[[i]])=="yuima" || class(x[[i]])=="yuima.data"){
      cce.result[[i]] <- cce(x[[i]])
    }else{
      stop("Objects in list-class must be yuima-class or yuima.data-class.")
    }
  }  

  ##:: return result
  names(cce.result) <- names(x)
  return(cce.result)
  
})
