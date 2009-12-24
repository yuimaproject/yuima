
##:: function subsampling
##:: takes any yuima object with data or a yuima.data object and
##:: performs subsampling according to some method

# poisson.random.sampling
# returns sample of data using poisson sampling

setGeneric("subsampling", 
function(x, sampling=NULL, Initial, Terminal, delta, 
 grid = as.numeric(NULL), random = FALSE, sdelta=as.numeric(NULL), 
 sgrid=as.numeric(NULL), interpolation="none") 
 standardGeneric("subsampling")
)

setMethod("subsampling","yuima", 
function(x, sampling=NULL, Initial, Terminal, delta, 
 grid = as.numeric(NULL), random = FALSE, sdelta=as.numeric(NULL), 
 sgrid=as.numeric(NULL), interpolation="none")
 return(subsampling(x@data, sampling=sampling, Initial = Initial, 
  Terminal = Terminal, delta = delta, 
  grid = grid, random = random, sdelta=sdelta, 
  sgrid=sgrid, interpolation=interpolation))
)

setMethod("subsampling", "yuima.data",
function(x, sampling=sampling, Initial, Terminal, delta, 
 grid = as.numeric(NULL), random = FALSE, sdelta=as.numeric(NULL), 
 sgrid=as.numeric(NULL), interpolation="none"){

 tmpsamp <- NULL
 if(missing(sampling)){
	tmpsamp <- setSampling(Initial = Initial, Terminal = Terminal, 
				delta = delta, grid = grid, random = random, sdelta=sdelta, 
				sgrid=sgrid, interpolation=interpolation)
 } else {
	tmpsamp <- sampling
 }

 
 Data <- get.zoo.data(x)
 tmpgrid <- NULL

# random sampling
	 if(is.list(tmpsamp@random)){
		 rdist <- c(tmpsamp@random$rdist)
		 if(is.null(rdist))
			stop("provide at least `rdist' argument for random sampling")
		 n.rdist <- length(rdist)
		 n.data <- length(Data)
		 tmpgrid <- vector(n.data, mode="list")
		 r.gen <- rep( rdist, n.data) # eventually reciclying arguments
		 r.gen <- r.gen[1:n.data]
		 for(i in 1:n.data){
			 tmptime <- start(Data[[i]])	
			 T <- end(Data[[i]])
			 while(	sum( tmptime ) < T )
				tmptime <- c(tmptime, r.gen[[i]](1))
			 tmpgrid[[i]] <- cumsum(tmptime)
			 if(tail(tmpgrid[[i]],1)>T)
				tmpgrid[[i]] <- tmpgrid[[i]][-length(tmpgrid[[i]])]
		 }
	 }

	 otime <- vector(n.data, mode="list")
	 
	 for(i in 1:n.data){
	  otime[[i]] <- time(Data[[i]]) 
	  idx <- as.numeric(sapply(tmpgrid[[i]], function(x) max(which(otime[[i]] <= x))))
	  x@zoo.data[[i]] <- zoo(as.numeric(Data[[i]][idx]), order.by=tmpgrid[[i]])	 
	 }
	 return(x)
 } ### end method
)


