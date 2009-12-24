
##Constructor and Initializer of class 'sampling'

# we convert objects to "zoo" internally

 
setMethod("initialize", "yuima.sampling",
           function(.Object, Initial, Terminal, n, delta, grid, random, 
				 regular, sdelta, sgrid, oindex, interpolation){  
				.Object@sdelta <- as.numeric(NULL) 	 
				.Object@sgrid <- as.numeric(NULL) 	 
				.Object@oindex <- as.numeric(NULL) 	 
				.Object@interpolation <- interpolation 	 
			   if(length(grid)>0){
				   testInitial<-(min(grid)==Initial)
				   testTerminal<-(max(grid)==Terminal)
				   testn<-(abs(n-diff(range(grid))/mean(diff(grid))+1)<10^(-10))
				   testdelta<-(abs(delta-mean(diff(grid)))<10^(-10))
				   testregular<-all(abs(diff(diff(grid)))<10^(-10))
				   
				   if(!testInitial){
					   cat("\n Start time has been set with the grid \n")
				   }
				   if(!testTerminal){
					   cat("\n Terminal time has been set with the grid \n")
				   }
				   if(!testn){
					   cat("\n Division has been set with the grid \n")
				   }
				   if(!testdelta){
					   cat("\n delta has been set with the grid \n")
				   }
				   if(testregular){
				.Object@n <- diff(range(grid))/mean(diff(grid))+1
				.Object@delta    <- mean(diff(grid))
			    .Object@regular  <- TRUE
			       }else{
			    .Object@n <- length(grid)-1
			    .Object@delta    <- as.numeric(NULL)
			    .Object@regular  <- FALSE  
				   }
				.Object@grid     <- grid
				.Object@Initial  <- min(grid)
				.Object@Terminal <- max(grid)   
				.Object@random   <- random   
			   }else{ 
			 # There is no grid
			 eqn <- length(Terminal)
             if(length(Terminal)==length(n)){
               .Object@Initial  <- Initial
			   .Object@Terminal <- Terminal
               .Object@n <- n
			   .Object@delta    <- delta
			   .Object@grid     <- grid
			   .Object@random   <- random
			   .Object@regular  <- regular
             }else if(length(Terminal)==1){               
               .Object@Initial  <- Initial
			   .Object@Terminal <- rep(Terminal, length(n))
			   .Object@n <- n	 
			   .Object@delta    <- delta
			   .Object@grid     <- grid
			   .Object@random   <- random
			   .Object@regular  <- regular	 
             }else if(length(n)==1){
               .Object@Initial  <- Initial
			   .Object@Terminal <- Terminal
               .Object@n <- rep(n, length(Terminal))
			   .Object@delta    <- delta
			   .Object@grid     <- grid
			   .Object@random   <- random
			   .Object@regular  <- regular
             }else{
               cat("\nDimension missmatch.\n")
               return(NULL)
             }}
             return(.Object)
           })

setSampling <- function(Initial=0, Terminal=1, n=100, delta=0.1, 
 grid=as.numeric(NULL), random=FALSE, sdelta=as.numeric(NULL), 
 sgrid=as.numeric(NULL), interpolation="pt" ){
  return(new("yuima.sampling", Initial=Initial, Terminal=Terminal, 
	n=n, delta=delta, grid=grid, random=random, 
			 regular=TRUE, sdelta=sdelta, sgrid=sgrid,
			 interpolation=interpolation))
}

