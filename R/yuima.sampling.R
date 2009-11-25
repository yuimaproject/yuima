
##Constructor and Initializer of class 'sampling'

# we convert objects to "zoo" internally

 
setMethod("initialize", "yuima.sampling",
           function(.Object, Initial, Terminal, division, delta, grid, random, regular){  
             
			   if(length(grid)>0){
				   
				   testInitial<-(min(grid)==Initial)
				   testTerminal<-(max(grid)==Terminal)
				   testdivision<-(abs(division-diff(range(grid))/mean(diff(grid))+1)<10^(-10))
				   testdelta<-(abs(delta-mean(diff(grid)))<10^(-10))
				   testregular<-all(abs(diff(diff(grid)))<10^(-10))
				   
				   if(!testInitial){
					   cat("\n Start time has been set with the grid \n")
				   }
				   if(!testTerminal){
					   cat("\n Terminal time has been set with the grid \n")
				   }
				   if(!testdivision){
					   cat("\n Division has been set with the grid \n")
				   }
				   if(!testdelta){
					   cat("\n delta has been set with the grid \n")
				   }
				   if(testregular){
				.Object@division <- diff(range(grid))/mean(diff(grid))+1
				.Object@delta    <- mean(diff(grid))
			    .Object@regular  <- TRUE
			       }else{
			    .Object@division <- length(grid)-1
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
             if(length(Terminal)==length(division)){
               .Object@Initial  <- Initial
			   .Object@Terminal <- Terminal
               .Object@division <- division
			   .Object@delta    <- delta
			   .Object@grid     <- grid
			   .Object@random   <- random
			   .Object@regular  <- regular
             }else if(length(Terminal)==1){               
               .Object@Initial  <- Initial
			   .Object@Terminal <- rep(Terminal, length(division))
			   .Object@division <- division	 
			   .Object@delta    <- delta
			   .Object@grid     <- grid
			   .Object@random   <- random
			   .Object@regular  <- regular	 
             }else if(length(division)==1){
               .Object@Initial  <- Initial
			   .Object@Terminal <- Terminal
               .Object@division <- rep(division, length(Terminal))
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

setSampling <- function(Initial=0, Terminal=1, division=100, delta=0.1, grid=as.numeric(NULL), random=FALSE){
  return(new("yuima.sampling", Initial=Initial, Terminal=Terminal, division=division, delta=delta, grid=grid, random=random, regular=TRUE))
}

# accessors
## setGeneric("setSampling",
##            function(x,Terminal,division)
##            standardGeneric("setSampling")
##            )

## setMethod("setSampling", signature(x="yuima.model"),
##           function(x,Terminal=1,division=100){
##             if(length(x@equation.number) == 0) {
##               cat("model is empty!\n")
##               return(NULL)
##             }
##             d.size <- x@equation.number
##             if(length(Terminal)==1) {
##               Terminal <- rep(Terminal,d.size)
##             }
##             if(length(division)==1) {
##               division <- rep(division,d.size)
##             }
##             if(length(Terminal)!=d.size) {
##               cat("length of Terminal is not matched with model!\n")
##               return(NULL)
##             }
##             if(length(division)!=d.size) {
##               cat("length of division is not matched with model!\n")
##               return(NULL)
##             }
##             return(new("yuima.sampling", Terminal=Terminal, division=division ))
##           })

## setMethod("setSampling", signature(x="yuima.data"),
##           function(x,Terminal=1,division=100){
## 		    d.size <- dim(get.zoo.data(x))[2]
## 		    if(is.null(d.size)) {
## 			  cat("data is empty!\n")
## 			  return(NULL)
## 			}
			
## 			if(length(Terminal)==1) {
## 			  Terminal <- rep(Terminal,d.size)
## 			}
## 			if(length(division)==1) {
## 			  division <- rep(division,d.size)
## 			}
## 			if(length(Terminal)!=d.size) {
## 			  cat("length of Terminal is not matched with model!\n")
## 			  return(NULL)
## 			}
## 			if(length(division)!=d.size) {
## 			  cat("length of division is not matched with model!\n")
## 			  return(NULL)
## 			}
##             return(new("yuima.sampling", Terminal=Terminal, division=division ))
##           }) 
		  
## setMethod("setSampling", signature(x="yuima"),
##           function(x,Terminal=1,division=100){
## 		    if(length(x@model@EQNS) == 0) {
## 			  return(setSampling(x@data,Terminal,division))
## 			}
## 			return(setSampling(x@model,Terminal,division))
##           }) 
