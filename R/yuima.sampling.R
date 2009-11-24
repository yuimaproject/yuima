
##Constructor and Initializer of class 'sampling'

# we convert objects to "zoo" internally

 
setMethod("initialize", "yuima.sampling",
           function(.Object, Terminal, division, Initial, grid, random){  
             eqn <- length(Terminal)
             if(length(Terminal)==length(division)){
               .Object@Terminal <- Terminal
               .Object@division <- division
			   .Object@Initial  <- Initial
			   .Object@grid     <- grid
			   .Object@random   <- random 
             }else if(length(Terminal)==1){
               .Object@division <- division
               .Object@Terminal <- rep(Terminal, length(division))
			   .Object@Initial  <- Initial
			   .Object@grid     <- grid
			   .Object@random   <- random 
             }else if(length(division)==1){
               .Object@Terminal <- Terminal
               .Object@division <- rep(division, length(Terminal))
			   .Object@Initial  <- Initial
			   .Object@grid     <- grid
			   .Object@random   <- random 
             }else{
               cat("\nDimension missmatch.\n")
               return(NULL)
             }
             return(.Object)
           })

setSampling <- function(Terminal, division, Initial=0, grid=as.numeric(NULL), random=FALSE){
  return(new("yuima.sampling", Terminal=Terminal, division=division, Initial=Initial, grid=grid, random=random))
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
