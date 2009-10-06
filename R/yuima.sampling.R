
##Constructor and Initializer of class 'sampling'

# we convert objects to "zoo" internally

 
setMethod("initialize", "yuima.sampling",
           function(.Object, Terminal, division){
             eqn <- length(Terminal)
             if(length(Terminal)==length(division)){
               .Object@Terminal <- Terminal
               .Object@division <- division
             }else if(length(Termianl)==1){
               .Object@division <- division
               .Object@Terminal <- rep(Terminal, length(division))
             }else if(length(division)==1){
               .Object@Terminal <- Terminal
               .Object@division <- rep(division, length(Terminal))
             }else{
               cat("\nDimension missmatch.\n")
               return(NULL)
             }
             return(.Object)
           })

setSampling <- function(Terminal, division){
  return(new("yuima.sampling", Terminal=Terminal, division=division))
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
