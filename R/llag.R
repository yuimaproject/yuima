#lead-lag estimation

#x:data plot:T or F
setGeneric( "llag", function(x) standardGeneric("llag") )
setMethod( "llag", "yuima", function(x) llag(x@data ))
setMethod( "llag", "yuima.data", function(x) {

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

##  dat <- get.zoo.data(x)
##  dat <- x@data

  data1 <- dat@zoo.data[[1]]
  data2 <- dat@zoo.data[[2]]
##  data1 <- as.numeric(dat[[1]])
##  data2 <- as.numeric(dat[[2]])
  time1 <- time(data1)
  time2 <- time(data2)  

  lagcce <- function(theta){
    stime <- time2 + theta  #shifted time
    time(dat@zoo.data[[2]]) <- stime
    return(cce(dat))
  }
  

# calculate the maximum of correlation by substituting theta to lagcce

#n:=2*length(data)

  n <- round(2*max(length(data1),length(data2)))+1

# maximum value of lagcce
  theta <- as.numeric(time2[1])-as.numeric(time1[1]) # time lag (sec)

  num1 <- as.numeric(time1[length(time1)])-as.numeric(time1[1]) # elapsed time for series 1
  num2 <- as.numeric(time2[length(time2)])-as.numeric(time2[1]) # elapsed time for series 2

  lagmax <- function(n){
    y <- seq(-num2-theta,num1-theta,length=n)
    tmp <- real(n)
    for(i in 2:(n-1)){
      tmp[i] <- lagcce(y[i])$covmat[1,2]
    }

    mat <- cbind(y[2:(n-1)],tmp[2:(n-1)])
    return(mat)
  }

 
  mat <- lagmax(n)

  opt <- mat[,1][abs(mat[,2])==max(abs(mat[,2]))][1] # make the first timing of max or min

  covmat <- lagcce(opt)$covmat
  cormat <- lagcce(opt)$cormat
  
  return(list(lagcce=opt,covmat=covmat,cormat=cormat,mat=mat))
})

