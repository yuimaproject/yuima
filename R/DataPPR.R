# The following example is necessary for the structure of
# the PPR Data
# library(yuima)
# Terminal <- 10
# samp <- setSampling(T=Terminal,n=1000)
# 
# # Ex 1. (Simple homogeneous Poisson process)
# mod1 <- setPoisson(intensity="lambda", df=list("dconst(z,1)"))
# set.seed(123)
# y1 <- simulate(mod1, true.par=list(lambda=1),sampling=samp)
# y1@data@original.data
# simprvKern->yuimaPPR
get.counting.data<-function(yuimaPPR,type="zoo"){
  count <- yuimaPPR@PPR@counting.var
  dimCount <- length(count)
  Time_Arrivals_Grid <- index(yuimaPPR@data@zoo.data[[1]])
  if(dimCount==1){
    cond <- yuimaPPR@model@solve.variable %in% count 
    Count.Var <- yuimaPPR@data@original.data[,cond]
  }
  
  Time_Arrivals <- t(Time_Arrivals_Grid[1]) 
  colnames(Time_Arrivals) <- "Time"
  Obser <- t(yuimaPPR@data@original.data[1,])
  colnames(Obser) <- yuimaPPR@model@solve.variable
  cond <- diff(Count.Var)==1
  Time_Arrivals <- rbind(Time_Arrivals,as.matrix(Time_Arrivals_Grid[-1][cond]))
  if(is(yuimaPPR@data@original.data,"matrix")){
    Obser <- rbind(Obser,yuimaPPR@data@original.data[-1,][cond,])
  }else{
    Obser <- rbind(Obser,as.matrix(yuimaPPR@data@original.data[-1][cond]))
  }  
  
  Data <- zoo(x = Obser,order.by = Time_Arrivals) 
 # plot(Data)
  if(type=="zoo"){
    Data <- Data
    return(Data)
  }
  if(type=="yuima.PPR"){
    yuimaPPR@data@original.data <- Data
    return(yuimaPPR)
  }
  if(type=="matrix"){
    Data <- cbind(Time_Arrivals, Obser)
    return(Data)
  }
  yuima.stop("type is not supported.")
}

DataPPR <- function(CountVar, yuimaPPR, samp){
  if(!is(yuimaPPR,"yuima.PPR")){
    yuima.stop("...")
  }
  if(!is(samp,"yuima.sampling")){
    yuima.stop("...")
  }
  if(!is(CountVar,"zoo")){
    yuima.stop("...")
  }
  names <- colnames(CountVar)
  if(all(names %in% yuimaPPR@model@solve.variable) && all(yuimaPPR@model@solve.variable %in% names)){
    TimeArr <- index(CountVar)
    grid <- samp@grid[[1]]
    dimData <- dim(CountVar)[2]
    DataOr <- NULL
    for(i in c(1:dimData)){
      prv <- approx(x=TimeArr, y = CountVar[,i], xout = grid, method = "constant")$y
      DataOr <- cbind(DataOr,prv)
    }
    DataFin <- zoo(DataOr, order.by = grid)
    colnames(DataFin) <- names
    myData <- setData(DataFin)
    yuimaPPR@data <- myData
    yuimaPPR@sampling <- samp
    return(yuimaPPR)
  }else{
    dummy <- paste0("The labels of variables should be ", paste0(yuimaPPR@model@solve.variable,collapse= ", "),collapse = " ")
    yuima.stop(dummy)
  }  
}