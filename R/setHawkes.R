setHawkes <- function(lower.var="0", upper.var = "t", var.dt = "s",
  process = "N", dimension = 1, intensity = "lambda",
  ExpKernParm1="c", ExpKernParm2 ="a",
  const = "nu", measure = NULL, measure.type = NULL){
  if(dimension==1){
    if(is.null(measure)){
      my.rPPR1 <- function(n){
        res <- t(t(rep(1,n)))
        return(res)
      }


      my.dPPR1 <-function(x){
        res<-1
        return(res)
      }

      # yuima.law for the underlying source of randomness

      Law.PPR1 <- setLaw(rng = my.rPPR1, density = my.dPPR1)
      measure <- list(df = Law.PPR1)
      measure.type <- "code"
    }
  }else{
    if(is.null(measure)){
      yuima.stop("Missing argument measure: see setHawkes documentation")
    }
  }

  PROCESS <-  paste0(process,c(1:dimension))
  leng <- length(PROCESS)

  mod1 <- setModel(drift = rep("0",leng),
    diffusion = matrix("0",leng,leng),
    jump.coeff = diag("1",leng,leng),
    measure = measure, measure.type = measure.type,
    solve.variable = PROCESS)

  INTENSITY <- as.list(paste0(intensity,c(1:dimension)))

  gFun <- paste0(const,c(1:dimension))

  Ccoeff<-as.character(MatrCoeff(ExpKernParm1, dimension))
  Acoeff<-as.character(MatrCoeff(ExpKernParm2, dimension))
  Kernelpar<-c(Acoeff,Ccoeff)

  Kernel<- matrix(paste0(Ccoeff,"*exp(-",Acoeff,"*(","t-",var.dt,"))"),dimension, dimension)

  res <- aux.setPPR(yuima = mod1, counting.var=PROCESS,
    gFun, Kernel,
    var.dx = PROCESS, var.dt = var.dt, lambda.var = INTENSITY,
    lower.var=lower.var, upper.var = upper.var,
    nrow =dimension ,ncol=dimension, general = FALSE)

    return(res)
}

MatrCoeff<-function(lett, dimension){
  c1<-paste0(lett,c(1:dimension))

  cMatrix<-matrix(NA,dimension,dimension)
  for(i in c(1:dimension)){
    cMatrix[i,]<-paste0(c1[i],c(1:dimension))
  }
  return(cMatrix)
}
