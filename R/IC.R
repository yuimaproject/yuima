## information criteria

IC <- function(yuima, data = NULL, start, lower, upper, joint = FALSE, rcpp = FALSE,...){
  if(missing(yuima)) stop("yuima object is missing.")
  
  if(!is(yuima,"yuima")) stop("This function is for yuima-class.")
  
  state.num <- length(yuima@model@state.variable)
  if(is.null(yuima@data@zoo.data) == TRUE){
    if(is.matrix(data) == FALSE){
      n <- length(data)
      sub.zoo.data <- list(zoo(x = data, order.by = yuima@sampling@grid[[1]]))
      names(sub.zoo.data)[1] <- "Series 1"
    }else{
      sub.zoo.data <- list()
      if(ncol(data)-state.num != 0){
        data <- t(data)
      }
      n <- nrow(data)
      for(i in 1:ncol(data)){
        sub.zoo.data <- c(sub.zoo.data, list(zoo(x = data[,i], order.by = yuima@sampling@grid[[1]])))
        names(sub.zoo.data)[i] <- paste("Series", i)
      }
    }
    yuima@data@zoo.data <- sub.zoo.data
  }else{
    n <- yuima@sampling@n[1]
  }
  n <- n-1
  #alpha <- yuima@model@parameter@drift
  #beta <- yuima@model@parameter@diffusion
  Terminal <- yuima@sampling@Terminal[1]
  
  para.num.init  <- match(yuima@model@parameter@all, names(start))
  para.num.low  <- match(yuima@model@parameter@all, names(lower))
  para.num.upp  <- match(yuima@model@parameter@all, names(upper))
  para.start <- NULL
  para.lower <- NULL
  para.upper <- NULL
  for(j in 1:length(yuima@model@parameter@all)){
    para.start <- c(para.start, list(start[[para.num.init[j]]]))
    para.lower <- c(para.lower, list(lower[[para.num.low[j]]]))
    para.upper <- c(para.upper, list(upper[[para.num.upp[j]]]))
  }
  names(para.start) <- yuima@model@parameter@all
  names(para.lower) <- yuima@model@parameter@all
  names(para.upper) <- yuima@model@parameter@all
  
  mle <- qmle(yuima, start = para.start, lower = para.lower, upper = para.upper, method = "L-BFGS-B", joint = joint, rcpp = rcpp)
  hess <- list(mle@details$hessian)
  hess.diff <- subset(hess[[1]], rownames(hess[[1]])%in%yuima@model@parameter@diffusion, select=yuima@model@parameter@diffusion)
  hess.drif <- subset(hess[[1]], rownames(hess[[1]])%in%yuima@model@parameter@drift, select=yuima@model@parameter@drift)
  
  esti <- coef(mle)
  names(esti) <- c(yuima@model@parameter@diffusion, yuima@model@parameter@drift)
  cic <- summary(mle)@m2logL+2*(length(yuima@model@parameter@drift)+length(yuima@model@parameter@diffusion))
  bic <- summary(mle)@m2logL+length(yuima@model@parameter@drift)*log(Terminal)+length(yuima@model@parameter@diffusion)*log(n)
  if(det(hess.diff) > 0 && det(hess.drif) > 0){
    qbic <- summary(mle)@m2logL+log(det(hess.diff))+log(det(hess.drif))
  }else{
    qbic <- summary(mle)@m2logL+length(yuima@model@parameter@drift)*log(Terminal)+length(yuima@model@parameter@diffusion)*log(n)
  }
  
  final.res <- list(par = esti, BIC = bic, QBIC = qbic, CIC = cic)
  return(final.res)
}

