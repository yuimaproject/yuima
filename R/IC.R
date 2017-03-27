## information criteria

IC <- function(yuima, data, start, lower, upper, joint = FALSE, rcpp = FALSE,...){
  if(is(yuima,"yuima") == TRUE){
    
    if(is.matrix(data) == FALSE){
      sub.zoo.data <- list(zoo(x = data, order.by = yuima@sampling@grid[[1]]))
      names(sub.zoo.data)[1] <- "Series 1"
    }else{
      sub.zoo.data <- list()
      for(i in 1:ncol(data)){
        sub.zoo.data <- c(sub.zoo.data, list(zoo(x = data[,i], order.by = yuima@sampling@grid[[1]])))
        names(sub.zoo.data)[i] <- paste("Series", i)
      }
    }
    yuima@data@zoo.data <- sub.zoo.data
    alpha <- yuima@model@parameter@drift
    beta <- yuima@model@parameter@diffusion
    Terminal <- yuima@sampling@Terminal
    n <- yuima@sampling@n
    
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
    hess.diff <- subset(hess[[1]], rownames(hess[[1]])%in%beta, select=beta)
    hess.drif <- subset(hess[[1]], rownames(hess[[1]])%in%alpha, select=alpha)
    
    esti <- coef(mle)
    names(esti) <- c(beta, alpha)
    cic <- summary(mle)@m2logL+2*(length(alpha)+length(beta))
    bic <- summary(mle)@m2logL+length(alpha)*log(Terminal[1])+length(beta)*log(n[1])
    if(det(hess.diff) > 0 && det(hess.drif) > 0){
      qbic <- summary(mle)@m2logL+log(det(hess.diff))+log(det(hess.drif))
    }else{
      warning(cat("det(hessian of diffusion coefficient)<=0 or det(hessian of drift coefficient)<=0 \n"))
      qbic <- summary(mle)@m2logL+length(alpha)*log(Terminal)+length(beta)*log(n)
    }
    
    final.res <- list(par = esti, BIC = bic, QBIC = qbic, CIC = cic)
    return(final.res)
  }else{
    yuima.stop("This function is for yuima-class.")
  }
}




