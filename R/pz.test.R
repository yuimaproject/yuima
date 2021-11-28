
# Podolskij-Ziggel's Jump Test

pz.test <- function(yuima, p = 4, threshold = "local", tau = 0.05){
  
  data <- get.zoo.data(yuima)
  d.size <- length(data)
  
  adx <- vector(d.size,mode="list")
  
  for(d in 1:d.size){
    
    X <- as.numeric(data[[d]])
    idt <- which(is.na(X))
    if(length(idt>0)){
      X <- X[-idt]
    }
    if(length(X)<2) {
      stop("length of data (w/o NA) must be more than 1")
    }
    adx[[d]] <- abs(diff(X))
  }
  
  n <- sapply(adx, "length")
  
  if(threshold == "local"){
    thr <- local_univ_threshold(data)
  }else if(threshold == "PZ"){
    thr <- 2.3^2 * mpv(yuima, r = c(1, 1)) * n^(-0.8)
  }else{
    thr <- threshold
  }
  
  # core part
  result <- vector(d.size, mode = "list")
  
  for(d in 1:d.size){
    # thresholding
    tr.idx <- (adx[[d]]^2 <= thr[[d]])
    
    eta <- sample(c(1 - tau, 1 + tau), size = n[d], replace = TRUE)
    
    numer <- sum(adx[[d]]^p * (1 - eta * tr.idx))
    denom <- tau * sqrt(sum(adx[[d]]^(2*p) * tr.idx))
    
    PZ <- numer/denom
    
    pval <- pnorm(PZ,lower.tail=FALSE)
    result[[d]] <- list(statistic=c(PZ=PZ),p.value=pval,
                        method="Podolskij-Ziggel jump test",
                        data.names=paste("x",d,sep=""))
    class(result[[d]]) <- "htest"
    
  }
  
  return(result)
}
