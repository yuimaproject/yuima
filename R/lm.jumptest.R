
# Lee-Mykland's jump test

lm.jumptest <- function(yuima, K){
  
  data <- get.zoo.data(yuima)
  d.size <- length(data)
  n <- length(yuima) - 1
  
  if(missing(K)) K <- pmin(floor(sqrt(252*n)), n)
  
  K <- rep(K, d.size)
  
  cn <- sqrt(2*log(n)) - 
    (log(pi) + log(log(n)))/(2*sqrt(2*log(n)))
  sn <- 1/sqrt(2*log(n))
   
  result <- vector(d.size,mode="list") 
  
  for(d in 1:d.size){
    
    x <- data[[d]]
    adx <- abs(diff(as.numeric(x)))
    
    v.hat <- (pi/2) * rollmeanr(adx[-n[d]] * adx[-1], 
                                k = K[d] - 2, fill = "e")
    v.hat <- append(v.hat, v.hat[1], after = 0)
    
    stat <- (max(adx/sqrt(v.hat)) - cn[d])/sn[d]
    p <- 1 - exp(-exp(-stat))
    
    result[[d]] <- list(statistic = c(LM = stat),
                        p.value = p,
                        method = "Lee and Mykland jump test",
                        data.names = paste("x",d,sep = ""))
    class(result[[d]]) <- "htest"
    
  }
  
  return(result)
}
