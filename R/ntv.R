
###################################################### 
# Volatility estimation and jump test using
# nearest neighbor truncation
######################################################

######## Volatility estimation ########

## MinRV
minrv <- function(yuima){
  
  x <- get.zoo.data(yuima)
  d <- length(x)
  out <- double(d)
  
  for(j in 1:d){
    dx2 <- diff(as.numeric(x[[j]]))^2
    out[j] <- length(dx2) * mean(pmin(dx2[-length(dx2)], dx2[-1]))
  }
  
  return((pi/(pi - 2)) * out)
}

## MedRV
medrv <- function(yuima){
  
  x <- get.zoo.data(yuima)
  d <- length(x)
  out <- double(d)
  
  for(j in 1:d){
    dx2 <- diff(as.numeric(x[[j]]))^2
    out[j] <- length(dx2) * mean(rollmedian(dx2, k = 3))
  }
  
  return((pi/(6 - 4*sqrt(3) + pi)) * out)
}


######## Jump test ########

## function to compute log-type test statistics
logstat <- function(rv, iv, iq, n, theta, adj){
  
  avar <- iq/iv^2
  
  if(adj) avar[avar < 1] <- 1
  
  return(log(rv/iv)/sqrt(theta * avar/n))
}

## function to compute ratio-type test statistics
ratiostat <- function(rv, iv, iq, n, theta, adj){
  
  avar <- iq/iv^2
  
  if(adj) avar[avar < 1] <- 1
  
  return((1 - iv/rv)/sqrt(theta * avar/n))
}

## Jump test based on minRV
minrv.test <- function(yuima, type = "ratio", adj = TRUE){
  
  data <- get.zoo.data(yuima)
  d.size <- length(data)
  
  rv <- double(d.size)
  iv <- double(d.size)
  iq <- double(d.size)
  n <- integer(d.size)
  
  for(d in 1:d.size){
    
    dx2 <- diff(as.numeric(data[[d]]))^2
    n[d] <- length(dx2)
    
    rv[d] <- sum(dx2)
    
    obj <- pmin(dx2[-n[d]], dx2[-1])
    iv[d] <- (pi/(pi - 2)) * n[d] * mean(obj)
    iq[d] <- (pi/(3*pi - 8)) * n[d]^2 * mean(obj^2)
    
  }
  
  stat <- switch(type,
                 "standard" = (rv - iv)/sqrt(1.81 * iq/n),
                 "ratio" = ratiostat(rv, iv, iq, n, 1.81, adj),
                 "log" = logstat(rv, iv, iq, n, 1.81, adj))
  
  p <- pnorm(stat, lower.tail = FALSE)
  
  result <- vector(d.size,mode="list") 
  
  for(d in 1:d.size){
    result[[d]] <- list(statistic = c(ADS = stat[d]),
                        p.value = p[d],
                        method = "Andersen-Dobrev-Schaumburg jump test based on MinRV",
                        data.names = paste("x",d,sep = ""))
    class(result[[d]]) <- "htest"
  }
  
  return(result)
}

## Jump test based on medRV
medrv.test <- function(yuima, type = "ratio", adj = TRUE){
  
  data <- get.zoo.data(yuima)
  d.size <- length(data)
  
  rv <- double(d.size)
  iv <- double(d.size)
  iq <- double(d.size)
  n <- integer(d.size)
  
  for(d in 1:d.size){
    
    dx2 <- diff(as.numeric(data[[d]]))^2
    n[d] <- length(dx2)
    
    rv[d] <- sum(dx2)
    
    obj <- rollmedian(dx2, k = 3)
    iv[d] <- (pi/(6 - 4*sqrt(3) + pi)) * n[d] * mean(obj)
    iq[d] <- (3*pi/(9*pi + 72 - 52*sqrt(3))) * n[d]^2 * mean(obj^2)
    
  }
  
  stat <- switch(type,
                 "standard" = (rv - iv)/sqrt(0.96 * iq/n),
                 "ratio" = ratiostat(rv, iv, iq, n, 0.96, adj),
                 "log" = logstat(rv, iv, iq, n, 0.96, adj))
  
  p <- pnorm(stat, lower.tail = FALSE)
  
  result <- vector(d.size,mode="list") 
  
  for(d in 1:d.size){
    result[[d]] <- list(statistic = c(ADS = stat[d]),
                        p.value = p[d],
                        method = "Andersen-Dobrev-Schaumburg jump test based on MedRV",
                        data.names = paste("x",d,sep = ""))
    class(result[[d]]) <- "htest"
  }
  
  return(result)
}

