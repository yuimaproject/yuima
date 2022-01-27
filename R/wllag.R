# Scale-by-scale lead-lag estimation by wavelets

## function to compute Daubechies' extremal phase wavelet filter
## The function is based on implementation of wavethresh package
daubechies.wavelet <- function(N){
  
  switch (N,
    "1" = {
      H <- rep(0, 2)
      H[1] <- 1/sqrt(2)
      H[2] <- H[1]
    },
    "2" = {
      H <- rep(0, 4)
      H[1] <- 0.482962913145
      H[2] <- 0.836516303738
      H[3] <- 0.224143868042
      H[4] <- -0.129409522551
    },
    "3" = {
      H <- rep(0, 6)
      H[1] <- 0.33267055295
      H[2] <- 0.806891509311
      H[3] <- 0.459877502118
      H[4] <- -0.13501102001
      H[5] <- -0.085441273882
      H[6] <- 0.035226291882
    },
    "4" = {
      H <- rep(0, 8)
      H[1] <- 0.230377813309
      H[2] <- 0.714846570553
      H[3] <- 0.63088076793
      H[4] <- -0.027983769417
      H[5] <- -0.187034811719
      H[6] <- 0.030841381836
      H[7] <- 0.032883011667
      H[8] <- -0.010597401785
    },
    "5" = {
      H <- rep(0, 10)
      H[1] <- 0.160102397974
      H[2] <- 0.603829269797
      H[3] <- 0.724308528438
      H[4] <- 0.138428145901
      H[5] <- -0.242294887066
      H[6] <- -0.032244869585
      H[7] <- 0.07757149384
      H[8] <- -0.006241490213
      H[9] <- -0.012580752
      H[10] <- 0.003335725285
    },
    "6" = {
      H <- rep(0, 12)
      H[1] <- 0.11154074335
      H[2] <- 0.494623890398
      H[3] <- 0.751133908021
      H[4] <- 0.315250351709
      H[5] <- -0.226264693965
      H[6] <- -0.129766867567
      H[7] <- 0.097501605587
      H[8] <- 0.02752286553
      H[9] <- -0.031582039318
      H[10] <- 0.000553842201
      H[11] <- 0.004777257511
      H[12] <- -0.001077301085
    },
    "7" = {
      H <- rep(0, 14)
      H[1] <- 0.077852054085
      H[2] <- 0.396539319482
      H[3] <- 0.729132090846
      H[4] <- 0.469782287405
      H[5] <- -0.143906003929
      H[6] <- -0.224036184994
      H[7] <- 0.071309219267
      H[8] <- 0.080612609151
      H[9] <- -0.038029936935
      H[10] <- -0.016574541631
      H[11] <- 0.012550998556
      H[12] <- 0.000429577973
      H[13] <- -0.001801640704
      H[14] <- 0.0003537138
    },
    "8" = {
      H <- rep(0, 16)
      H[1] <- 0.054415842243
      H[2] <- 0.312871590914
      H[3] <- 0.675630736297
      H[4] <- 0.585354683654
      H[5] <- -0.015829105256
      H[6] <- -0.284015542962
      H[7] <- 0.000472484574
      H[8] <- 0.12874742662
      H[9] <- -0.017369301002
      H[10] <- -0.044088253931
      H[11] <- 0.013981027917
      H[12] <- 0.008746094047
      H[13] <- -0.004870352993
      H[14] <- -0.000391740373
      H[15] <- 0.000675449406
      H[16] <- -0.000117476784
    },
    "9" = {
      H <- rep(0, 18)
      H[1] <- 0.038077947364
      H[2] <- 0.243834674613
      H[3] <- 0.60482312369
      H[4] <- 0.657288078051
      H[5] <- 0.133197385825
      H[6] <- -0.293273783279
      H[7] <- -0.096840783223
      H[8] <- 0.148540749338
      H[9] <- 0.030725681479
      H[10] <- -0.067632829061
      H[11] <- 0.000250947115
      H[12] <- 0.022361662124
      H[13] <- -0.004723204758
      H[14] <- -0.004281503682
      H[15] <- 0.001847646883
      H[16] <- 0.000230385764
      H[17] <- -0.000251963189
      H[18] <- 3.934732e-05
    },
    "10" = {
      H <- rep(0, 20)
      H[1] <- 0.026670057901
      H[2] <- 0.188176800078
      H[3] <- 0.527201188932
      H[4] <- 0.688459039454
      H[5] <- 0.281172343661
      H[6] <- -0.249846424327
      H[7] <- -0.195946274377
      H[8] <- 0.127369340336
      H[9] <- 0.093057364604
      H[10] <- -0.071394147166
      H[11] <- -0.029457536822
      H[12] <- 0.033212674059
      H[13] <- 0.003606553567
      H[14] <- -0.010733175483
      H[15] <- 0.001395351747
      H[16] <- 0.001992405295
      H[17] <- -0.000685856695
      H[18] <- -0.000116466855
      H[19] <- 9.358867e-05
      H[20] <- -1.3264203e-05
    }
  )
  
  return(H)
}

## function to compute autocorrelation wavelets
autocorrelation.wavelet <- function(J, N){
  
  h <- daubechies.wavelet(N)
  Phi1 <- convolve(h, h, type = "o")
  
  out <- vector(mode = "list", J)
  out[[1]] <- (-1)^((-2*N+1):(2*N-1)) * Phi1
  
  if(J > 1){
    
    Phi1 <- Phi1[seq(1,4*N-1,by=2)]
    
    for(j in 2:J){
      
      Lj <- (2^j - 1) * (2 * N - 1) + 1
      out[[j]] <- double(2*Lj - 1)
      out[[j]][seq(2*N,2*Lj-2*N,by=2)] <- out[[j-1]]
      out[[j]][seq(1,2*Lj-1,by=2)] <- 
        convolve(out[[j-1]], Phi1, type = "o")
      
    }
  }
  
  return(out)
}


## main function
wllag <- function(x, y, J = 8, N = 10, #family = "DaubExPhase", 
                  tau = 1e-3, from = -to, to = 100, 
                  verbose = FALSE, in.tau = FALSE, tol = 1e-6){
  
  time1 <- as.numeric(time(x))
  time2 <- as.numeric(time(y))
  
  grid <- seq(from, to, by = 1) * tau
  
  Lj <- (2^J - 1) * (2 * N - 1) + 1
  #Lj <- (2^J - 1) * (length(wavethresh::filter.select(N, family)$H) - 1) + 1
  grid2 <- seq(from - Lj + 1, to + Lj - 1, by = 1) * tau
  
  dx <- diff(as.numeric(x))
  dy <- diff(as.numeric(y))
  
  tmp <- .C("HYcrosscov2",
            as.integer(length(grid2)),
            as.integer(length(time2)),
            as.integer(length(time1)),
            as.double(grid2/tol),
            as.double(time2/tol),
            as.double(time1/tol),
            as.double(dy),
            as.double(dx),
            value=double(length(grid2)),
            PACKAGE = "yuima")$value
  
  #if(missing(J)) J <- floor(log2(length(grid)))
  
  #if(J < 2) stop("J must be larger than 1")
  
  #acw <- wavethresh::PsiJ(-J, filter.number = N, family = "DaubExPhase")
  #acw <- wavethresh::PsiJ(-J, filter.number = N, family = family)
  acw <- autocorrelation.wavelet(J, N)
  
  theta <- double(J)
  #covar <- double(J)
  #LLR <- double(J)
  corr <- double(J)
  crosscor <- vector("list", J)
  
  for(j in 1:J){
    
    wcov <- try(stats::filter(tmp, filter = acw[[j]], method = "c", 
                              sides = 2)[Lj:(length(grid) + Lj - 1)],
                silent = TRUE)
    #Mj <- (2^J - 2^j) * (2 * N - 1)
    #wcov <- try(convolve(tmp, acw[[j]], conj = FALSE, type = "filter")[(Mj + 1):(length(grid) + Mj)],
    #            silent = TRUE)
    
    if(class(wcov) == "try-error"){
      
      theta[j] <- NA
      crosscor[[j]] <- NA
      corr[j] <- NA
      
    }else{
      
      #tmp.grid <- grid[-attr(wcov, "na.action")]
      crosscor[[j]] <- zoo(wcov, grid)
      
      obj <- abs(wcov)
      idx1 <- which(obj == max(obj, na.rm = TRUE))
      idx <- idx1[which.max(abs(grid[idx1]))]
      # if there are multiple peaks, take the lag farthest from zero
      theta[j] <- grid[idx]
      corr[j] <- crosscor[[j]][idx]
      
    }
    
  }
  
  if(verbose == TRUE){
    
    #obj0 <- tmp[(Lj + 1):(length(grid) + Lj)]
    obj0 <- tmp[Lj:(length(grid) + Lj - 1)]/sqrt(sum(dx^2)*sum(dy^2))
    obj <- abs(obj0)
    idx1 <- which(obj == max(obj, na.rm = TRUE))
    idx <- idx1[which.max(abs(grid[idx1]))]
    # if there are multiple peaks, the lag farthest from zero
    theta.hy <- grid[idx]
    corr.hy <- obj0[idx]
    
    if(in.tau == TRUE){
      theta <- round(theta/tau)
      theta.hy <- round(theta.hy/tau)
    }
    
    result <- list(lagtheta = theta, obj.values = corr, 
                   obj.fun = crosscor, theta.hry = theta.hy, 
                   cor.hry = corr.hy, ccor.hry = zoo(obj0, grid))
    
    class(result) <- "yuima.wllag"
    
  }else{
    if(in.tau == TRUE){
      result <- round(theta/tau)
    }else{
      result <- theta
    }
  }
  
  return(result)
}

# print method for yuima.wllag-class
print.yuima.wllag <- function(x, ...){
  
  cat("Estimated scale-by-scale lead-lag parameters\n")
  print(x$lagtheta, ...)
  cat("Corresponding values of objective functions\n")
  print(x$obj.values, ...)
  cat("Estimated lead-lag parameter in the HRY sense\n")
  print(x$theta.hry, ...)
  cat("Corresponding correlation coefficient\n")
  print(x$cor.hry, ...)
  
}

# plot method for yuima.wllag class
plot.yuima.wllag <- function(x, selectJ = NULL, xlab = expression(theta), 
                             ylab = "", ...){
  
  J <- length(x$lagtheta)
  
  if(is.null(selectJ)) selectJ <- 1:J
  
  for(j in selectJ){
    #plot(x$ccor[[j]], main=paste("j=",j), xlab=expression(theta),
    #     ylab=expression(U[j](theta)), type = type, pch = pch, ...)
    plot(x$obj.fun[[j]], main = paste("j = ", j, sep =""), 
         xlab = xlab, ylab = ylab, ...)
    abline(0, 0, lty = "dotted")
  }
  
}
