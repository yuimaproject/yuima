# Scale-by-scale lead-lag estimation by wavelets

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
  
  acw <- wavethresh::PsiJ(-J, filter.number = N, family = "DaubExPhase")
  #acw <- wavethresh::PsiJ(-J, filter.number = N, family = family)
  
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
