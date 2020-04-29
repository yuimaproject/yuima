
# sinc function
sinc <- function(x){
  
  out <- rep(1, length(x))
  
  idx <- abs(x) > 1e-7
  out[idx] <- sin(x[idx])/x[idx] 
  
  return(out)
}

# Littlewood-Paley wavelet function
psi.lp <- function(x) 
  2 * sinc(2 * pi * x) - sinc(pi * x)

# coefficient matrices for circulant embedding
simBmllag.coef <- function(n, J, rho, theta, delta = 1/2^(J+1)){
  
  if(length(rho) < J + 1) 
    rho <- append(rho, double(J + 1 - length(rho)))
  
  if(length(theta) < J + 1) 
    theta <- append(theta, double(J + 1 - length(theta)))
  
  m <- 3^ceiling(log(2 * n - 2, base = 3))
  M <- as.integer((m - 1)/2)
  
  tl <- ((-M):M) * delta
  
  c12 <- double(m)
  for(j in 1:(J + 1)){
    c12 <- c12 + 2^(J - j + 1) * delta^2 * rho[j] * psi.lp(2^(J - j + 1) * (tl - theta[j]))
  }
  
  c12 <- c(c12[-(1:M)], c12[1:M])
  
  A12 <- fft(c12)
  
  s <- sqrt(delta^2 - Mod(A12)^2)
  t <- sqrt(2 * (delta + s))
  
  a <- array(0, dim = c(m, 2, 2))
  
  a[ ,1,1] <- (delta + s)/t
  a[ ,2,2] <- a[ ,1,1]
  a[ ,1,2] <- A12/t
  a[ ,2,1] <- Conj(a[ ,1,2])
  
  return(a)
}

# main function
simBmllag <- function(n, J, rho, theta, delta = 1/2^(J+1),
                      imaginary = FALSE){
  
  a <- simBmllag.coef(n, J, rho, theta, delta)
  m <- dim(a)[1]
  
  z1 <- rnorm(m) + 1i * rnorm(m)
  z2 <- rnorm(m) + 1i * rnorm(m)
  
  y1 <- a[ ,1,1] * z1 + a[ ,1,2] * z2
  y2 <- a[ ,2,1] * z1 + a[ ,2,2] * z2
  
  dW <- mvfft(cbind(y1, y2))[1:n, ]/sqrt(m)
  
  if(imaginary == TRUE){
    return(dW)
  }else{
    return(Re(dW))
  }
  
}
