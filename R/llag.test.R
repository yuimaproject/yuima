
llag.test <- function(x, from = -Inf, to = Inf, division = FALSE, 
                      grid, R = 999, parallel = "no", 
                      ncpus = getOption("boot.ncpus", 1L), 
                      cl = NULL, tol = 1e-6){
  
  x <- get.zoo.data(x)
  
  d <- length(x)
  d.size <- d*(d-1)/2
  
  # allocate memory
  ser.times <- vector(d, mode="list") # time index in 'x'
  ser.diffX <- vector(d, mode="list") # difference of data
  vol <- double(d)
  
  # treatment of the grid (2016-07-04: we implement this before the NA treatment)
  if(missing(grid)) 
    grid <- make.grid(d, d.size, x, from, to, division)
  
  # Set the tolerance to avoid numerical erros in comparison
  #tol <- 1e-6
  
  for(i in 1:d){
    
    # NA data must be skipped
    idt <- which(is.na(x[[i]]))
    if(length(idt>0)){
      x[[i]] <- x[[i]][-idt]
    }
    if(length(x[[i]])<2) {
      stop("length of data (w/o NA) must be more than 1")
    }
    
    # set data and time index
    ser.times[[i]] <- as.numeric(time(x[[i]]))/tol
    # set difference of the data 
    ser.diffX[[i]] <- diff( as.numeric(x[[i]]) )
    vol[i] <- sum(ser.diffX[[i]]^2)
  }
  
  pval <- matrix(0,d,d)
  mcov <- diag(vol)
  mcor <- diag(d)
  
  rownames(pval) <- names(x)
  rownames(mcov) <- names(x)
  rownames(mcor) <- names(x)
  colnames(pval) <- names(x)
  colnames(mcov) <- names(x)
  colnames(mcor) <- names(x)
  
  # treatment of the grid
  #if(missing(grid)) 
  #  grid <- make.grid(d, d.size, x, from, to, division)
  
  if(is.list(grid)){
    G <- relist(unlist(grid)/tol, grid)
  }else{
    if(is.numeric(grid)){
      G <- data.frame(matrix(grid/tol,length(grid),d.size))
      grid <- data.frame(matrix(grid,length(grid),d.size))
    }else{
      print("llag:invalid grid")
      return(NULL)
    }
  }
  
  # core part
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      
      time1 <- ser.times[[i]]
      time2 <- ser.times[[j]] 
      
      num <- d*(i-1) - (i-1)*i/2 + (j-i)
      
      hry <- function(data){
        
        dx <- data$dx
        dy <- data$dy
        
        U <- .C("HYcrosscov2",
                as.integer(length(G[[num]])),
                as.integer(length(time1)),
                as.integer(length(time2)),
                as.double(G[[num]]),
                as.double(time1),
                as.double(time2),
                #double(length(time2)),
                as.double(dx),
                as.double(dy),
                value=double(length(G[[num]])),
                PACKAGE = "yuima")$value
        
        return(max(abs(U)))
      }
      
      ran.gen <- function(data, mle){
        dx <- data$dx
        dy <- data$dy
        list(dx = (2*rbinom(length(dx),1,0.5)-1) * dx, 
             dy = (2*rbinom(length(dy),1,0.5)-1) * dy)
      }
      
      out <- boot(list(dx = ser.diffX[[i]], dy = ser.diffX[[j]]),
                  hry, R, sim = "parametric", ran.gen = ran.gen,
                  parallel = parallel, ncpus = ncpus, cl = cl)
      
      pval[i,j] <- mean(out$t > out$t0)
      mcov[i,j] <- out$t0
      mcor[i,j] <- out$t0/sqrt(vol[i]*vol[j])
      
      pval[j,i] <- pval[i,j]
      mcov[j,i] <- mcov[i,j]
      mcor[j,i] <- mcor[i,j]
      
    }
  }
  
  return(list(p.values = pval, max.cov = mcov, max.corr = mcor))
}
