## Simulate Cox-Ingersoll-Ross process with parameters alpha, beta and gamma at times specified via time.points
simCIR <- function (time.points, n, h, alpha, beta, gamma, equi.dist=FALSE ) {
  
  # generate an equidistant time vector of length n+1 and distant h between observations 
  if (equi.dist==TRUE) {time.points <- 0:n*h }
  
  # must start in t=0, otherwise t_vec is adjusted
  if ( time.points[1] != 0 ) { time.points <- c(0, time.points) }
  
  # define auxiliary variables, following notation of Malham and Wiese
  nu <- 4 * alpha /  gamma   # degrees of freedom
  eta_vec <- 4 * beta * exp(-beta * diff(time.points) ) /   # auxiliary vector for the computation of the
              (gamma  * (1 - exp(-beta * diff(time.points) )) ) # non-centrality parameter in each step
  
  # sample X_0 from stationary distribution
  X <- rgamma(1, scale = gamma / (2 * beta), shape = 2 * alpha / gamma)
  
  # compute X_t iteratively, using Prop. 1 from Malham and Wiese (2012)
  for ( i in seq_along(eta_vec) ) {
    lambda <- tail(X, 1) * eta_vec[i] # non-centrality parameter of the conditional distribution
    X <- c(X, rchisq(1, df = nu, ncp = lambda) * exp(-beta * diff(time.points)[i]) / eta_vec[i]) # calculate
                                                    # next step of the CIR
  }
  
  # return data
  return(rbind(t = time.points, X = X)) # first row: time points, second row: CIR at time point
}

