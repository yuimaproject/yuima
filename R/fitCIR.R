## function to provide the preliminary explicit estimator
get_preliminaryEstimators <- function(data){ # use data returned from calling simulate_CIR()
  n <- dim(data)[2]-1 # we have observations at t_j for j=0,...,n, this equals n+1 observations in total
                      # -> therefore, n is the number of observations minus 1.
  h <- as.numeric(data[1,2]) # as the vector containing the t always starts with 0, this equals diff(data[1,])[1]
                 # and thus gives h if the t_j are chosen to be equidistant.
  X_major <- data[2,-1] # the observations of the CIR process starting at t_1=h.
  X_minor <- head(data[2,], -1) # the observations of the CIR process discarding j=n.
  X_mean_major <- mean(X_major) # mean of all observations without j=0.
  X_mean_minor <- mean(X_minor) # mean of all observations without j=n.
  beta_0n <- -1/h * log( sum((X_minor - X_mean_minor) * (X_major - X_mean_major)) / # calculate estimate for beta from Eq. (2.6)
              sum( (X_minor - X_mean_minor) ^ 2 ) ) 
  alpha_0n <- (X_mean_major - exp(-beta_0n * h) * X_mean_minor) * beta_0n / (1 - exp(-beta_0n * h) ) # calculate estimate for alpha from Eq. (2.5)
  to_gamma_numerator <- (X_major - exp(-beta_0n * h) * X_minor - alpha_0n / beta_0n * (1 - exp(-beta_0n * h))) ^ 2 
  to_gamma_denominator <- (1 - exp(-beta_0n * h)) / beta_0n * (exp(-beta_0n * h) * X_minor + alpha_0n * (1 - exp(-beta_0n * h)) / (2 * beta_0n) )
  gamma_0n <- 1 / n * sum(to_gamma_numerator / to_gamma_denominator) # calculate estimate for gamma from Eq. (2.7), where first the numerator and denominator
                                                                     # are computed separately. 
  return(as.matrix(c(alpha_0n, beta_0n, gamma_0n))) # return the estimated parameter vector.
}


### ONE STEP IMPROVEMENT
# first we need a few auxiliary functions

#the mean of the CIR with parameter (alpha,beta,gamma) conditioned on X_h=x
mu<- function(alpha,beta,gamma,x,h){
  return(exp(-beta*h)*x+alpha/beta*(1-exp(-beta*h)))
}

#derivatives of the mean of CIR with parameter (alpha,beta,gamma) conditioned on X_h=x
mu_alpha <-function(alpha,beta,gamma,x,h)
{
  return ((1-exp(-beta*h))/beta)
}

mu_alpha_alpha<-function(alpha,beta,gamma,x,h)
{
  return (0)
}

mu_alpha_beta <- function(alpha,beta,gamma,x,h)
{
  return ( (exp(-beta*h)*(beta*h+1)-1)/beta^2)
}

mu_alpha_gamma <- function(alpha,beta,gamma,x,h)
{
  return (0)
}

mu_beta <- function(alpha,beta,gamma,x,h)
{
  return(-h*exp(-beta*h)*x+alpha*(exp(-beta*h)*(beta*h+1)-1)/beta^2)
}

mu_beta_beta <- function(alpha,beta,gamma,x,h)
{
  return(h^2*exp(-beta*h)*x+alpha*(exp(-beta*h)*(-beta^2*h^2-2*beta*h-2)+2)/beta^3)
}

mu_beta_gamma <- function(alpha,beta,gamma,x,h)
{
  return(0)
}

mu_gamma <- function(alpha,beta,gamma,x,h)
{
  return(0)
}

mu_gamma_gamma <- function(alpha,beta,gamma,x,h)
{
  return(0)
}

#the variance of CIR with parameter (alpha,beta,gamma) conditioned on X_h=x

sigma_mod <- function(alpha,beta,gamma,x,h){  ## HM 2022/6/21: sigma --> sigma_mod
  return(gamma/beta*(1-exp(-beta*h))*(exp(-beta*h)*x+alpha/(2*beta)*(1-exp(-beta*h))))
}

#derivatives of the variance of CIR with parameter (alpha,beta,gamma) conditioned on X_h=x

sigma_alpha <- function(alpha,beta,gamma,x,h)
{
  return(gamma*(1-exp(-beta*h))^2/(2*beta^2))
}

sigma_alpha_alpha <- function(alpha,beta,gamma,x,h)
{
  return(0)
}

sigma_alpha_beta <- function(alpha,beta,gamma,x,h)
{
  return(- gamma*exp(-2*beta*h)*(exp(beta*h)-1)*(-beta*h+exp(beta*h)-1)/beta^3)
}

sigma_alpha_gamma <- function(alpha,beta,gamma,x,h)
{
  return((1-exp(-beta*h))^2/(2*beta^2))
}

sigma_beta <- function(alpha, beta, gamma, x, h) 
{
  exp(-2 * h * beta) / beta ^3 * ( x * beta * (2 * h * beta - exp(h * beta) * (h * beta + 1) + 1) -
                                     alpha * (exp(beta * h) - 1) * (-h * beta + exp(h * beta) - 1) )
}

sigma_beta_beta <- function(alpha,beta,gamma,x,h)
{
  exp(-2 * beta * h) / (beta ^ 4) * (
    alpha * (2 * h * beta * (h * beta + 2) + 3 * exp(2 * beta * h) - exp(h * beta) * (h * beta * (h * beta + 4) + 6) + 3) +
      x * beta * (exp(beta * h) * (h ^ 2 * beta ^ 2 + 2 * h * beta + 2) - 2 * (2 * h ^ 2 * beta ^ 2 + 2 * h * beta + 1))
  )
}

sigma_beta_gamma <- function(alpha,beta,gamma,x,h)
{
  return(h*exp(-beta*h)*(exp(-beta*h)*x/beta+alpha*(1-exp(-beta*h))/(2*beta^2))+(1-exp(-beta*h))*
           (-exp(-beta*h)*(beta*h+1)/beta^2*x+alpha*exp(-beta*h)*(beta*h-2*exp(beta*h)+2)/(2*beta^3)))
}

sigma_gamma <- function(alpha,beta,gamma,x,h)
{
  return((1-exp(-beta*h))*(exp(-beta*h)*x+alpha/(2*beta)*(1-exp(-beta*h)))/beta)
}

sigma_gamma_gamma <- function(alpha,beta,gamma,x,h)
{
  return(0)
}

#the inverse of the Hessian matrix of H_n
H_n_Hessian <- function(alpha,beta,gamma,x,h)
{
  mu <- mu(alpha,beta,gamma,x,h)
  
  mu_1 <- mu_alpha(alpha,beta,gamma,x,h)
  mu_2 <- mu_beta(alpha,beta,gamma,x,h)
  mu_3 <- mu_gamma(alpha,beta,gamma,x,h)
  
  mu_11 <- mu_alpha_alpha(alpha,beta,gamma,x,h)
  mu_12 <- mu_alpha_beta(alpha,beta,gamma,x,h)
  mu_13 <- mu_alpha_gamma(alpha,beta,gamma,x,h)
  mu_22 <- mu_beta_beta(alpha,beta,gamma,x,h)
  mu_23 <- mu_beta_gamma(alpha,beta,gamma,x,h)
  mu_33 <- mu_gamma_gamma(alpha,beta,gamma,x,h)
  
  sigma <- sigma_mod(alpha,beta,gamma,x,h)  ## HM 2022/6/21: sigma --> sigma_mod
  
  sigma_1 <- sigma_alpha(alpha,beta,gamma,x,h)
  sigma_2 <- sigma_beta(alpha,beta,gamma,x,h)
  sigma_3 <- sigma_gamma(alpha,beta,gamma,x,h)
  
  sigma_11 <- sigma_alpha_alpha(alpha,beta,gamma,x,h)
  sigma_12 <- sigma_alpha_beta(alpha,beta,gamma,x,h)
  sigma_13 <- sigma_alpha_gamma(alpha,beta,gamma,x,h)
  sigma_22 <- sigma_beta_beta(alpha,beta,gamma,x,h)
  sigma_23 <- sigma_beta_gamma(alpha,beta,gamma,x,h)
  sigma_33 <- sigma_gamma_gamma(alpha,beta,gamma,x,h)
  
  H_11 <- -0.5*sum( (sigma_11*sigma-sigma_1*sigma_1)/sigma^2-(sigma_11*sigma-2*sigma_1*sigma_1)/sigma^3*(x-mu)^2+
                      2*sigma_1/sigma^2*(x-mu)*mu_1-
                      2*(mu_11*sigma-mu_1*sigma_1)*(x-mu)/sigma^2+2*mu_1*mu_1/sigma)
  
  
  H_22 <- -0.5*sum( (sigma_22*sigma-sigma_2*sigma_2)/sigma^2-(sigma_22*sigma-2*sigma_2*sigma_2)/sigma^3*(x-mu)^2+
                      2*sigma_2/sigma^2*(x-mu)*mu_2-
                      2*(mu_22*sigma-mu_2*sigma_2)*(x-mu)/sigma^2+2*mu_2*mu_2/sigma)
  
  H_33 <- -0.5*sum( (sigma_33*sigma-sigma_3*sigma_3)/sigma^2-(sigma_33*sigma-2*sigma_3*sigma_3)/sigma^3*(x-mu)^2+
                      2*sigma_3/sigma^2*(x-mu)*mu_3-
                      2*(mu_33*sigma-mu_3*sigma_3)*(x-mu)/sigma^2+2*mu_3*mu_3/sigma)
  
  
  H_12 <- -0.5*sum( (sigma_12*sigma-sigma_1*sigma_2)/sigma^2-(sigma_12*sigma-2*sigma_1*sigma_2)/sigma^2*(x-mu)^2+
                      2*sigma_1/sigma^2*(x-mu)*mu_2-
                      2*(mu_12*sigma-mu_1*sigma_2)*(x-mu)/sigma^2+2*mu_1*mu_2/sigma)
  
  H_13 <- -0.5*sum( (sigma_13*sigma-sigma_1*sigma_3)/sigma^2-(sigma_13*sigma-2*sigma_1*sigma_3)/sigma^3*(x-mu)^2+
                      2*sigma_1/sigma^2*(x-mu)*mu_3-
                      2*(mu_13*sigma-mu_1*sigma_3)*(x-mu)/sigma^2+2*mu_1*mu_3/sigma)
  
  H_23 <- -0.5*sum( (sigma_23*sigma-sigma_2*sigma_3)/sigma^2-(sigma_23*sigma-2*sigma_2*sigma_3)/sigma^3*(x-mu)^2+
                      2*sigma_2/sigma^2*(x-mu)*mu_3-
                      2*(mu_23*sigma-mu_2*sigma_3)*(x-mu)/sigma^2+2*mu_2*mu_3/sigma)
  
  
  H_det <- H_11*H_22*H_33+H_12*H_23*H_13+H_13*H_12*H_23-H_13*H_22*H_13-H_11*H_23*H_23-H_12*H_12*H_33
  
  return(matrix(1/H_det*c(H_22*H_33-H_23^2,H_13*H_23-H_12*H_33, H_12*H_23-H_13*H_22,
                          H_13*H_23-H_12*H_33,H_11*H_33-H_13^2,H_13*H_12-H_11*H_23,
                          H_12*H_23-H_13*H_22,H_13*H_12-H_11*H_23, H_11*H_22-H_12^2),nrow=3))
  
}

#derivatives of the function H_n(theta) 

H_n_alpha <- function(alpha,beta,gamma,x,h){
  x_minor <- head(x, -1)
  x_major <- x[-1]
  sigma_deriv <- sigma_alpha(alpha,beta,gamma,x_minor,h)
  sigma <- sigma_mod(alpha,beta,gamma,x_minor,h)  ## HM 2022/6/21: sigma --> sigma_mod
  mu <- mu(alpha,beta,gamma,x_minor,h)
  mu_deriv <- mu_alpha(alpha,beta,gamma,x_minor,h)
  return (sum(-0.5*(sigma_deriv/sigma-sigma_deriv/sigma^2*(x_major-mu)^2-2/sigma*(x_major-mu)*mu_deriv)))
}

H_n_beta <- function(alpha,beta,gamma,x,h){
  x_minor <- head(x, -1)
  x_major <- x[-1]
  sigma_deriv <- sigma_beta(alpha,beta,gamma,x_minor,h)
  sigma <- sigma_mod(alpha,beta,gamma,x_minor,h)  ## HM 2022/6/21: sigma --> sigma_mod
  mu <- mu(alpha,beta,gamma,x_minor,h)
  mu_deriv <- mu_beta(alpha,beta,gamma,x_minor,h)
  return (sum(-0.5*(sigma_deriv/sigma-sigma_deriv/sigma^2*(x_major-mu)^2-2/sigma*(x_major-mu)*mu_deriv)))
}

H_n_gamma <- function(alpha,beta,gamma,x,h){
  x_minor <- head(x, -1)
  x_major <- x[-1]
  sigma_deriv <- sigma_gamma(alpha,beta,gamma,x_minor,h)
  sigma <- sigma_mod(alpha,beta,gamma,x_minor,h)  ## HM 2022/6/21: sigma --> sigma_mod
  mu <- mu(alpha,beta,gamma,x_minor,h)
  return (sum(-0.5*(sigma_deriv/sigma-sigma_deriv/sigma^2*(x_major-mu)^2)))
}

get_finalEstimators_scoring <- function(data) {
  prelim_estim <- get_preliminaryEstimators(data) #We need the preliminary estimator for the one step improvement
  alpha <- prelim_estim[1]
  beta <- prelim_estim[2]
  gamma <- prelim_estim[3]
  T <- tail(data[1,],1) #T is the time of the last observation n*h
  n<- length(data[1,])-1 #number of observations minus 1
  x<-data[2,-1]
  h<-data[1,2]
  
  H_n<-c(H_n_alpha(alpha,beta,gamma,x,h),H_n_beta(alpha,beta,gamma,x,h),H_n_gamma(alpha,beta,gamma,x,h)) #calculate estimate from Eq. 
  return(c(alpha,beta,gamma)+diag(c(1/T,1/T,1/n))%*%
           matrix(c(alpha*(2*alpha-gamma)/beta, 2*alpha-gamma,0,2*alpha-gamma,2*beta,0,0,0,2*gamma^2),nrow=3,byrow=TRUE)%*%
           H_n) 
}

get_finalEstimators_newton <- function(data) {
  prelim_estim <- get_preliminaryEstimators(data)
  alpha <- prelim_estim[1]
  beta <- prelim_estim[2]
  gamma <- prelim_estim[3]
  T <- tail(data[1,],1)
  n<- length(data[1,])-1
  x<-data[2,-1]
  h<-data[1,2]
  
  H_n<-c(H_n_alpha(alpha,beta,gamma,x,h),H_n_beta(alpha,beta,gamma,x,h),H_n_gamma(alpha,beta,gamma,x,h))
  return(c(alpha,beta,gamma)-H_n_Hessian(alpha,beta,gamma,x,h)%*%H_n)
}

fitCIR<- function(data){
  if( (max(diff(data[1,]))/ min(diff(data[1,])) -1) > 1e-10 ) stop('Please use equidistant sampling points')
  # fittedCIR.yuima <- setModel(drift="alpha - beta*x", diffusion="sqrt(gamma*x)", solve.variable=c("x"))
  return(list(get_preliminaryEstimators(data),get_finalEstimators_newton(data),get_finalEstimators_scoring(data)
              # ,fittedCIR.yuima
              ))
}


# ###Examples of usage
# ##If the sampling points are not equidistant, there will be an corresponding output. 
# data <- sim.CIR(alpha=3,beta=1,gamma=1,time.points = c(0,0.1,0.2,0.25,0.3))
# fit.CIR(data)
# ##Otherwise it calculate the three estimators
# data <- sim.CIR(alpha=3,beta=1,gamma=1,n=1000,h=0.1,equi.dist=TRUE)
# fit.CIR(data)
