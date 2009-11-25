\name{simulate}
\alias{simulate}
\title{Simulator function for multi-dimensional stochastic processes}
\description{Simulate multi-dimensional stochastic processes.}
\usage{
simulate(object, nsim, seed, xinit, true.parameter, space.discretized = FALSE, 
 increment.W = NULL, increment.L = NULL, methodfGn = "Cholesky")
}
\arguments{
  \item{object}{an \code{yuima} object.}
  \item{xinit}{initial value vector of state variables.}
  \item{true.parameter}{named list of parameters.}
  \item{space.discretized}{flag to switch to space-discretized Euler
	Maruyama method.}
  \item{increment.W}{to specify Wiener increment for each time tics in advance.}
  \item{increment.L}{to specify Levy increment for each time tics in advance.}
  \item{nsim}{Not used yet. Included only to match the standard genenirc in package \code{stats}.}
  \item{seed}{Not used yet. Included only to match the standard genenirc in package \code{stats}.}
  \item{methodfGn}{simulation methods for fractional Gaussian noise.}
}
\details{
\code{simulate} is a function to solve SDE using the Euler-Maruyama
method. This function supports usual Euler-Maruyama method for
multidimensional SDE, and space
discretized Euler-Maruyama method for one dimensional SDE.
}
\value{
  \item{yuima}{A multi-dimensional \code{yuima} object(time series object).}
}
\author{YUIMA Project Team}
\note{There may be missing information in the simulate description.
Please contribute with suggestions and fixings.
}
\examples{
set.seed(123)

# Path-simulation for 1-dim diffusion process. 
# dXt = 3*Xt*dt + dWt
mod <- setModel(drift="3*y", diffusion=1, solve.variable=c("y"))
str(mod)

# Set the model in an `yuima' object with a sampling scheme. 
T <- 1
division <- 1000
samp <- setSampling(Terminal=T, division=division)
ou <- setYuima(model=mod, sampling=samp)

# Solve SDEs using Euler-Maruyama method. 
ou <- simulate(ou, xinit=1)
plot(ou)



# Path-simulation for 1-dim diffusion process. 
# dXt = theta*Xt*dt + dWt
mod1 <- setModel(drift="theta*y", diffusion=1, solve.variable=c("y"))
str(mod1)
ou1 <- setYuima(model=mod1, sampling=samp)

# Solve SDEs using Euler-Maruyama method. 
ou1 <- simulate(ou, xinit=1, true.p = list(theta=-1))
plot(ou1)



# A multi-dimensional (correlated) diffusion process. 
# To describe the following model: 
# X=(X1,X2,X3); dXt = U(t,Xt)dt + V(t)dWt
# For drift coeffcient
U <- c("-x1","-2*x2","-t*x3")
# For diffusion coefficient of X1 
v1 <- function(t) 0.5*sqrt(t)
# For diffusion coefficient of X2
v2 <- function(t) sqrt(t)
# For diffusion coefficient of X3
v3 <- function(t) 2*sqrt(t)
# correlation
rho <- function(t) sqrt(1/2)
# coefficient matrix for diffusion term
V <- matrix( c( "v1(t)",
                "v2(t) * rho(t)",
                "v3(t) * rho(t)",
                "",
                "v2(t) * sqrt(1-rho(t)^2)",
                "",
                "",
                "",
                "v3(t) * sqrt(1-rho(t)^2)" 
               ), 3, 3)
# Model sde using "setModel" function
cor.mod <- setModel(drift = U, diffusion = V,
state.variable=c("x1","x2","x3"), 
solve.variable=c("x1","x2","x3") )
str(cor.mod)

# Set the `yuima' object. 
cor.samp <- setSampling(Terminal=T, division=division)
cor <- setYuima(model=cor.mod, sampling=cor.samp)

# Solve SDEs using Euler-Maruyama method. 
set.seed(123)
cor <- simulate(cor)
plot(cor)

# solve SDEs using Space-discretized Euler-Maruyama method
v4 <- function(t,x){
  return(0.5*(1-x)*sqrt(t))
}
mod_sd <- setModel(drift = c("0.1*x1", "0.2*x2"),
                     diffusion = c("v1(t)","v4(t,x2)"),
                     solve.var=c("x1","x2")
                     )
samp_sd <- setSampling(Terminal=T, division=division)
sd <- setYuima(model=mod_sd, sampling=samp_sd)
sd <- simulate(sd, xinit=c(1,1), space.discretized=TRUE)
plot(sd)


## example of simulation by specifying increments
## Path-simulation for 1-dim diffusion process
## dXt = 3*Xt*dt + dWt

mod <- setModel(drift="3*y", diffusion=3,solve.variable=c("y"))
str(mod)

## Set the model in an `yuima' object with a sampling scheme. 
Terminal <- 1
division <- 500
mod.sampling <- setSampling(Terminal=Terminal, division=division)
yuima.mod <- setYuima(model=mod, sampling=mod.sampling)

##use original increment
delta <- Terminal/division
my.dW <- rnorm(division * yuima.mod@model@noise.number, 0, sqrt(delta))
my.dW <- t(matrix(my.dW, nrow=division, ncol=yuima.mod@model@noise.number))

## Solve SDEs using Euler-Maruyama method.
yuima.mod <- simulate(yuima.mod,
                      xinit=1,
                      space.discretized=FALSE,
                      increment.W=my.dW)
if( !is.null(yuima.mod) ){
  x11()
  plot(yuima.mod)
}

## A multi-dimensional (correlated) diffusion process. 
## To describe the following model: 
## X=(X1,X2,X3); dXt = U(t,Xt)dt + V(t)dWt
## For drift coeffcient
U <- c("-x1","-2*x2","-t*x3")
## For process 1
diff.coef.1 <- function(t) 0.5*sqrt(t)
## For process 2
diff.coef.2 <- function(t) sqrt(t)
## For process 3
diff.coef.3 <- function(t) 2*sqrt(t)
## correlation
cor.rho <- function(t) sqrt(1/2)
## coefficient matrix for diffusion term
V <- matrix( c( "diff.coef.1(t)",
               "diff.coef.2(t) * cor.rho(t)",
               "diff.coef.3(t) * cor.rho(t)",
               "",
               "diff.coef.2(t)",
               "diff.coef.3(t) * sqrt(1-cor.rho(t)^2)",
               "diff.coef.1(t) * cor.rho(t)",
               "",
               "diff.coef.3(t)" 
               ), 3, 3)
## Model sde using "setModel" function
cor.mod <- setModel(drift = U, diffusion = V,
                    solve.variable=c("x1","x2","x3") )
str(cor.mod)
## Set the `yuima' object.
set.seed(123)
obj.sampling <- setSampling(Terminal=Terminal, division=division)
yuima.obj <- setYuima(model=cor.mod, sampling=obj.sampling)

##use original dW
my.dW <- rnorm(division * yuima.obj@model@noise.number, 0, sqrt(delta))
my.dW <- t(matrix(my.dW, nrow=division, ncol=yuima.obj@model@noise.number))

## Solve SDEs using Euler-Maruyama method.
yuima.obj.path <- simulate(yuima.obj, space.discretized=FALSE, increment.W=my.dW)
if( !is.null(yuima.obj.path) ){
  x11()
  plot(yuima.obj.path)
}


##:: sample for Levy process ("CP" type)
obj.model <- setModel(drift=c("-theta*x"), diffusion="sigma",
jump.coeff="z", measure=list(intensity="1", df=list("dnorm(z, 0, 1)")),
measure.type="CP", solve.variable="x")

##:: Parameters
lambda <- 3
theta <- 6
sigma <- 1
xinit <- runif(1)
n <- 500
h <- n^(-0.7)
eps <- h/50
division <- 50*n
T <- n*h

obj.sampling <- setSampling(Terminal=T, division=division)
obj.yuima <- setYuima(model=obj.model, sampling=obj.sampling)
X <- simulate(obj.yuima, xinit=xinit, true.parameter=list(theta=theta, sigma=sigma))
plot(X)


##:: sample for Levy process ("code" type)
## dX_{t} = -x dt + dZ_t
obj.model <- setModel(drift="-x", xinit=1, jump.coeff="1", measure.type="code", measure=list(df="rIG(z, 1, 0.1)"))
obj.sampling <- setSampling(Terminal=10, division=10000)
obj.yuima <- setYuima(model=obj.model, sampling=obj.sampling)
result <- simulate(obj.yuima)
plot(result)

##:: sample for multidimensional Levy process ("code" type)
## dX = (theta - A X)dt + dZ,
##    theta=(theta_1, theta_2) = c(1,.5)
##    A=[a_ij], a_11 = 2, a_12 = 1, a_21 = 1, a_22=2
##  xinit <- c(1,1)
##  beta <- c(.1,.1)
##  mu <- c(0,0)
##  Lambda <- matrix(c(1,0,0,1),2,2)
##  obj.model <- setModel(drift=c("1 - 2*x1-x2",".5-x1-2*x2"), xinit=c(1,1), solve.variable=c("x1","x2"), jump.coeff="1", measure.type="code", measure=list(df="rIG(z, alpha=1, beta=beta, mu=mu, Lambda=Lambda)"))
##  obj.sampling <- setSampling(Terminal=10, division=10000)
##  obj.yuima <- setYuima(model=obj.model, sampling=obj.sampling)
##  result <- simulate(obj.yuima)	   
## plot(result)

}

\keyword{ts}
