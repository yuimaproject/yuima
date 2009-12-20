\name{ql}
\alias{ql}
\alias{ml.ql}
%\alias{ml.ql2}
\alias{rql}
%\alias{ql,ANY-method} 
%\alias{ml.ql,ANY-method} 
%\alias{ml.ql2,ANY-method} 
%\alias{rql,ANY-method} 
\title{Calculate quasi-likelihood and ML estimator}
\description{Calculate the quasi-likelihood and estimate of the parameters of the
  stochastic differential equation by the maximum likelihood method.}
\usage{
ml.ql(yuima,theta2,theta1,h,theta2.lim=matrix(c(0,1),1,2),theta1.lim=matrix(c(0,1),1,2),print=FALSE,BFGS=FALSE,param,interval)
ql(yuima,theta2,theta1,h,print=FALSE,param)
rql(yuima,theta2,theta1,ptheta2,ptheta1,h,print=FALSE,param,prevparam)
}
\arguments{
  \item{yuima}{a yuima object.}
  \item{theta2,theta1}{parameters of the sdeModel.}
  \item{h}{ time span of observations.}
  \item{theta2.lim, theta1.lim}{matrixes to specify the domains of the
    parameters. Vector can be available only if theta is a scalar.}
  \item{ptheta2,ptheta1}{}
  \item{print}{you can see a progress of the estimation when print is TRUE.}
  \item{BFGS}{}
  \item{param}{}
  \item{interval}{}
  \item{prevparam}{}
}
\details{
  A function ql calculate the quasi-likelihood of a time series data X with any
  parameters. A function ml.pl estimates parameters of the sdeModel by
  maximizing the quasi-likelihood.
}
\value{
  \item{QL}{a real value.}
  \item{opt}{a list with components the same as 'optim' function.}
}
\author{The YUIMA Project Team}
\note{
The function ml.ql uses the function optim internally.
}
\examples{
#dXt^e = -theta2 * Xt^e * dt + theta1 * dWt
diff.matrix <- matrix(c("theta1"), 1, 1)
ymodel <- setModel(drift=c("(-1)*theta2*x"), diffusion=diff.matrix, time.variable="t", state.variable="x", solve.variable="x")
n <- 100

ysamp <- setSampling(Terminal=(n)^(1/3), n=n) 
yuima <- setYuima(model=ymodel, sampling=ysamp)
set.seed(123)
yuima <- simulate(yuima, xinit=1, true.parameter=c(0.3, 0.1))
QL <- ql(yuima, 0.8, 0.7, h=1/((n)^(2/3)))
QL

## another way of parameter specification
##param <- list(theta2=0.8, theta1=0.7)
##QL <- ql(yuima, h=1/((n)^(2/3)), param=param)
##QL

system.time(
opt <- ml.ql(yuima, 0.8, 0.7, h=1/((n)^(2/3)), c(0, 1), c(0, 1))
)
print("True param")
print("theta2 = .3, theta1 = .1")
print("ML estimator")
opt$par

## another way of parameter specification
##interval <- list(theta2.lim=c(0,1), theta1.lim=c(0,1))
##system.time(
##opt <- ml.ql(yuima, h=1/((n)^(2/3)), param=param, interval=interval)
##)
##print("True param")
##print("theta2 = .3, theta1 = .1")
##print("ML estimator")
##opt$par

system.time(
opt <- ml.ql(yuima, 0.8, 0.7, h=1/((n)^(2/3)), c(0, 1), c(0, 1), BFGS=TRUE)
)
print("True param")
print("theta2 = .3, theta1 = .1")
print("ML estimator")
opt$par

##multidimension case
##dXt^e = - drift.matrix * Xt^e * dt + diff.matrix * dWt
diff.matrix <- matrix(c("theta1.1","theta1.2", "1", "1"), 2, 2)

drift.c <- c("-theta2.1*x1", "-theta2.2*x2", "-theta2.2", "-theta2.1")
drift.matrix <- matrix(drift.c, 2, 2)

ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, time.variable="t",
                   state.variable=c("x1", "x2"), solve.variable=c("x1", "x2"))
n <- 100
ysamp <- setSampling(Terminal=(n)^(1/3), n=n)
yuima <- setYuima(model=ymodel, sampling=ysamp)
set.seed(123)

##xinit=c(x1, x2) #true.parameter=c(theta2.1, theta2.2, theta1.1, theta1.2)
yuima <- simulate(yuima, xinit=c(1, 1), true.parameter=c(0.5, 0.3, 0.6, 0.2))

theta2 <- c(0.8, 0.2) #c(theta2.1, theta2.2)
theta1 <- c(0.7, 0.1) #c(theta1.1, theta1.2)
QL <- ql(yuima, theta2, theta1, h=1/((n)^(2/3)))
QL

## another way of parameter specification
#param <- list(theta2=theta2, theta1=theta1)
#QL <- ql(yuima, h=1/((n)^(2/3)), param=param)
#QL

theta2.1.lim <- c(0, 1)
theta2.2.lim <- c(0, 1)
theta1.1.lim <- c(0, 1)
theta1.2.lim <- c(0, 1)
theta2.lim <- t( matrix( c(theta2.1.lim, theta2.2.lim), 2, 2) )
theta1.lim <- t( matrix( c(theta1.1.lim, theta1.2.lim), 2, 2) )

system.time(
opt <- ml.ql(yuima, theta2, theta1, h=1/((n)^(2/3)), theta2.lim, theta1.lim)
)
opt$par

## another way of parameter specification
#interval <- list(theta2.1.lim, theta2.2.lim, theta1.1.lim, theta1.2.lim)
#system.time(
#opt <- ml.ql(yuima, h=1/((n)^(2/3)), param=param, interval=interval)
#)
#opt$par

system.time(
opt <- ml.ql(yuima, theta2, theta1, h=1/((n)^(2/3)), theta2.lim, theta1.lim, BFGS=TRUE)
)
opt$par

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
