\name{qmle}
%\alias{ql}
%\alias{ml.ql}
\alias{qmle}
\alias{quasilogl}
%\alias{LSE}
%\alias{ml.ql2}
\alias{rql}
\alias{lse}
%\alias{ql,ANY-method} 
%\alias{ml.ql,ANY-method} 
%\alias{ml.ql2,ANY-method} 
%\alias{rql,ANY-method} 
\title{Calculate quasi-likelihood and ML estimator of least squares estimator}
\description{Calculate the quasi-likelihood and estimate of the parameters of the
  stochastic differential equation by the maximum likelihood method or least squares estimator
  of the drift parameter.}
\usage{
%ml.ql(yuima,theta2,theta1,h,theta2.lim=matrix(c(0,1),1,2),theta1.lim=matrix(c(0,1),1,2),print=FALSE,method,param,interval)
%ql(yuima,theta2,theta1,h,print=FALSE,param)
%rql(yuima,theta2,theta1,ptheta2,ptheta1,h,print=FALSE,param,prevparam)
qmle(yuima, start, method="BFGS", fixed = list(), print=FALSE, lower, upper, joint=FALSE, ...)
quasilogl(yuima, param, print=FALSE)
}
\arguments{
  \item{yuima}{a yuima object.}
%  \item{theta2,theta1}{parameters of the sdeModel.}
%  \item{h}{ time span of observations.}
%  \item{theta2.lim, theta1.lim}{matrixes to specify the domains of the
%    parameters. Vector can be available only if theta is a scalar.}
%  \item{ptheta2,ptheta1}{}
  \item{print}{you can see a progress of the estimation when print is TRUE.}
  \item{method}{see Details.}
  \item{param}{\code{list} of parameters for the  quasi loglikelihood.}
%  \item{interval}{}
%  \item{prevparam}{}
  \item{lower}{a named list for specifying lower bounds of parameters}
  \item{upper}{a named list for specifying upper bounds of parameters}
  \item{start}{initial values to be passed to the optimizer.}
  \item{fixed}{for conditional (quasi)maximum likelihood estimation.}
  \item{joint}{perform joint estimation or two stage estimation? by default \code{joint=FALSE}.}
  \item{...}{passed to \code{\link{optim}} method. See Examples.}
}
\details{
%   A function ql calculate the quasi-likelihood of a time series data X with any
%   parameters. A function ml.pl estimates parameters of the sdeModel by
%   maximizing the quasi-likelihood.
  \code{qmle} behaves more likely the standard \code{mle} function in \pkg{stats4} and
  argument \code{method} is one of the methods available in \code{\link{optim}}.

  \code{lse} calculates least squares estimators of the drift parameters. This is
  useful for initial guess of \code{qmle} estimation.
  
  \code{quasilogl} returns the valueof the  quasi loglikelihood for a given
  \code{yuima} object and list of parameters \code{coef}.
  
}
\value{
  \item{QL}{a real value.}
  \item{opt}{a list with components the same as 'optim' function.}
}
\author{The YUIMA Project Team}
\note{
  %The function ml.ql uses the function optim internally.
  The function qmle uses the function optim internally.
}
\examples{
#dXt^e = -theta2 * Xt^e * dt + theta1 * dWt
diff.matrix <- matrix(c("theta1"), 1, 1)
ymodel <- setModel(drift=c("(-1)*theta2*x"), diffusion=diff.matrix, time.variable="t", state.variable="x", solve.variable="x")
n <- 100

ysamp <- setSampling(Terminal=(n)^(1/3), n=n) 
yuima <- setYuima(model=ymodel, sampling=ysamp)
set.seed(123)
yuima <- simulate(yuima, xinit=1, true.parameter=list(theta1=0.3,
theta2=0.1))
QL <- quasilogl(yuima, param=list(theta2=0.8, theta1=0.7))
##QL <- ql(yuima, 0.8, 0.7, h=1/((n)^(2/3)))
QL

## another way of parameter specification
##param <- list(theta2=0.8, theta1=0.7)
##QL <- ql(yuima, h=1/((n)^(2/3)), param=param)
##QL


## old code
##system.time(
##opt <- ml.ql(yuima, 0.8, 0.7, h=1/((n)^(2/3)), c(0, 1), c(0, 1))
##)
##cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
##print(coef(opt))


system.time(
opt2 <- qmle(yuima, start=list(theta1=0.8, theta2=0.7), lower=list(theta1=0,theta2=0), 
 upper=list(theta1=1,theta2=1), method="L-BFGS-B")
)
cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
print(coef(opt2))

## initial guess for theta2 by least squares estimator
tmp <- lse(yuima, start=list(theta2=0.7), lower=list(theta2=0), upper=list(theta2=1))
tmp

system.time(
opt3 <- qmle(yuima, start=list(theta1=0.8, theta2=tmp), lower=list(theta1=0,theta2=0), 
 upper=list(theta1=1,theta2=1), method="L-BFGS-B")
)
cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
print(coef(opt3))


## perform joint estimation? Non-optimal, just for didactic purposes
system.time(
opt4 <- qmle(yuima, start=list(theta1=0.8, theta2=0.7), lower=list(theta1=0,theta2=0), 
 upper=list(theta1=1,theta2=1), method="L-BFGS-B", joint=TRUE)
)
cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
print(coef(opt4))

## old code
##system.time(
##opt <- ml.ql(yuima, 0.8, 0.7, h=1/((n)^(2/3)), c(0, 1), c(0, 1), method="Newton")
##)
##cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
##print(coef(opt))




###multidimension case
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
yuima <- simulate(yuima, xinit=c(1, 1), true.parameter=list(theta2.1=0.5, theta2.2=0.3, theta1.1=0.6, theta1.2=0.2))

## theta2 <- c(0.8, 0.2) #c(theta2.1, theta2.2)
##theta1 <- c(0.7, 0.1) #c(theta1.1, theta1.2)
## QL <- ql(yuima, theta2, theta1, h=1/((n)^(2/3)))
## QL

## ## another way of parameter specification
## #param <- list(theta2=theta2, theta1=theta1)
## #QL <- ql(yuima, h=1/((n)^(2/3)), param=param)
## #QL

## theta2.1.lim <- c(0, 1)
## theta2.2.lim <- c(0, 1)
## theta1.1.lim <- c(0, 1)
## theta1.2.lim <- c(0, 1)
## theta2.lim <- t( matrix( c(theta2.1.lim, theta2.2.lim), 2, 2) )
## theta1.lim <- t( matrix( c(theta1.1.lim, theta1.2.lim), 2, 2) )

## system.time(
## opt <- ml.ql(yuima, theta2, theta1, h=1/((n)^(2/3)), theta2.lim, theta1.lim)
## )
## opt@coef

system.time(
opt2 <- qmle(yuima, start=list(theta2.1=0.8, theta2.2=0.2, theta1.1=0.7, theta1.2=0.1),
 lower=list(theta1.1=.1,theta1.2=.1,theta2.1=.1,theta2.2=.1),
 upper=list(theta1.1=4,theta1.2=4,theta2.1=4,theta2.2=4), method="L-BFGS-B")
)
opt2@coef
summary(opt2)

## unconstrained optimization
system.time(
opt3 <- qmle(yuima, start=list(theta2.1=0.8, theta2.2=0.2, theta1.1=0.7, theta1.2=0.1))
)
opt3@coef
summary(opt3)


quasilogl(yuima, param=list(theta2.1=0.8, theta2.2=0.2, theta1.1=0.7, theta1.2=0.1))

##system.time(
##opt <- ml.ql(yuima, theta2, theta1, h=1/((n)^(2/3)), theta2.lim, theta1.lim, method="Newton")
##)
##opt@coef
##}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}