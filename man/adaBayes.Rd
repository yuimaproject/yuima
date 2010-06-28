\name{adaBayes}
\alias{adaBayes}
\alias{adaBayes,yuima-method}
\title{Adaptive Bayes estimator for the parameters in sde model}
\description{Adaptive Bayes estimator for the parameters in a specific type of sde.}
\usage{
adaBayes(yuima, print = FALSE, start, prior, propose, n.iter = 100, 
    lower, upper, n.burnin, method = "nomcmc", mhtype = "independent")
}
\arguments{
  \item{yuima}{a 'yuima' object.}
  \item{prior}{ a list of prior distributions for the parameters [function].}
  \item{start}{ initial suggestion for parameter values }
  \item{propose}{ ... }
  \item{print}{trace optimization?}
  \item{method}{ \code{nomcmc} requires package \code{adapt} }
  \item{mhtype}{ ... }
  \item{n.burnin}{ ... }
  \item{lower}{ ... }
  \item{upper}{ ... }
  \item{n.iter}{ ... }
}
\details{
  Calculate the values of the parameters theta1 and theta2. When the quasi-likelihood is too large for the given data, the integral in the Bayes estimator may diverges. For such a case, offset values for the parameters can be specified not to diverge.}
\value{
  \item{vector}{a vector of the paramter estimate}
}
\author{The YUIMA Project Team}
\note{
  This routine is not stable and accurate. Dr.Kamatani is now working for improvements.
}
\examples{
library(yuima)
library(mnormt)
set.seed(123)

ymodel <- setModel(drift="(-1)*theta2*x", diffusion="sqrt(theta1^2+1)",
                   time.variable="t", state.variable="x", solve.variable="x")
n <- 100#10^2
h <- 1/((n)^(2/3))
ysamp <- setSampling(Terminal=(n)^(1/3), n=n) 
yuima <- setYuima(model=ymodel, sampling=ysamp)
param.true <- list(theta2=0.3, theta1=0.5)
yuima <- simulate(yuima, xinit=1, true.parameter=param.true)

prior.theta1 <- function(theta2)
  1*(theta2 > 0 & theta2 < 1)

prior.theta2 <- function(theta1)
  1*(theta1 > 0 & theta1 < 1)

prior <- list(theta1=list(measure.type="density", density=prior.theta1,domain=c(-2,2)),
theta2=list(measure.type="density", density=prior.theta2,domain=c(-2,2)))

param.init <- list(theta2=0.35,theta1=0.52)

lower = c(0,0)
upper=c(1,1)

bayes1 <- adaBayes(yuima, start=param.init, lower=lower,upper=upper,prior=prior, method="nomcmc")

mle1 <- qmle(yuima, start=param.init, lower=lower, upper=upper, method="L-BFGS-B")
mle1@coef
}
\keyword{ts}
