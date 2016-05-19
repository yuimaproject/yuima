\name{adaBayes}
\alias{adaBayes}
\alias{adaBayes,yuima-method}
\title{Adaptive Bayes estimator for the parameters in sde model}
\description{Adaptive Bayes estimator for the parameters in a specific type of sde.}
\usage{
adaBayes(yuima, start, prior, lower, upper, method = "nomcmc")
}
\arguments{
  \item{yuima}{a 'yuima' object.}
  \item{start}{ initial suggestion for parameter values }
  \item{prior}{ a list of prior distributions for the parameters [function].}
  \item{lower}{a named list for specifying lower bounds of parameters}
  \item{upper}{a named list for specifying upper bounds of parameters}
  \item{method}{ \code{nomcmc} requires package \code{cubature} }
}
\details{
  Calculate the values of the parameters theta1 and theta2. When the quasi-likelihood is too large for the given data, the integral in the Bayes estimator may diverge. For such a case, offset values for the parameters can be specified not to diverge.}
\value{
  \item{vector}{a vector of the paramter estimate}
}
\author{The YUIMA Project Team}
\note{
  This routine is not stable and accurate. Dr.Kamatani is now working for improvements.
}
\references{
Yoshida, N. (2011). Polynomial type large deviation inequalities and quasi-likelihood analysis for stochastic differential equations. Annals of the Institute of Statistical Mathematics, 63(3), 431-479.

Uchida, M., & Yoshida, N. (2014). Adaptive Bayes type estimators of ergodic diffusion processes from discrete observations. Statistical Inference for Stochastic Processes, 17(2), 181-219.
}
\examples{
set.seed(123)

ymodel <- setModel(drift="(-1)*theta2*x", diffusion="sqrt(theta1^2+1)",
                   time.variable="t", state.variable="x", solve.variable="x")
n <- 500
h <- 1/((n)^(2/3))
ysamp <- setSampling(Terminal=(n)^(1/3), n=n) 
yuima <- setYuima(model=ymodel, sampling=ysamp)
param.true <- list(theta2=0.3, theta1=0.5)
yuima <- simulate(yuima, xinit=1, true.parameter=param.true)

prior.theta1 <- function(theta2)
  1*(theta2 > 0 & theta2 < 1)

prior.theta2 <- function(theta1)
  1*(theta1 > 0 & theta1 < 1)

prior <- list(theta2=list(measure.type="code",df="dunif(z,0,1)"), 
theta1=list(measure.type="code",df="dunif(z,0,1)"))

param.init <- list(theta2=0.35,theta1=0.52)

lower = c(0,0)
upper=c(1,1)

bayes1 <- adaBayes(yuima, start=param.init, prior=prior, method="nomcmc")
bayes1@coef

mle1 <- qmle(yuima, start=param.init, lower=list(theta1=0,theta2=0), 
     	upper=list(theta1=1,theta2=1), method="L-BFGS-B")
mle1@coef
}
\keyword{ts}
