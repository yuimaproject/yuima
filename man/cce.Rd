\name{cce}
\alias{cce}

\title{Nonsyncronous Cumulative Covariance Estimator}
\description{
\code{cce} is a function to estimate the covariance of some diffusion-type processes when they are observed only at discrete 
times in a nonsynchronous manner.}
\usage{
cce(x)
}
\arguments{
  \item{x}{An object of class 'yuima' or 'yuima.data'}
}
\details{
  This function is a method of class \code{yuima.data} or \code{yuima}.
  In the case of \code{yuima} object extracts the \code{data} object from
  it.
}
\value{
  \item{value}{returns a matrix of estimated covariance}
}
\references{
 Hayashi, T and Yoshida, N. (2005)
\emph{On covariance estimation of non-synchronously observed diffusion processes.}
\emph{Bernoulli} \bold{11}, no. 2, 359--379.
}
\author{The YUIMA Project Team}
\note{The example shows CLT for the non synchronous covariance estimator
}
\examples{
## Set a model
diff.coef.1 <- function(t, x1 = 0, x2 = 0) sqrt(1+t)
diff.coef.2 <- function(t, x1 = 0, x2 = 0) sqrt(1+t^2)
cor.rho <- function(t, x1 = 0, x2 = 0) sqrt(1/2)
diff.coef.matrix <- matrix(c("diff.coef.1(t,x1,x2)", "diff.coef.2(t,x1,x2) * cor.rho(t,x1,x2)", "", "diff.coef.2(t,x1,x2) * sqrt(1-cor.rho(t,x1,x2)^2)"), 2, 2) 
cor.mod <- setModel(drift = c("", ""), diffusion = diff.coef.matrix,solve.variable = c("x1", "x2")) 

set.seed(111) 

## We use a function poisson.random.sampling to get observation by Poisson sampling.
yuima.samp <- setSampling(Terminal = 1, division = 1200) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima) 
psample<- poisson.random.sampling(yuima, rate = c(0.2,0.3), n = 1000) 

## CCE takes the psample and returns an estimate of the quadratic covariation. 
cce(psample)[1, 2] 

## True value of the quadratic covariation.
cc.theta <- function(T, sigma1, sigma2, rho) { 
	tmp <- function(t) return(sigma1(t) * sigma2(t) * rho(t)) 
	integrate(tmp, 0, T) 
}

theta <- cc.theta(T = 1, diff.coef.1, diff.coef.2, cor.rho)$value 
cat(sprintf("theta =\%.5f\n", theta))

names(psample@zoo.data)
}
\seealso{\link{setModel}}
\keyword{ts}
