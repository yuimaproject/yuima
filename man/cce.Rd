\name{cce}
\alias{cce}

\title{Nonsynchronous Cumulative Covariance Estimator}
\description{
This function estimates the covariance between two Ito processes when they are
observed at discrete times possibly nonsynchronously. It can apply to 
irregularly sampled one-dimensional data as a special case.
}
\usage{
cce(x)
}
\arguments{
  \item{x}{an object of class 'yuima' or 'yuima.data'}
}
\details{
  This function is a method of class \code{\link{yuima.data-class}} and \code{\link{yuima-class}}.
  It extracts the \code{data} slot when applied to a an object of
  class \code{\link{yuima-class}}.
}
\value{
  \item{value}{estimate of the covariance matrix}
}
\note{The example shows the central limit theorem for the nonsynchronous 
covariance estimator.  Estimation of the asymptotic variance with p-values and
the second-order correction will be provided in a future version of the package. 
}
\references{
 Hayashi, T and Yoshida, N. (2005)
  On covariance estimation of non-synchronously observed diffusion processes,
  \emph{Bernoulli}, \bold{11}, no. 2, 359--379.

 Hayashi, T and Yoshida, N. (2008) Asymptotic normality of a covariance 
  estimator for  nonsynchronously observed diffusion processes, \emph{Annals 
  of the Institute of Statistical Mathematics}, \bold{60}, no. 2, 367--406.
}
\author{The YUIMA Project Team}
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
