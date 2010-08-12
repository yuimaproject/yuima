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
  \item{x}{an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.}
}
\details{
  This function is a method for objects of \code{\link{yuima.data-class}} and \code{\link{yuima-class}}.
  It extracts the \code{data} slot when applied to a an object of
   \code{\link{yuima-class}}.
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
diff.coef.matrix <- matrix(c("diff.coef.1(t,x1,x2)", "diff.coef.2(t,x1,x2) * cor.rho(t,x1,x2)", 
"", "diff.coef.2(t,x1,x2) * sqrt(1-cor.rho(t,x1,x2)^2)"), 2, 2) 
cor.mod <- setModel(drift = c("", ""), diffusion = diff.coef.matrix,solve.variable = c("x1", "x2")) 

set.seed(111) 

## We use a function poisson.random.sampling to get observation by Poisson sampling.
yuima.samp <- setSampling(Terminal = 1, n = 1200) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima) 
psample<- poisson.random.sampling(yuima, rate = c(0.2,0.3), n = 1000) 

## cce takes the psample and returns an estimate of the quadratic covariation. 
cce(psample)$covmat[1, 2] 

## True value of the quadratic covariation.
cc.theta <- function(T, sigma1, sigma2, rho) { 
	tmp <- function(t) return(sigma1(t) * sigma2(t) * rho(t)) 
	integrate(tmp, 0, T) 
}

theta <- cc.theta(T = 1, diff.coef.1, diff.coef.2, cor.rho)$value 
cat(sprintf("theta =\%.5f\n", theta))

names(psample@zoo.data)






# Example. A stochastic differential equation with nonlinear feedback. 

## Set a model
drift.coef.1 <- function(x1,x2) x2
drift.coef.2 <- function(x1,x2) -x1
drift.coef.vector <- c("drift.coef.1","drift.coef.2")
diff.coef.1 <- function(t,x1,x2) sqrt(abs(x1))*sqrt(1+t)
diff.coef.2 <- function(t,x1,x2) sqrt(abs(x2))
cor.rho <- function(t,x1,x2) 1/(1+x1^2)
diff.coef.matrix <- matrix(c("diff.coef.1(t,x1,x2)", "diff.coef.2(t,x1,x2) * cor.rho(t,x1,x2)","", "diff.coef.2(t,x1,x2) * sqrt(1-cor.rho(t,x1,x2)^2)"), 2, 2) 
cor.mod <- setModel(drift = drift.coef.vector, diffusion = diff.coef.matrix,solve.variable = c("x1", "x2"))

## Generate a path of the process
set.seed(111) 
yuima.samp <- setSampling(Terminal = 1, n = 10000) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima, xinit=c(2,3)) 
plot(yuima)


## The "true" value of the quadratic covariation.
cce(yuima)

## We use the function poisson.random.sampling to generate nonsynchronous observations by Poisson sampling.
psample<- poisson.random.sampling(yuima, rate = c(0.2,0.3), n = 3000) 

## cce takes the psample to return an estimated value  of the quadratic covariation. The off-diagonal elements are the value of the Hayashi-Yoshida estimator. 
cce(psample)




}
\seealso{\link{setModel}}
\keyword{ts}
