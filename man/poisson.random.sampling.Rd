\name{poisson.random.sampling}
\alias{poisson.random.sampling}
\title{Poisson random sampling method} 
\description{
Poisson random sampling method.}
\usage{
poisson.random.sampling(x, rate, n)
}
\arguments{
  \item{x}{an object of  \code{\link{yuima.data-class}} or \code{\link{yuima-class}}.}
  \item{rate}{a Poisson intensity or a vector of Poisson intensities.}
  \item{n}{a common multiplier to the Poisson intensities. The default value is 1.}
}
\details{
  It returns an object
  of type \code{\link{yuima.data-class}} which is a copy of the original input 
  data where observations are sampled according to
  the Poisson process. The unsampled data are set to \code{NA}.
}
\value{
  an object of  \code{\link{yuima.data-class}}.
}
\author{The YUIMA Project Team}
\seealso{\link{cce}}
\examples{
## Set a model
diff.coef.1 <- function(t, x1=0, x2) x2*(1+t)
diff.coef.2 <- function(t, x1, x2=0) x1*sqrt(1+t^2)
cor.rho <- function(t, x1=0, x2=0) sqrt((1+cos(x1*x2))/2)
diff.coef.matrix <- matrix(c("diff.coef.1(t,x1,x2)",
"diff.coef.2(t,x1,x2)*cor.rho(t,x1,x2)", "",
"diff.coef.2(t,x1,x2)*sqrt(1-cor.rho(t,x1,x2)^2)"),2,2)
cor.mod <- setModel(drift=c("",""), diffusion=diff.coef.matrix, 
solve.variable=c("x1", "x2"), xinit=c(3,2))
set.seed(111)

## We first simulate the two dimensional diffusion model
yuima.samp <- setSampling(Terminal=1, n=1200)
yuima <- setYuima(model=cor.mod, sampling=yuima.samp)
yuima.sim <- simulate(yuima)

## Then we use function poisson.random.sampling to get observations
## by Poisson sampling.
psample <- poisson.random.sampling(yuima.sim,  rate = c(0.2, 0.3), n=1000)
str(psample)
}
\keyword{ts}
