\name{subsampling}
\alias{subsampling}
\title{subsampling }
\description{subsampling}
\usage{
subsampling(x, sampling=NULL, Initial, Terminal, delta, 
 grid = as.numeric(NULL), random = FALSE, sdelta=as.numeric(NULL), 
 sgrid=as.numeric(NULL), interpolation="none") 
}
\arguments{
  \item{x}{an \code{\link{yuima-class}} or 
    \code{\link{yuima.model-class}} object.}
  \item{sampling}{a \code{\link{yuima.sampling-class}} object.}
  \item{Initial}{Initial time of the grid.}
  \item{Terminal}{Terminal time of the grid.}
  \item{delta}{mesh size in case of regular time grid.}
  \item{grid}{a grid of times for the simulation, possibly empty.}
  \item{random}{specify if it is random sampling. See Details.}
  \item{sdelta}{mesh size in case of regular space grid.}
  \item{sgrid}{a grid in space for the simulation, possibly empty.}
  \item{interpolation}{a rule of interpolation in case of subsampling. 
   By default, no interpolation. See Details.}
}
\value{
  \item{yuima}{a \code{yuima.data-class} object.}
}
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

plot(yuima.sim, plot.type="single")

## random sampling with exponential times

newsamp <- setSampling(
 random=list(rdist=c( function(x) rexp(x, rate=10), 
  function(x) rexp(x, rate=20))))
newdata <- subsampling(yuima.sim, sampling=newsamp)
points(get.zoo.data(newdata)[[1]],col="red")
points(get.zoo.data(newdata)[[2]],col="green")
}
\author{The YUIMA Project Team}
\keyword{ts}
