\name{setYuima}
\alias{setYuima}
\title{
Creates a "yuima" object by combining "model", "data", "sampling", "characteristic" 
and "functional"slots.
}
\description{
  \code{setYuima} constructs an object of  \code{\link{yuima-class}}.
}
\usage{
setYuima(data, model, sampling, characteristic, functional)
}
\arguments{
  \item{data}{an object of  \code{\link{yuima.data-class}}.}
  \item{model}{an object of  \code{\link{yuima.model-class}}.}
  \item{sampling}{an object of  \code{\link{yuima.sampling-class}}.}
  \item{characteristic}{an object of  \code{\link{yuima.characteristic-class}}.}
  \item{functional}{an object of class \code{\link{yuima.functional-class}}.}
}
\details{
The \code{yuima-class} object is the main object of the \pkg{yuima} package.
Some of the slots can be missing. 

The slot \code{data}  contains the data, either empirical or simulated. 

The slot \code{model}  contains the description of the 
(statistical) model which is used to generate the \code{data} via different 
simulation schemes, to draw inference from the \code{data} or both.

The \code{sampling} slot contains information on how the \code{data} have been
collected or how they should be simulated.

The slot \code{characteristic} contains information on \code{PLEASE FINISH THIS}.
The slot \code{functional} contains information on \code{PLEASE FINISH THIS}.

Please refer to the vignettes and the examples in the \pkg{yuimadocs}
package for more informations.
}
\value{
  an object of  \code{\link{yuima-class}}.
}
\author{The YUIMA Project Team}
\examples{
# Creation of a yuima object with all slots for a 
# stochastic differential equation
# dXt^e = -theta2 * Xt^e * dt + theta1 * dWt 
diffusion <- matrix(c("theta1"), 1, 1)
drift <- c("-1*theta2*x")
ymodel <- setModel(drift=drift, diffusion=diffusion)
n <- 100
ysamp <- setSampling(Terminal=1, n=n)

yuima <- setYuima(model=ymodel, sampling=ysamp)

str(yuima)
}
\keyword{ts}
