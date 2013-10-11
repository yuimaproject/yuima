\name{qgv}
\alias{qgv}
\title{qgv}
\description{Estimate the local Holder exponent with quadratic generalized variations method}
\usage{
qgv(yuima, filter.type = "Daubechies", order = 2, a = NULL)
}
\arguments{
  \item{yuima}{A \code{yuima} object.}
  \item{filter.type}{ The \code{filter.type} can be set to \code{"Daubechies"} or 
\code{"Classical"}.}
  \item{order}{The order of the filter \code{a} to be chosen}

  \item{a}{Any other filter}
}
\details{
Estimation of the Hurst index and the constant of the fractional Ornstein-Uhlenbeck process.
}
\value{an object of class qgv
}
\references{
Brouste, A., Iacus, S.M. (2013) Parameter estimation for the discretely observed fractional Ornstein-Uhlenbeck process and the Yuima R package, Computational Statistics, pp. 1129--1147.
}
\author{The YUIMA Project Team}
\seealso{
See also \code{\link{mmfrac}}.
}
\examples{
# Estimating both Hurst parameter and diffusion coefficient in fractional Ornstein-Uhlenbeck

model<-setModel(drift="-x*lambda",hurst=NA,diffusion="theta")
sampling<-setSampling(T=100,n=10000)
yui1<-simulate(model,true.param=list(theta=1,lambda=4),hurst=0.7,sampling=sampling)
qgv(yui1)


# Estimating Hurst parameter only in diffusion processes

model2<-setModel(drift="-x*lambda",hurst=NA,diffusion="theta*sqrt(x)")
sampling<-setSampling(T=1,n=10000)
yui2<-simulate(model2,true.param=list(theta=1,lambda=4),hurst=0.7,sampling=sampling,xinit=10)
qgv(yui2)
}
% R documentation directory.
\keyword{ts}