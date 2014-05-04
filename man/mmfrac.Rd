\name{mmfrac}
\alias{mmfrac}
\title{
mmfrac
}
\description{
Estimates the drift of a fractional Ornstein-Uhlenbeck and, if necessary, also the Hurst and diffusion parameters.
}
\usage{
mmfrac(yuima, ...)
}
\arguments{
  \item{yuima}{a \code{yuima} object.
}
  \item{\dots}{arguments passed to \code{\link{qgv}}.
}
}
\details{
Estimates the drift of s fractional Ornstein-Uhlenbeck and, if necessary, also the Hurst and diffusion parameters.
}
\value{an object of class mmfrac
}
\references{
Brouste, A., Iacus, S.M. (2013) Parameter estimation for the discretely observed fractional Ornstein-Uhlenbeck process and the Yuima R package, Computational Statistics, pp. 1129--1147.
}
\author{The YUIMA Project Team}
\seealso{
See also \code{\link{qgv}}.
}
\examples{
# Estimating all Hurst parameter, diffusion coefficient  and drift coefficient 
# in fractional Ornstein-Uhlenbeck

model<-setModel(drift="-x*lambda",hurst=NA,diffusion="theta")
sampling<-setSampling(T=100,n=10000)
yui1<-simulate(model,true.param=list(theta=1,lambda=4),hurst=0.7,sampling=sampling)
mmfrac(yui1)

}
% R documentation directory.
\keyword{ts}
