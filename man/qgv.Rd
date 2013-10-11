\name{qgv}
\alias{qgv}
\title{qgv}
\description{Estimate the local Holder exponent with quadratic generalized variations method}
\usage{
qgv(yuima, filter.type = "Daubechies", order = 2, a = NULL, ...)
}
\arguments{
  \item{yuima}{A \code{yuima} object.}
  \item{filter.type}{ The \code{filter.type} can be set to \code{"Daubechies"} or 
\code{"Classical"}.
}
  \item{order}{The order of the filter \code{a} to be chosen}
}
  \item{a}{Any other filter}
}
}
\details{
Estimation of the Husrt index and the constant of the fractional Ornstein-Uhlenmbeck process.
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
#Example
}
% R documentation directory.
\keyword{ts}