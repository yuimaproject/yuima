\name{poisson.random.sampling}
\alias{poisson.random.sampling}
\title{Poisson random sampling method} 
\description{
Poisson random sampling method. A method of class 'yuima.data' and 'yuima'.}
\usage{
poisson.random.sampling(x, rate, n)
}
\arguments{
  \item{x}{An object of class 'yuima.data', 'yuima'}
  \item{rate}{poisson rate}
  \item{n}{To specify intensity as n * rate}
}
\details{
  It returns an object
  of type \code{yuima.data} which is a copy of the original input \code{data}
  (or its \code{zoo.data} field) where observations are sampled according to
  the poisson process. The remaining observations are set to \code{NA}.
}
\value{
  An object of class 'yuima.data'.
}
\references{Yoshida, Hayashi, etc. ADD CORRECT REFERENCE}
\author{The YUIMA Project Team}
\seealso{\link{cce}}
\keyword{ts}
