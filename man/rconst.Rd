\name{rconst}
\alias{rconst}
\alias{dconst}
\title{
Fictitious rng for the constant random variable used to generate and describe
Poisson jumps.
}
\description{
Fictitious rng for the constant random variable used to generate and describe
Poisson jumps.
}
\usage{
rconst(n, k = 1)
dconst(x, k = 1)
}
\arguments{
  \item{n}{number of replications}
  \item{k}{the size of the jump}
  \item{x}{the fictitious argument}
}
\value{
returns a numeric vector
}
\author{The YUIMA Project Team}
\examples{
dconst(1,1)
dconst(2,1)
dconst(2,2)

rconst(10,3)
}
\keyword{ts}