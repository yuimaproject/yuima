\name{rng}
\alias{rng}
\alias{rIG}
\alias{rNIG}
\alias{rbgamma}
\alias{rngamma}
\alias{rstable}
\title{Random Number Generators for yuima package}
\description{\code{simulate} function uses specific random number generators to generate Levy paths.}
\usage{
rIG(x,delta,gamma)
rNIG(x,alpha=1,beta=0,delta=1,mu=0,Lambda)
rbgamma(x,delta.plus,gamma.plus,delta.minus,gamma.minus)
rngamma(x,lambda,alpha,beta,mu,Lambda)
rstable(x,alpha,beta,sigma,gamma)
}
\arguments{
  \item{x}{Number of R.Ns to be geneated.}
  \item{delta}{param}
  \item{delta.plus}{param}
  \item{delta.minus}{param}
  \item{gamma}{param}
  \item{gamma.plus}{param}
  \item{gamma.minus}{param}
  \item{alpha}{param}
  \item{beta}{parame}
  \item{mu}{param}
  \item{lambda}{param}
  \item{Lambda}{param}
  \item{sigma}{param}
}
\details{
  to be written
}
\value{
  \item{rIG}{A vector of random numbers.}
  \item{rNIG}{A vector of random numbers.}
  \item{rbgamma}{A vector of random numbers.} 
}
\author{The YUIMA Project Team}
\note{
  Random number generators may be rewritten with C language to speed up.
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
