\name{yuima.law-class}
\docType{class}
\alias{yuima.law-class}
\alias{yuima.law}
\alias{initialize,yuima.law-method}
\alias{rand,yuima.law-method}
\alias{dens,yuima.law-method}
\alias{cdf,yuima.law-method}
\alias{quant,yuima.law-method}
\alias{char,yuima.law-method}

\title{\code{yuima law-class}: A Mathematical Description for the Noise}

\description{A yuima class that contains all information on the noise. This class is a  bridge between a \code{\link{yuima.model-class}} and a noise constructed by users.}

\section{Slots}{
\describe{
  \item{rng}{
  A user defined function that generates the noise sample.}
    \item{density}{
  A user defined function that is the density of the noise at time \code{t}.}
    \item{cdf}{
  A user defined function that is the cumulative distribution function of the noise at time \code{t}.}
    \item{quantile}{A user defined function that is the quantile of the noise at time \code{t}.}
    \item{characteristic}{
  A user defined function that is the characteristic function of the noise at time \code{t}.
  }
  \item{param.measure}{A \code{character} object that contains the parameters of the noise.}
    \item{time.var}{
  the label of the time variable.}
    \item{dim}{
  Dimension of the noise}
}
}

\section{Methods}{
  \describe{
    \item{rand}{\code{signature(object = "yuima.law", n = "numeric",
	  param = "list", ...)}: This method returns a sample of the noise, \code{n} is the sample size.}
  \item{dens}{\code{signature(object = "yuima.law", x = "numeric",
	  param = "list", log = FALSE, ...)}: This method returns the density of the noise, \code{x} is the vector of the support.}
	\item{cdf}{\code{signature(object = "yuima.law", q = "numeric",
	  param = "list", ...)}: This method returns the cdf of the noise, \code{q} is the vector of the support.}
	\item{quant}{\code{signature(object = "yuima.law", p = "numeric",
	  param = "list", ...)}: This method returns the quantile of the noise, \code{p} is the vector of the support.}
	\item{char}{\code{signature(object = "yuima.law", u = "numeric",
	  param = "list", ...)}: This method returns the characteristic function of the noise, \code{u} is the vector of the support.}
  }
}
\author{The YUIMA Project Team

Contacts: Lorenzo Mercuri \email{lorenzo.mercuri@unimi.it}
}
\keyword{classes}

