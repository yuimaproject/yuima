\name{yuima.th-class}
\docType{class}
\alias{yuima.th-class}
\alias{yuima.th}
\alias{t-Levy process}
\alias{initialize,yuima.th-method}
\alias{rand,yuima.th-method}
\alias{dens,yuima.th-method}
\alias{cdf,yuima.th-method}
\alias{quant,yuima.th-method}
\alias{char,yuima.th-method}


\title{\code{yuima.th-class}: A mathematical description for the t-Levy process.}

\description{A yuima class that contains all information on the noise for t-Levy process. This class extends \code{\link{yuima.law-class}} and contains info on the numerical method used for the inversion of the characteristic function. Three inversion methods are available: \code{cos}, \code{Laguerre} and \code{FFT}.

An object of this class can be created by calls of the function \code{\link{setLaw_th}}.
}

\section{Methods}{
  \describe{
    \item{rand}{\code{signature(object = "yuima.th", n = "numeric",
	  param = "list", ...)}: This method returns a sample of the noise, \code{n} is the sample size.}
  \item{dens}{\code{signature(object = "yuima.th", x = "numeric",
	  param = "list", log = FALSE, ...)}: This method returns the density of the noise, \code{x} is the vector of the support.}
	\item{cdf}{\code{signature(object = "yuima.th", q = "numeric",
	  param = "list", ...)}: This method returns the cdf of the noise, \code{q} is the vector of the support.}
	\item{quant}{\code{signature(object = "yuima.th", p = "numeric",
	  param = "list", ...)}: This method returns the quantile of the noise, \code{p} is the vector of the support.}
	\item{char}{\code{signature(object = "yuima.th", u = "numeric",
	  param = "list", ...)}: This method returns the characteristic function of the noise, \code{u} is the vector of the support.}
  }
}
\author{The YUIMA Project Team

Contacts: Lorenzo Mercuri \email{lorenzo.mercuri@unimi.it}
}
\keyword{classes}
