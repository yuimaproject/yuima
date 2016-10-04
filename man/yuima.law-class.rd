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

\title{Class of yuima law}

\description{Insert Description Here}

\section{Slots}{...:
\describe{
  \item{rng}{
  function}
    \item{density}{
  function}
    \item{cdf}{
  function}
    \item{quantile}{function}
    \item{characteristic}{
  function
  }
  \item{param.measure}{...}
    \item{time.var}{
  label}
    \item{dim}{
  number}
}
}

\section{Methods}{
  \describe{
    \item{rand}{\code{signature(object = "yuima.law", n = "numeric",
	  param = "list", ...)}: INSERT HERE DESCRIPTION}
  \item{dens}{\code{signature(object = "yuima.law", x = "numeric",
	  param = "list", log = FALSE, ...)}: INSERT HERE DESCRIPTION}
	\item{cdf}{\code{signature(object = "yuima.law", q = "numeric",
	  param = "list", ...)}: INSERT HERE DESCRIPTION}
	\item{quant}{\code{signature(object = "yuima.law", p = "numeric",
	  param = "list", ...)}: INSERT HERE DESCRIPTION}
	\item{char}{\code{signature(object = "yuima.law", u = "numeric",
	  param = "list", ...)}: INSERT HERE DESCRIPTION}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}

