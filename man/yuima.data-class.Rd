\name{yuima.data-class}
\docType{class}
\alias{yuima.data-class}

\alias{get.zoo.data,yuima.data-method}
\alias{plot,yuima.data-method}
%\alias{start,yuima.data-method}
%\alias{end,yuima.data-method}
%\alias{time,yuima.data-method}
\alias{dim,yuima.data-method}
\alias{length,yuima.data-method}
\alias{cce,yuima.data-method}
\alias{initialize,yuima.data-method}
%\alias{setSampling,yuima.data-method}
\alias{poisson.random.sampling,yuima.data-method}

\title{Classes for stochastic differential equations data object}
\description{
  The \code{yuima.data} class is a class of the \pkg{yuima} package.  
}
\section{Slots}{
  \describe{
    \item{\code{original.data}:}{The original data.}
	\item{\code{zoo.data}:}{A list of \code{zoo} format data.}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(x = "yuima.data")}: sets both slot of \code{yuima.data} class.}
    \item{get.zoo.data}{\code{signature(x = "yuima.data")}: returns slot \code{zoo.data}.}
    \item{plot}{\code{signature(x = "yuima.data")}: same as \code{zoo} or \code{ts}.}
   % \item{start}{\code{signature(x = "yuima.data")}: same as \code{zoo} or \code{ts}.}
   % \item{end}{\code{signature(x = "yuima.data")}: same as \code{zoo} or \code{ts}.}
   % \item{time}{\code{signature(x = "yuima.data")}: same as \code{zoo} or \code{ts}.}
    \item{dim}{\code{signature(x = "yuima.data")}: same as \code{zoo} or \code{ts}.}
    \item{length}{\code{signature(x = "yuima.data")}: same as \code{zoo} or \code{ts}.}
    \item{cce}{\code{signature(x = "yuima.data")}: Cumulative Covariance Estimator.}
  }
}
\author{YUIMA Project}
\keyword{classes}
