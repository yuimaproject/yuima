\name{yuima-class}
\docType{class}
\alias{yuima-class}

\alias{get.zoo.data,yuima-method}
\alias{plot,yuima-method}
%\alias{start,yuima-method}
%\alias{end,yuima-method}
%\alias{time,yuima-method}
\alias{dim,yuima-method}
\alias{length,yuima-method}
\alias{cce,yuima-method}
\alias{initialize,yuima-method}
%\alias{setSampling,yuima-method}
\alias{simulate,yuima-method}
\alias{poisson.random.sampling,yuima-method}
\alias{ql,yuima-method}
\alias{rql,yuima-method}
\alias{ml.ql,yuima-method}
\alias{limiting.gamma,yuima-method}

\title{Class for stochastic differential equations}
\description{
  The \code{yuima} class is a class of the \pkg{yuima} package.  
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{add explanation}
	\item{\code{model}:}{add explanation}
	\item{\code{sampling}:}{add explanation}
    \item{\code{characteristic}:}{add explanation}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(x = "yuima")}: intialization method for \code{yuima} object.}
    \item{get.data}{\code{signature(x = "yuima")}: returns \code{data}.}
    \item{plot}{\code{signature(x = "yuima")}: same as for the list of \code{zoo} or \code{ts}.}
    %\item{start}{\code{signature(x = "yuima")}: same as \code{zoo} or \code{ts}.}
    %\item{end}{\code{signature(x = "yuima")}: same as \code{zoo} or \code{ts}.}
    %\item{time}{\code{signature(x = "yuima")}: same as \code{zoo} or \code{ts}.}
    \item{dim}{\code{signature(x = "yuima")}: the number of SDEs in the
	  \code{yuima} object.}
    \item{length}{\code{signature(x = "yuima")}: a vector of length of
	  each SDE described in the \code{yuima} object.}
    \item{cce}{\code{signature(x = "yuima")}: Cumulative Covariance Estimator.}
  }
}
\author{YUIMA Project}
\keyword{classes}
