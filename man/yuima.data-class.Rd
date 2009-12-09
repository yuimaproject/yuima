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
  The \code{yuima.data-class} is a class of the \pkg{yuima} package.  
}
\section{Slots}{
  \describe{
    \item{\code{original.data}:}{The original data.}
	\item{\code{zoo.data}:}{A list of \code{zoo} format data.}
  }
}
\details{
This object is one of the slots of the \code{\link{yuima.data-class}}.
The object contains two slots. 

The slot \code{original.data} contains, as the name suggests,
a copy of the original data passed by the user to the \code{initialize} method 
and/or the \code{\link{setData}} function. It is intended for backup purposes.

Internally, the \pkg{yuima} package stores and operates on  
\code{\link{zoo}}-type objects.
The \code{\link{zoo}}-fied data are contained in the slot \code{zoo.data}.

If data are obtained by simulation, the \code{original.data} slot 
is usually empty.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(x = "yuima.data")}: sets both slot of \code{\link{yuima.data-class}}.}
    \item{get.zoo.data}{\code{signature(x = "yuima.data")}: returns slot \code{zoo.data}.}
    \item{plot}{\code{signature(x = "yuima.data")}: same as \code{\link{zoo}}.}
    \item{dim}{\code{signature(x = "yuima.data")}: same as \code{\link{zoo}}.}
    \item{length}{\code{signature(x = "yuima.data")}: same as \code{\link{zoo}}.}
    \item{cce}{\code{signature(x = "yuima.data")}: Cumulative Covariance Estimator.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
