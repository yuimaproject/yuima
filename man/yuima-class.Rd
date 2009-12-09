\name{yuima-class}
\docType{class}
\alias{yuima-class}

\alias{get.zoo.data,yuima-method}
\alias{plot,yuima-method}
\alias{dim,yuima-method}
\alias{length,yuima-method}
\alias{cce,yuima-method}
\alias{initialize,yuima-method}
\alias{simulate,yuima-method}
\alias{poisson.random.sampling,yuima-method}
\alias{ql,yuima-method}
\alias{rql,yuima-method}
\alias{ml.ql,yuima-method}
\alias{limiting.gamma,yuima-method}

\title{Class for stochastic differential equations}
\description{
  The \code{yuima} S4 class is a class of the \pkg{yuima} package.  
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{an object of class \code{\link{yuima.data-class}}}
	\item{\code{model}:}{an object of class \code{\link{yuima.model-class}}}
	\item{\code{sampling}:}{an object of class 
	 \code{\link{yuima.sampling-class}}}
    \item{\code{characteristic}:}{an object of class 
	 \code{\link{yuima.characteristic-class}}}
  }
}
\details{
The \code{yuima-class} object is the main object of the \pkg{yuima} package.
Some of the slots may be missing. 

The \code{data} slot contain the data, either empirical or simulated. 

The \code{model} contains the description of the 
(statistical) model which is used to generate the \code{data} via different 
simulation schemes, to draw inference from the \code{data} or both.

The \code{sampling} slot contains information on how the \code{data} have been
collected or how they should be generated.

The slot \code{characteristic} contains information on 
 \code{PLEASE FINISH THIS}.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(x = "yuima")}: intialization method for 
	 \code{yuima-class} object.}
    \item{get.data}{\code{signature(x = "yuima")}: returns \code{data}.}
    \item{plot}{\code{signature(x = "yuima")}: same as for the list of 
	 \code{\link{zoo}}.}
    \item{dim}{\code{signature(x = "yuima")}: the number of SDEs in the
	  \code{yuima} object.}
    \item{length}{\code{signature(x = "yuima")}: a vector of length of
	  each SDE described in the \code{yuima} object.}
    \item{cce}{\code{signature(x = "yuima")}: Cumulative Covariance Estimator.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
