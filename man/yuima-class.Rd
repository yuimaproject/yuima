\name{yuima-class}
\docType{class}
\alias{yuima-class}
\alias{get.zoo.data,yuima-method}
\alias{plot,yuima,ANY-method}
\alias{dim,yuima-method}
\alias{length,yuima-method}
\alias{cce,yuima-method}
\alias{llag,yuima-method}
\alias{initialize,yuima-method}
\alias{simulate,yuima-method}
\alias{poisson.random.sampling,yuima-method}
\alias{subsampling,yuima-method}
\alias{ql,yuima-method}
\alias{rql,yuima-method}
\alias{ml.ql,yuima-method}
\alias{limiting.gamma,yuima-method}
\alias{LSE,yuima-method}
\alias{cbind,yuima,ANY-method}

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
    \item{\code{functional}:}{an object of class
	 \code{\link{yuima.functional-class}}}
  }
}
\details{
The \code{yuima-class} object is the main object of the \pkg{yuima} package.
Some of the slots may be missing.

The \code{data} slot contains the data, either empirical or simulated.

The \code{model} contains the description of the
(statistical) model which is used to generate the \code{data} via different
simulation schemes, to draw inference from the \code{data} or both.

The \code{sampling} slot contains information on how the \code{data} have been
collected or how they should be generated.

The slot \code{characteristic} contains information on
 \code{PLEASE FINISH THIS}.
The slot \code{functional} contains information on
 \code{PLEASE FINISH THIS}.
}
\section{Methods}{
  \describe{
    \item{new}{\code{signature(x = "yuima", data = "yuima.data",
	  model = "yuima.model", sampling = "yuima.sampling",
	  characteristic = "yuima.characteristic"}: the function makes a copy of
	  the prototype object from the class definition of
	  \code{\link{yuima-class}}, then calls the \code{initialize} method
	  passing as arguments the newly created object and the remaining
	  arguments.}
	\item{initialize}{\code{signature(x = "yuima", data = "yuima.data",
	  model = "yuima.model", sampling = "yuima.sampling",
	  characteristic = "yuima.characteristic"}:  makes a copy of each argument
	  in the corresponding slots of the object \code{x}.}
    \item{get.data}{\code{signature(x = "yuima")}: returns the content of the
	  slot \code{data}.}
    \item{plot}{\code{signature(x = "yuima", \dots)}: calls
	 \code{\link{plot}} from the \code{\link{zoo}} package with argument
	 \code{x@data@zoo.data}. Additional arguments \code{\dots} are passed
	 as is to the \code{\link{plot}} function.}
	\item{dim}{\code{signature(x = "yuima")}: the number of SDEs in the
	  \code{yuima} object.}
    \item{length}{\code{signature(x = "yuima")}: a vector of length of
	  each SDE described in the \code{yuima} object.}
    \item{cce}{\code{signature(x = "yuima")}: calculates the asyncronous
	 covariance estimator on the data contained in \code{x@data@zoo.data}.
	 For more details see \code{\link{cce}}.}
   \item{llag}{\code{signature(x = "yuima")}: calculates the lead lag estimate
	 r on the data contained in \code{x@data@zoo.data}.
	 For more details see \code{\link{llag}}.}
    \item{simulate}{simulation method. For more information see
	  \code{\link{simulate}}.}
    \item{cbind}{\code{signature(x = "yuima")}: bind yuima.data object.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
