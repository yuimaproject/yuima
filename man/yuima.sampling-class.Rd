\name{yuima.sampling-class}
\docType{class}
\alias{yuima.sampling-class}
\alias{initialize,yuima.sampling-method}
\title{Classes for stochastic differential equations sampling scheme}
\description{
  The \code{yuima.sampling} class is a class of the \pkg{yuima} package.  
}
\section{Slots}{
  \describe{
	\item{\code{Initial}:}{initial time of the grid.}
	\item{\code{Terminal}:}{terminal time fo the grid.}
	\item{\code{n}:}{the number of observations - 1.}
	\item{\code{delta}:}{in case of a regular grid is the mesh.}
	\item{\code{grid}:}{the grid of times.}
	\item{\code{random}:}{either \code{FALSE} or the distribution 
	 of the random times.}
	\item{\code{scale}:}{by default 1. NOT YET IN USE.}
	\item{\code{regular}:}{indicator of whether the grid is regular or not.
	 For internal use only.}
  }
}
\details{
This object is created by \code{\link{setSampling}} or as a result of a 
simulation scheme by the \code{\link{simulate}} function.
}
\author{The YUIMA Project Team}
\keyword{classes}
