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
    \item{Initial}{Initial time of the data and/or model.}
    \item{Terminal}{Horizon of the stochastic processes.}
    \item{division}{Time division for the data and/or model.}
    \item{grid}{a grid of time points on which data are simulated.}
    \item{random}{whether simulation consists of random times.}
	\item{delta}{the mesh size in case of regular sampling.}
	\item{regular}{whether the ampling is regular. For internal use only.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
