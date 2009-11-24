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
    \item{\code{Terminal}:}{Horizon of the stochastic processes.}
    \item{\code{division}:}{Time division for the data and/or model.}
    \item{\code{Initial}:}{Initial time of the data and/or model.}
    \item{\code{grid}:}{a grid of time points on which data are simulated.}
    \item{\code{random}:}{whether simulation consists of random times.}
  }
}
\author{YUIMA Project}
\keyword{classes}
