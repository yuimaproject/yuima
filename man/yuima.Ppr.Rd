\name{yuima.PPR}
\docType{class}
\alias{yuima.PPR}
\alias{yuima.PPR-class}
\alias{initialize,yuima.PPR-method}
\alias{simulate,yuima.PPR-method}

\title{Class for a mathematical description of a Point Process}
\description{
The \code{yuima.PPR} class is a class of the \pkg{yuima} package that extends the \code{\link{yuima-class}}. The object of this class contains all the information about the Point Process Regression Model.
}
\section{Objects from the Class}{
Objects can be created by calls of the function \code{\link{setPPR}}.
}
\section{Slots}{
  \describe{
     \item{\code{PPR}:}{is an object of class \code{\link{info.PPR}}.}
    \item{\code{gFun}:}{is an object of class \code{\link{info.Map}}.}
    \item{\code{Kernel}:}{is an object of class \code{\link{Integral.sde}}.}
    \item{\code{data}:}{is an object of class \code{\link{yuima.data-class}}. The slot contain either true data or simulated data}
    \item{\code{model}:}{is an object of class \code{\link{yuima.model-class}}. The slot contains all the information about the covariates}
    \item{\code{sampling}:}{is an object of class \code{\link{yuima.sampling-class}}.}
    \item{\code{characteristic}:}{is an object of class \code{\link{yuima.characteristic-class}}.}
    \item{\code{model}:}{is an object of class \code{\link{yuima.functional-class}}.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}

