\name{yuima.CP.qmle-class}
\docType{class}
\alias{yuima.CP.qmle-class}
\alias{yuima.qmle-class}
\alias{plot,yuima.CP.qmle,ANY-method}
\alias{qmle.CP}
\alias{CP.qmle}

\title{Class for Quasi Maximum Likelihood Estimation of Compound Poisson-based and SDE models}
\description{
  The \code{yuima.CP.qmle} class is a class of the  \pkg{yuima} package that extends the \code{mle-class} of the \pkg{stats4} package.
}
\section{Slots}{
  \describe{
    \item{\code{Jump.times}:}{a vector which contains the estimated time of jumps.}
    \item{\code{Jump.values}:}{a vector which contains the jumps.}
    \item{\code{X.values}:}{the value of the process at the jump times.}
    \item{\code{model}:}{is an object of of \code{\link{yuima.model-class}}.}
    \item{\code{call}:}{is an object of class \code{language}. }
    \item{\code{coef}:}{is an object of class \code{numeric} that contains estimated parameters.}
    \item{\code{fullcoef}:}{is an object of class \code{numeric} that contains estimated and fixed parameters.}
    \item{\code{vcov}:}{is an object of class \code{matrix}.}
    \item{\code{min}:}{is an object of class \code{numeric}.}
    \item{\code{minuslogl}:}{is an object of class \code{function}.}
    \item{\code{method}:}{is an object of class \code{character}.}
    \item{\code{model}:}{is an object of class \code{yuima.model-class}.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{Plot method for plotting the jump times.}
    \item{Methods mle}{All methods for \code{mle-class} are available.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
