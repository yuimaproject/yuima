\name{yuima.qmleLevy.incr}
\docType{class}
\alias{yuima.qmleLevy.incr-class}
\alias{yuima.qmleLevy.incr,ANY-method}
\alias{initialize,yuima.qmleLevy.incr-method}
\alias{qmleLevy.incr}
\alias{incr.qmleLevy}

\title{Class for Quasi Maximum Likelihood Estimation of Levy SDE model}
\description{
  The \code{yuima.qmleLevy.incr-class} is a class of the  \pkg{yuima} package that extends the \code{mle-class} of the \pkg{stats4} package.  
}
\section{Slots}{
  \describe{
    \item{\code{Incr.Lev}:}{is an object of class \code{\link{yuima.data-class}} that contains the estimated increments of the noise.}
    \item{\code{logL.Incr}:}{an numeric object that represents the value of the loglikelihood for the estimated Levy increments.}
    \item{\code{minusloglLevy}:}{an R function that evaluates the loglikelihood of the estimated Levy increments. The function is used internally in \code{\link{qmleLevy}} for the estimation of the Levy measure parameters.}
    \item{\code{Levydetails}:}{a \code{list} containing additional information about the optimization procedure in the estimation of the Levy measure parameters. See \code{\link{optim}} help for the meaning of the components of this \code{list}.}
      \item{\code{Data}:}{is an object of \code{\link{yuima.data-class}} containing observation data.}
    \item{\code{model}:}{is an object of of \code{\link{yuima.carma-class}}.}
    \item{\code{call}:}{is an object of class \code{language}. }
    \item{\code{coef}:}{is an object of class \code{numeric} that contains estimated parameters.}
    \item{\code{fullcoef}:}{is an object of class \code{numeric} that contains estimated and fixed parameters.}
    \item{\code{vcov}:}{is an object of class \code{matrix}.}
    \item{\code{min}:}{is an object of class \code{numeric}.}
    \item{\code{minuslogl}:}{is an object of class \code{function}.}
    \item{\code{nobs}:}{an object of class \code{numeric}.}
    \item{\code{method}:}{is an object of class \code{character}.}
  }
}
\section{Methods}{
  \describe{
    \item{Methods mle}{All methods for \code{mle-class} are available.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}


