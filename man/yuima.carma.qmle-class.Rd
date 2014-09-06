\name{yuima.carma.qmle-class}
\docType{class}
\alias{yuima.carma.qmle-class}
\alias{plot,yuima.carma.qmle,ANY-method}
\alias{qmle.carma}
\alias{carma.qmle}
%%\alias{setSampling,yuima.carma-method}

\title{Class for Quasi Maximum Likelihood Estimation of CARMA(p,q) model}
\description{
  The \code{yuima.carma.qmle} class is a class of the  \pkg{yuima} package that extends the \code{mle-class} of the \pkg{stats4} package.  
}
\section{Slots}{
  \describe{
    \item{\code{Incr.Lev}:}{is an object of class \code{\link{zoo}} that contains the estimated increments of the noise obtained using \code{\link{CarmaNoise}}.}
    \item{\code{model}:}{is an object of of \code{\link{yuima.carma-class}}.}
    \item{\code{logL.Incr}:}{is an object of class \code{numeric} that contains the value of the log-likelihood for estimated Levy increments.}
    \item{\code{call}:}{is an object of class \code{language}. }
    \item{\code{coef}:}{is an object of class \code{numeric} that contains estimated parameters.}
    \item{\code{fullcoef}:}{is an object of class \code{numeric} that contains estimated and fixed parameters.}
    \item{\code{vcov}:}{is an object of class \code{matrix}.}
    \item{\code{min}:}{is an object of class \code{numeric}.}
    \item{\code{minuslogl}:}{is an object of class \code{function}.}
    \item{\code{method}:}{is an object of class \code{character}.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{Plot method for estimated increment of the noise.}
    \item{Methods mle}{All methods for \code{mle-class} are available.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
