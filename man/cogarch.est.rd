\name{cogarch.est.-class}
\docType{class}
\alias{cogarch.est-class}
\alias{plot,cogarch.est.,ANY-method}
%%\alias{setSampling,yuima.carma-method}

\title{Class for Generalized Method of Moments Estimation for COGARCH(p,q) model}
\description{
  The \code{cogarch.est} class is a class of the  \pkg{yuima} package that contains estimated parameters obtained by the function \code{\link{gmm}} or \code{\link{qmle}}.
}
\section{Slots}{
  \describe{
    \item{\code{yuima}:}{is an object of of \code{\link{yuima-class}}.}
    \item{\code{objFun}:}{is an object of class \code{character} that indicates the objective function used in the minimization problem. See the documentation of the function \code{\link{gmm}} or \code{\link{qmle}} for more details.}
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
    \item{Methods mle}{All methods for \code{mle-class} are available.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
