\name{cogarch.gmm.-class}
\docType{class}
\alias{cogarch.gmm-class}
\alias{plot,cogarch.gmm.,ANY-method}
%%\alias{setSampling,yuima.carma-method}

\title{Class for Generalized Method of Moments Estimation for COGARCH(p,q) model}
\description{
  The \code{cogarch.gmm} class is a class of the  \pkg{yuima} package that contains estimated parameters obtained by the function \code{\link{gmm}}.  
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{is an object of of \code{\link{yuima.cogarch-class}}.}
    \item{\code{objFun}:}{is an object of class \code{character} that indicates the objective function used in the minimization problem. See the documentation of the function \code{\link{gmm}} for more details.}
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
