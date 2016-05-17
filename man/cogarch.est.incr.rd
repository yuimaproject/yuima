\name{cogarch.est.incr-class}
\docType{class}
\alias{cogarch.est.incr-class}
\alias{plot,cogarch.est.incr,ANY-method}
\alias{est.cogarch.incr-class}
\alias{cogarch.est.incr-class}
\alias{simulate,cogarch.est.incr-method}
%%\alias{setSampling,yuima.carma-method}

\title{Class for Estimation of COGARCH(p,q) model with underlying increments}
\description{
  The \code{cogarch.est.incr} class is a class of the  \pkg{yuima} package that extends the \code{\link{cogarch.est-class}} and is filled by the function \code{\link{gmm}} or \code{\link{qmle}}.
}
\section{Slots}{
  \describe{
    \item{\code{Incr.Lev}:}{is an object of class \code{\link{zoo}} that contains the estimated increments of the noise obtained using \code{\link{cogarchNoise}}.}
    \item{\code{yuima}:}{is an object of of \code{\link{yuima-class}}.}
    \item{\code{logL.Incr}:}{is an object of class \code{numeric} that contains the value of the log-likelihood for estimated Levy increments.}
    \item{\code{objFun}:}{is an object of class \code{character} that indicates the objective function used in the minimization problem. See the documentation of the function \code{\link{gmm}} or \code{\link{qmle}}  for more details.}
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
   \item{simulate}{simulation method. For more information see \code{\link{simulate}}.}
    \item{plot}{Plot method for estimated increment of the noise.}
    \item{Methods mle}{All methods for \code{mle-class} are available.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
