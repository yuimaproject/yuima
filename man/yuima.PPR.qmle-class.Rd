\name{Class for Quasi Maximum Likelihood Estimation of Point Process Regression Models}
\docType{class}
\alias{yuima.PPR.qmle-class}
\alias{yuima.PPR.qmle,ANY-method}
\alias{qmle.PPR}
\alias{PPR.qmle}

\title{Class for Quasi Maximum Likelihood Estimation of Point Process Regression Models}
\description{ The \code{yuima.PPR.qmle} class is a class of the  \pkg{yuima} package that extends the \code{mle-class} of the \pkg{stats4} package.
}
\section{Slots}{
  \describe{
     \item{\code{call}:}{is an object of class \code{language}. }
    \item{\code{coef}:}{is an object of class \code{numeric} that contains estimated parameters.}
    \item{\code{fullcoef}:}{is an object of class \code{numeric} that contains estimated and fixed parameters.}
    \item{\code{vcov}:}{is an object of class \code{matrix}.}
    \item{\code{min}:}{is an object of class \code{numeric}.}
    \item{\code{minuslogl}:}{is an object of class \code{function}.}
    \item{\code{method}:}{is an object of class \code{character}.}
    \item{\code{model}:}{is an object of class \code{\link{yuima.PPR-class}}.}
  }
}
\section{Methods}{
  \describe{
    \item{Methods mle}{All methods for \code{mle-class} are available.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}

