\name{cogarch.info-class}
\docType{class}
\alias{cogarch.info-class}

\title{Class for information about CoGarch(p,q)}
\description{
  The \code{cogarch.info-class} is a class of the  \pkg{yuima} package
}
\section{Slots}{
  \describe{
    \item{\code{p}:}{Number of autoregressive coefficients in the variance process. }
    \item{\code{q}:}{Number of moving average coefficients in the variance process. }
    \item{\code{ar.par}:}{Label of autoregressive coefficients. }
    \item{\code{ma.par}:}{Label of moving average coefficients. }
    \item{\code{loc.par}:}{Label of location coefficient in the variance process. }
    \item{\code{Cogarch.var}:}{Label of the observed process. }
    \item{\code{V.var}:}{Label of the variance process. }
    \item{\code{Latent.var}:}{Label of the latent process in the state representation of the variance.}
    \item{\code{XinExpr}:}{Logical variable. If \code{XinExpr=FALSE}, the starting condition of \code{Latent.var} is zero otherwise each component of Latent.var has a parameter as a starting point. }
    \item{\code{measure}:}{Levy measure for jump and quadratic part. }
    \item{\code{measure.type}:}{Type specification for Levy measure. }
  }
}

\note{
The \code{cogarch.info-class} object cannot be directly specified by the user
but it is built when the \code{\link{yuima.cogarch-class}} object is
constructed via \code{\link{setCogarch}}.
}

\author{The YUIMA Project Team}
\keyword{classes}
