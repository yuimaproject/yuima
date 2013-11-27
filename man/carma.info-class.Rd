\name{carma.info-class}
\docType{class}
\alias{carma.info-class}

\title{Class for information about CARMA(p,q) model}
\description{
  The \code{carma.info-class} is a class of the  \pkg{yuima} package.  
}
\section{Slots}{
  \describe{
    \item{\code{p}:}{Number of autoregressive coefficients.}
    \item{\code{q}:}{Number of moving average coefficients.}
    \item{\code{loc.par}:}{Label of location coefficient.}
    \item{\code{scale.par}:}{Label of scale coefficient.}
    \item{\code{ar.par}:}{Label of autoregressive coefficients.}
    \item{\code{ma.par}:}{Label of moving average coefficients.}
    \item{\code{lin.par}:}{Label of linear coefficients.}
    \item{\code{Carma.var}:}{Label of the observed process.}
    \item{\code{Latent.var}:}{Label of the unobserved process.}
    \item{\code{XinExpr}:}{Logical variable. If \code{XinExpr=FALSE}, the starting condition of \code{Latent.var} is zero otherwise each component of Latent.var has a parameter as a starting point. }
  }
}
\details{
The \code{carma.info-class} object cannot be directly specified by the user
but it is constructed when the \code{\link{yuima.carma-class}} object is
constructed via \code{\link{setCarma}}.
}
\author{The YUIMA Project Team}
\keyword{classes}
