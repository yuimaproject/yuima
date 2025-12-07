\name{carmaHawkes.info-class}
\docType{class}
\alias{carmaHawkes.info-class}

\title{Class for Information on the Hawkes Process with a CARMA(p,q) Intensity}

\description{
  The \code{carmaHawkes.info-class} is a class of the  \pkg{yuima} package.  
}
\section{Slots}{
  \describe{
    \item{\code{p}:}{Number of autoregressive coefficients.}
    \item{\code{q}:}{Number of moving average coefficients.}
    \item{\code{Counting.Process}:}{Label of Counting process.}
    \item{\code{base.Int}:}{Label of baseline Intensity parameter.}
    \item{\code{ar.par}:}{Label of autoregressive coefficients.}
    \item{\code{ma.par}:}{Label of moving average coefficients.}
    \item{\code{Intensity.var}:}{Label of the Intensity process.}
    \item{\code{Latent.var}:}{Label of the unobserved process.}
    \item{\code{XinExpr}:}{Logical variable. If \code{XinExpr=FALSE}, the starting condition of \code{Latent.var} is zero otherwise each component of Latent.var has a parameter as a starting point. }
    \item{\code{Type.Jump}:}{Logical variable. If \code{XinExpr=TRUE}, the jump size is deterministic}
  }
}

\details{
The \code{carmaHawkes.info-class} object cannot be directly specified by the user
but it is constructed when the \code{\link{yuima.carmaHawkes-class}} object is
constructed via \code{\link{setCarmaHawkes}}.
}
\author{The YUIMA Project Team

Contacts: Lorenzo Mercuri \email{lorenzo.mercuri@unimi.it}
}
\keyword{classes}


