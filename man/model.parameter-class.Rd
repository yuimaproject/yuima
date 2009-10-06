\name{model.parameter-class}
\docType{class}
\alias{model.parameter-class}

\title{Classes for the parameter description of stochastic differential equations}
\description{
  The \code{model.parameter-class} is a class of the  \pkg{yuima} package.  
}
\section{Slots}{
  \describe{
    \item{\code{drift}:}{Parameters in drift coefficients in sde.}
    \item{\code{diffusion}:}{Parameters in diffusion coefficients in sde.}
    \item{\code{jump}:}{Parameters in jump coefficients in sde.}
    \item{\code{measure}:}{Parameters in measure coefficients in sde.}
    \item{\code{all}:}{All parameters in sde.}
    \item{\code{common}:}{Common parameters in drift, diffusion, jump
	 and measure term.}
  }
}
\author{YUIMA Project}
\keyword{classes}
