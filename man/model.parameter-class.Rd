\name{model.parameter-class}
\docType{class}
\alias{model.parameter-class}

\title{Class for the parameter description of stochastic differential equations}
\description{
  The \code{model.parameter-class} is a class of the  \pkg{yuima} package.  
}
\section{Slots}{
  \describe{
    \item{\code{drift}:}{A vector of names belonging to 
		the drift coefficient.}
    \item{\code{diffusion}:}{A vector of names of parameters belonging to 
		the diffusion coefficient.}
    \item{\code{jump}:}{A vector of names of parameters belonging to 
		the jump coefficient.}
    \item{\code{measure}:}{A vector of names of parameters belonging to 
		the Levy measure.}
    \item{\code{xinit}:}{A vector of names of parameters belonging to
		the initial condition.}
    \item{\code{all}:}{A vector of names of all the parameters found in the
	components of the model.}
    \item{\code{common}:}{A vector of names of the parameters in common among
	 drift, diffusion, jump and measure term.}
  }
}
\details{
The \code{model.parameter-class} object cannot be directly specified by the user
but it is constructed when the \code{\link{yuima.model-class}} object is
constructed via \code{\link{setModel}}.
All the terms which are not in the list of \emph{solution, state, time, jump} variables
are considered as parameters. These parameters are
identified in the different components of the model (drift, diffusion and
jump part).
This information is later used to draw inference jointly or separately for
the different parameters depending on the model in hands. 
}
\author{The YUIMA Project Team}
\keyword{classes}
