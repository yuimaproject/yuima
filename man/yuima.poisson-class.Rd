\name{yuima.poisson-class}
\docType{class}
\alias{yuima.poisson-class}
\alias{initialize,yuima.poisson-method}

\title{Class for the mathematical description of Compound Poisson processes}
\description{
  The \code{yuima.poisson} class is a class of the  \pkg{yuima} package that extends the \code{\link{yuima.model-class}}.
}
\section{Slots}{
  \describe{
    \item{\code{drift}:}{always \code{expression((0))}.}
    \item{\code{diffusion}:}{a list of \code{expression((0))}.}
    \item{\code{hurst}:}{always \code{h=0.5}, but ignored for this model.}
    \item{\code{jump.coeff}:}{set according to \code{scale} in \code{\link{setPoisson}}.}
    \item{\code{measure}:}{a list containting the intensity measure and the jump distribution.}
    \item{\code{measure.type}:}{always \code{"CP"}.}
    \item{state.variable}{a vector of names identifying the names used to 
     denote the state variable in the drift and diffusion specifications.}
    \item{\code{parameter}:}{which is a short name for ``parameters'', is an 
	 object of class \code{\link{model.parameter-class}}. For more details see 
	 \code{\link{model.parameter-class}} documentation page.}
    \item{\code{state.variable}:}{identifies the state variables in the \R{} 
     expression.}
    \item{\code{jump.variable}:}{identifies the variable for the jump 
	 coefficient.}
    \item{\code{time.variable}:}{the time variable.}
    \item{\code{noise.number}:}{denotes the number of sources of noise.}
    \item{\code{equation.number}:}{denotes the dimension of the stochastic 
     differential equation.}
    \item{\code{dimension}:}{the dimensions of the parameter given in the 
     \code{parameter} slot.}
    \item{\code{solve.variable}:}{identifies the variable with respect to which 
	 the stochastic differential equation has to be solved.}
    \item{\code{xinit}:}{contains the initial value of the stochastic 
	 differential equation.}
    \item{\code{J.flag}:}{wheather jump.coeff include jump.variable.}
  }
}
\section{Methods}{
  \describe{
    \item{simulate}{simulation method. For more information see
	  \code{\link{simulate}}.}
%    \item{toLatex}{This method converts an object of \code{yuima.poisson-class} to character vectors with LaTeX markup.}
    \item{qmle}{Quasi maximum likelihood estimation procedure. For more information see \code{\link{qmle}}. }
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
