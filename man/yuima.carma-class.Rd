\name{yuima.carma-class}
\docType{class}
\alias{yuima.carma-class}
\alias{initialize,yuima.carma-method}
\alias{initialize,carma.info-method}
\alias{limiting.gamma,yuima.carma-method}
\alias{simulate,yuima.carma-method}
%%\alias{setSampling,yuima.carma-method}

\title{Class for the mathematical description of CARMA(p,q) model}
\description{
  The \code{yuima.carma} class is a class of the  \pkg{yuima} package that extends the \code{\link{yuima.model-class}}.  
}
\section{Slots}{
  \describe{
    \item{\code{info}:}{is an \code{\link{carma.info-class}} object that describes the structure of the CARMA(p,q) model.}
    \item{\code{drift}:}{is an \R{} expression which specifies the drift 
	 coefficient (a vector).}
    \item{\code{diffusion}:}{is an \R{} expression which specifies the diffusion 
     coefficient (a matrix).}
    \item{\code{hurst}:}{the Hurst parameter of the gaussian noise. If 
	 \code{h=0.5}, the process is Wiener otherwise it is fractional Brownian 
	 motion with that precise value of the Hurst index. Can be set to \code{NA} for further specification.}
    \item{\code{jump.coeff}:}{a vector of \code{expression}s for the jump 
	 component.}
    \item{\code{measure}:}{Levy measure for jump variables.}
    \item{\code{measure.type}:}{Type specification for Levy measures.}
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
    \item{\code{noise.number}:}{denotes the number of sources of noise. 
	 Currently only for the Gaussian part.}
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
    \item{toLatex}{This method converts an object of \code{yuima.carma-class} to character vectors with LaTeX markup.}
    \item{CarmaNoise}{Recovering underlying Levy. For more information see \code{\link{CarmaNoise}}. }
    \item{qmle}{Quasi maximum likelihood estimation procedure. For more information see \code{\link{qmle}}. }
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
