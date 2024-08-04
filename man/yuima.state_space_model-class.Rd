\name{yuima.state_space_model-class}
\docType{class}
\alias{yuima.state_space_model-class}
\alias{yuima.state_space_model}
\alias{initialize,yuima.state_space_model-method}
\alias{initialize,state_space_model.parameter-method}
\alias{limiting.gamma,yuima.state_space_model-method}
\alias{simulate,yuima.state_space_model-method}

\title{Class for the Mathematical Description of State Space Models}
\description{
  The \code{yuima.state_space_model} class is a class in the \pkg{yuima} package
  for the mathematical description and simulation of state space models.
}
\section{Slots}{
  \describe{
    \item{\code{drift}:}{An \R{} expression specifying the drift coefficient (a vector).}
    \item{\code{diffusion}:}{An \R{} expression specifying the diffusion coefficient (a matrix).}
    \item{\code{hurst}:}{The Hurst parameter of the Gaussian noise. If \code{hurst=0.5}, 
    the process is a Wiener process; otherwise, it is fractional Brownian motion 
    with the specified Hurst index. Can be set to \code{NA} for further specification.}
    \item{\code{jump.coeff}:}{A matrix of \R{} expressions for the jump component.}
    \item{\code{measure}:}{The Levy measure for jump variables.}
    \item{\code{measure.type}:}{The type specification for Levy measures.}
    \item{\code{state.variable}:}{A vector of names identifying the state variables 
    used in the drift and diffusion specifications.}
    \item{\code{parameter}:}{An object of class \code{\link{model.parameter-class}}, 
    representing the model parameters. 
    For more details, see the \code{\link{model.parameter-class}} documentation page.}
    \item{\code{jump.variable}:}{The variable for the jump coefficient.}
    \item{\code{time.variable}:}{The time variable.}
    \item{\code{noise.number}:}{The number of sources of noise, currently only for the Gaussian part.}
    \item{\code{equation.number}:}{The dimension of the stochastic differential equation.}
    \item{\code{dimension}:}{The dimensions of the parameter given in the \code{parameter} slot.}
    \item{\code{solve.variable}:}{The variable with respect to which the stochastic differential equation is solved.}
    \item{\code{xinit}:}{The initial value of the stochastic differential equation.}
    \item{\code{J.flag}:}{Indicates whether \code{jump.coeff} includes \code{jump.variable}.}
    \item{\code{is.observed}:}{Indicates whether each state variable is observed (i.e., data is given) or not.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
