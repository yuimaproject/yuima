\name{setModel}
\alias{setModel}

\title{
  Basic description of stochastic differential equations (SDE)
}

\description{
  'setModel' is used to give a mathematical description of 
  diffusion processes. All functions relying on the \pkg{yuima} 
  package should assume to get as much information as possible 
  from the different slots of the \code{\link{yuima-class}} structure
  without replicating the same code twice. 
  If there are missing pieces of information, some default values 
  can be assumed.
}

\usage{
setModel(drift = NULL, diffusion = NULL, hurst = 0.5, jump.coeff = character(),
measure = list(), measure.type = character(), state.variable = "x",
jump.variable = "z", time.variable = "t", solve.variable, xinit)
}

\arguments{
  \item{drift}{a vector of \code{expression}s (the default value is 0 when 
    \code{drift=NULL}).}
  \item{diffusion}{a matrix of \code{expression}s (the default value is 1 when
   \code{diffusion=NULL}).}
  \item{hurst}{the Hurst parameter of the gaussian noise. If \code{h=0.5}, the default,
  the process is Wiener otherwise it is fractional Brownian motion with that precise value
   of the Hurst index.}
  \item{jump.coeff}{a matrix of \code{expression}s for the jump component.}
  \item{measure}{Levy measure for jump variables.}
  \item{measure.type}{Type specification for Levy measures.}
  \item{state.variable}{a vector of names identifying the names used to 
    denote the state variable in the drift and diffusion
    specifications.}
  \item{jump.variable}{a vector of names identifying the names used to denote
    the jump variable in the jump specifications.}
  \item{time.variable}{a vector of names identifying the names used to 
    denote the time variable in the drift and diffusion specifications.}  
  \item{solve.variable}{a vector of names identifying the names used to 
    denote the left-hand-side of each equation in the model; this equals 
    \code{state.variable} as long as we have no exogenous variable other than 
    statistical parameters in the coefficients (\code{drift} and
    \code{diffusion}).}  
  \item{xinit}{a vector of numbers identifying the initial value of model.}
}

\details{
Please refer to the vignettes and the examples or to the \pkg{yuimadocs}
package.

An object of \code{\link{yuima.model-class}} contains several slots:

\describe{		
\item{\code{drift}:}{is an \R{} expression which specifies the drift coefficient 
(a vector).}
\item{\code{diffusion}:}{is an \R{} expression which specifies the diffusion 
coefficient (a matrix).}
\item{\code{jump.coeff}:}{is coefficients of the jump term.}
\item{\code{measure}:}{is a measure of the Levy process}
\item{\code{measure.type}:}{is a type specification of measure, such as 
 \code{CP}, \code{code} or \code{density}. See below.}
\item{\code{parameter}:}{which is a short name for ``parameters'', is an object 
 of class \code{\link{model.parameter-class}} which is a list of vectors of 
 names of parameter belonging to the single components of the model (drift, 
 diffusion, jump and measure), the names of common parameters and the names 
 of all parameters. For more details see \code{\link{model.parameter-class}} 
 documentation page.}
\item{\code{solve.variable}:}{is a vector of variable names, each element 
 corresponds to the name of the solution variable (left-hand-side) of each 
 equation in the model, in the corresponding order.}	  
\item{\code{state.variable}:}{used to identify the state variables in the \R{} 
 expression. By default, it is assumed to be \code{x}, otherwise the user can 
 choose any other model specification.}
\item{\code{jump.variable}:}{the variable for the jump coefficient. By default,
 it is assumed to be \code{z}, otherwise the user can choose any other model
 specification.}
\item{\code{time}:}{the time variable. By default, it is assumed to be \code{t}, 
 otherwise the user can choose any other model specification.}
\item{\code{solve.variable}:}{used to identify the solution variables in the
 \R{} expression, i.e. the variable with respect to which the stochastic 
 differential equation has to be solved. By default, it is assumed to be 
 \code{x}, otherwise the user can choose any other model specification.}
\item{\code{noise.number}:}{denotes the number of sources of noise. Currently 
 only for the Gaussian part.}
\item{\code{equation.number}:}{denotes the dimension of the stochastic 
 differential equation.}
\item{\code{dimension}:}{the dimensions of the parameter given in the 
 \code{parameter} slot.}
\item{\code{xinit}:}{denotes the initial value of the stochastic differential 
 equation.}
}

The \code{\link{yuima.model-class}} structure assumes that the user either uses the default
names for \code{state.variable}, \code{jump.variable}, \code{solution.variable} and 
\code{time.variable} or specify his own names. 
All the rest of the terms in the \R{} expressions are considered as parameters 
and identified accordingly in the \code{parameter} slot.
}
\value{
  \item{model}{an object of class \code{\link{yuima.model-class}}.}
}

\author{The YUIMA Project Team}
\note{There may be missing information in the model description. 
Please contribute with suggestions and fixings.
}

\examples{
# Ex 1. (One-dimensional homogeneous diffusion)
# To describe
# dXt = -3*Xt*dt + (1/(1+Xt^2))dWt,
# we set
mod1 <- setModel(drift = "-3*x", diffusion = "1/(1+x^{2})", solve.variable = c("x"))
# We may omit the solve.variable; then the default variable x is used
mod1 <- setModel(drift = "-3*x", diffusion = "1/(1+x^{2})")
# Look at the model structure by
str(mod1)

# Ex 2. (Two-dimensional homogeneous diffusion)
# To describe
# dX1t = -3*X1t*dt         + dW1t,
# dX2t = -(X1t + 2*X2t)*dt + 0.5*dW1t + dW2t,
# we set the drift coefficient
a <- c("-3*x1","-x1-2*x2")
# and also the diffusion coefficient
b <- matrix(c("1","0.5","0","1"),2,2)
# Then set
mod2 <- setModel(drift = a, diffusion = b, solve.variable = c("x1","x2"))
# Look at the model structure by
str(mod2)
# The noise.number is automatically determined by inputting the diffusion matrix expression.
# If the dimensions of the drift differs from the number of the rows of the diffusion, 
# the error message is returned.
}
\keyword{ts}


