\name{setModel}
\alias{setModel}

\title{
  Basic description of stochastic differential equations (SDE)
}

\description{
  'setModel' gives a description of 
  stochastic differential equation with or without jumps of the following
  form:

  \code{dXt = a(t,Xt, alpha)dt + b(t,Xt,beta)dWt + c(t,Xt,gamma)dZt,  X0=x0}  

  All functions relying on the \pkg{yuima} 
  package will get as much information as possible 
  from the different slots of the \code{\link{yuima-class}} structure
  without replicating the same code twice. 
  If there are missing pieces of information, some default values 
  can be assumed.
}

\usage{
setModel(drift = NULL, diffusion = NULL, hurst = 0.5, jump.coeff = NULL,
measure = list(), measure.type = character(), state.variable = "x",
jump.variable = "z", time.variable = "t", solve.variable, xinit)
}

\arguments{
  \item{drift}{a vector of \code{expression}s (the default value is 0 when 
    \code{drift=NULL}).}
  \item{diffusion}{a matrix of \code{expression}s (the default value is 0 when
   \code{diffusion=NULL}).}
  \item{hurst}{the Hurst parameter of the gaussian noise. If \code{h=0.5}, the default,
  the process is Wiener otherwise it is fractional Brownian motion with that precise value
   of the Hurst index. Can be set to \code{NA} for further specification.}
  \item{jump.coeff}{a matrix of \code{expression}s for the jump component.}
  \item{measure}{Levy measure for jump variables.}
  \item{measure.type}{type specification for Levy measures.}
  \item{state.variable}{a vector of names of the state variables in the drift 
	and diffusion coefficients.}
  \item{jump.variable}{a vector of names of the jump variables in the jump 
    coefficient.}
  \item{time.variable}{the name of the time variable.}  
  \item{solve.variable}{a vector of names of the variables in the left-hand-side 
    of the equations in the model; \code{solve.variable} equals 
    \code{state.variable} as long as we have no exogenous variable other than 
    statistical parameters in the coefficients (\code{drift} and
    \code{diffusion}).}  
  \item{xinit}{a vector of numbers identifying the initial value of the 
  \code{solve.variable}.}
}

\details{
Please refer to the vignettes and the examples or to the \pkg{yuimadocs}
package.

An object of \code{\link{yuima.model-class}} contains several slots:

\describe{		
\item{\code{drift}:}{an \R expression which specifies the drift coefficient 
(a vector).}
\item{\code{diffusion}:}{an \R expression which specifies the diffusion 
coefficient (a matrix).}
\item{\code{jump.coeff}:}{coefficient of the jump term.}
\item{\code{measure}:}{the Levy measure of the driving Levy process.}
\item{\code{measure.type}:}{specifies the type of the measure, such as 
 \code{CP}, \code{code} or \code{density}. See below.}
\item{\code{parameter}:}{a short name for ``parameters''. It  is an object 
 of  \code{\link{model.parameter-class}} which is a list of vectors of 
 names of parameters belonging to the single components of the model (drift, 
 diffusion, jump and measure), the names of common parameters and the names 
 of all parameters. For more details see \code{\link{model.parameter-class}} 
 documentation page.}
\item{\code{solve.variable}:}{a vector of variable names, each element 
 corresponds to the name of the solution variable (left-hand-side) of each 
 equation in the model, in the corresponding order.}	  
\item{\code{state.variable}:}{identifies the state variables in the \R 
 expression. By default, it is assumed to be \code{x}.}
\item{\code{jump.variable}:}{the variable for the jump coefficient. By default,
 it is assumed to be \code{z}.}
\item{\code{time}:}{the time variable. By default, it is assumed to be \code{t}.}
\item{\code{solve.variable}:}{used to identify the solution variables in the
 \R expression, i.e. the variable with respect to which the stochastic 
 differential equation has to be solved. By default, it is assumed to be 
 \code{x}, otherwise the user can choose any other model specification.}
\item{\code{noise.number}:}{denotes the number of sources of noise. Currently 
 only for the Gaussian part.}
\item{\code{equation.number}:}{denotes the dimension of the stochastic 
 differential equation.}
\item{\code{dimension}:}{the dimensions of the parameters in the 
 \code{parameter} slot.}
\item{\code{xinit}:}{denotes the initial value of the stochastic differential 
 equation.}
}

The \code{\link{yuima.model-class}} structure assumes that the user either uses the default
names for \code{state.variable}, \code{jump.variable}, \code{solution.variable} and 
\code{time.variable} or specifies his/her own names. 
All the rest of the terms in the \R expressions are considered as parameters 
and identified accordingly in the \code{parameter} slot.
}
\value{
  \item{model}{an object of  \code{\link{yuima.model-class}}.}
}

\author{The YUIMA Project Team}
\note{There may be missing information in the model description. 
Please contribute with suggestions and fixings.
}

\examples{
# Ex 1. (One-dimensional diffusion process)
# To describe
# dXt = -3*Xt*dt + (1/(1+Xt^2+t))dWt,
# we set
mod1 <- setModel(drift = "-3*x", diffusion = "1/(1+x^2+t)", solve.variable = c("x"))
# We may omit the solve.variable; then the default variable x is used
mod1 <- setModel(drift = "-3*x", diffusion = "1/(1+x^2+t)")
# Look at the model structure by
str(mod1)

# Ex 2. (Two-dimensional diffusion process with three factors)
# To describe
# dX1t = -3*X1t*dt + dW1t +X2t*dW3t,
# dX2t = -(X1t + 2*X2t)*dt + X1t*dW1t + 3*dW2t,
# we set the drift coefficient
a <- c("-3*x1","-x1-2*x2")
# and also the diffusion coefficient
b <- matrix(c("1","x1","0","3","x2","0"),2,3)
# Then set
mod2 <- setModel(drift = a, diffusion = b, solve.variable = c("x1","x2"))
# Look at the model structure by
str(mod2)
# The noise.number is automatically determined by inputting the diffusion matrix expression.
# If the dimensions of the drift differs from the number of the rows of the diffusion, 
# the error message is returned.

# Ex 3. (Process with jumps (compound Poisson process))
# To describe
# dXt = -theta*Xt*dt+sigma*dZt
mod3 <- setModel(drift=c("-theta*x"), diffusion="sigma",
jump.coeff="1", measure=list(intensity="1", df=list("dnorm(z, 0, 1)")),
measure.type="CP", solve.variable="x")
# Look at the model structure by
str(mod3)

# Ex 4. (Process with jumps (stable process))
# To describe
# dXt = -theta*Xt*dt+sigma*dZt
mod4 <- setModel(drift=c("-theta*x"), diffusion="sigma",
jump.coeff="1", measure.type="code",measure=list(df="rstable(z,1,0,1,0)"), solve.variable="x")
# Look at the model structure by
str(mod4)
# See rng about other candidate of Levy noises.

# Ex 5. (Two-dimensional stochastic differenatial equation with Levy noise)
# To describe
# dX1t = (1 - X1t - X2t)*dt+dZ1t 
# dX2t = (0.5 - X1t - X2t)*dt+dZ2t
beta<-c(.5,.5)
mu<-c(0,0)
Lambda<-matrix(c(1,0,0,1),2,2)
mod5 <- setModel(drift=c("1 - x1-x2",".5 - x1-x2"), 
solve.variable=c("x1","x2"), jump.coeff=Lambda, measure.type="code",
measure=list(df="rNIG(z, alpha, beta, delta0, mu, Lambda)"))
# Look at the model structure by
str(mod5)

# Ex 6. (Process with fractional Gaussian noise)
# dYt = 3*Yt*dt + dWt^h
mod6 <- setModel(drift="3*y", diffusion=1, hurst=0.3, solve.variable=c("y"))
# Look at the model structure by
str(mod6)
}
\keyword{ts}

