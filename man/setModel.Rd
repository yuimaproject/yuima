\name{setModel}
\alias{setModel}

\title{
  Basic description of stochastic differential equations (SDE)
}

\description{
  'setModel' is used to give a mathematical description of 
  diffusion processes (we will incorporate more general SDEs 
  including a jump component in time). All functions in the yuima 
  package should assume to get as much information as possible 
  from 'setModel', without replicating the same code in each place. 
  If there are missing pieces of information, some default values 
  are used.
}

\usage{
setModel(drift=NULL, diffusion=NULL, jump.coeff=character(),
measure=list(), measure.type=character(), state.variable="x",
jump.variable="z", time.variable="t", solve.variable, xinit)
}

\arguments{
  \item{drift}{a vector of \code{expression}s (the default value is 0).}
  \item{diffusion}{a matrix of \code{expression}s (the default value is 1).}
  \item{jump.coeff}{a matrix of \code{expression}x.}
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
  \item{xinit}{a vector of numbers identifying the initial figures of data.}
}

\details{
You should look at the vignette and the examples.

	An object of \code{yuima.model} contains several ``slots'' below: 
	to see inside its structure we use the \R{} command \code{str}.
	
	\code{drift} is an \R{} expression which specifies 
	the drift coefficient (a vector).
	
	\code{diffusion} is an \R{} expression which specifies 
	the diffusion coefficient (a matrix).
	
	\code{JUMP.COEFF} is coefficients of the jump term.
	
	\code{MEASURE} is a measure of the Levy process.

	\code{MEASURE.TYPE} is a type specification of measure,
	such as \code{CP}, \code{code} or \code{density}.

	\code{parameter} which is a short name for ``parameters'' which is 
	a list of objects:
	
	\code{all} contains the names of all the parameters found in the 
	diffusion and drift coefficients;
	
	\code{common} contains the names of the parameters that both the drift 
	and the diffusion coefficients contain;

	\code{diffusion} contains the parameters belonging only to the diffusion coefficient;
	
	\code{drift} contains the parameters belonging only to the drift
	coefficient.

	\code{jump} contains the parameters belonging only to the jump
	coefficient.

	\code{measure} contains the parameters belonging only to the measure
	coefficient.	

	\code{solve.variable} contains a vector of variable names,
	each element corresponds to the name of the solution variable (left-hand-side) of 
	each equation in the model, in the corresponding order;	  
	  
	next we have the names of the \code{state.variable}, \code{jump.variable} and \code{time.variable} variables. 
	By default, they are assumed to be $x$ and $t$ but the user can freely choose them.
	The \code{yuima.model} function assumes that the user either use default
	names for \code{state.variable}, \code{jump.variable} and \code{time.variable} variables or specify his own names. 
	All the rest of the symbols are considered parameters and distributed accordingly 
	in the \code{parameter} slot.

	\code{noise.number} denotes the number of sources of noise; currently, this is 
	the number of independent Brownian motions, but jump noise will be included sometime 
	in the future.

	\code{equation.number} denotes the dimension of the SDEs.

	\code{dimension} reports the dimensions of the parameter given in \code{parameter}.

	\code{xinit} denotes the initial numbers of SDEs.
}

\value{
  \item{model}{an object of class 'yuima.model'}
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


