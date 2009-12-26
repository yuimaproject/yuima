\name{setFunctional}
\alias{setFunctional}
\alias{setFunctional,yuima.model-method}
\alias{setFunctional,yuima-method}
\alias{getF}
\alias{getf}
\alias{getxinit}
\alias{gete}

\title{Description of a functional associated with a perturbed stochastic differential equation}
\description{This function is used to give a description
of the stochastic differential equation. The functional represent the price of the option in financial economics, for example.
}
\usage{
setFunctional(model, F, f, xinit,e)
}
\arguments{
  \item{model}{yuima or yuima.model object.}
  \item{F}{function of $X_t$ and $epsilon$}
  \item{f}{list of functions of $X_t$ and $epsilon$}
  \item{xinit}{initial values of state variable.}
  \item{e}{epsilon parameter}
}
\details{
You should look at the vignette and examples. 

	The object foi contains several ``slots''. To see inside its 
	structure we use the \R{} command str.
	f and Fare \R{} (list of) expressions which contains the functional of interest specification. 
	e is a small parameter on which we conduct asymptotic expansion of the functional.
}
\value{
  \item{yuima}{an object of class 'yuima' containing object of class 'functional'. If yuima object was given as 'model' argument, the result is just added and the other slots of the object are maintained.}
}
\author{The YUIMA Project Team}
\note{There may be missing information in the model description.
Please contribute with suggestions and fixings.
}
\examples{
set.seed(123)
# to the Black-Scholes economy:
# dXt^e = Xt^e * dt + e * Xt^e * dWt
diff.matrix <- matrix( c("x*e"), 1,1)
model <- setModel(drift = c("x"), diffusion = diff.matrix)
# call option is evaluated by averating
# max{ (1/T)*int_0^T Xt^e dt, 0}, the first argument is the functional of interest:
Terminal <- 1
xinit <- c(1)
f <- list( c(expression(x/Terminal)), c(expression(0)))
F <- 0
division <- 1000
e <- .3
yuima <- setYuima(model = model,sampling = setSampling(Terminal = Terminal, n = division))
yuima <- setFunctional( model = yuima, xinit=xinit, f=f,F=F,e=e)
# look at the model structure
str(yuima@functional)
}
\keyword{ts}


