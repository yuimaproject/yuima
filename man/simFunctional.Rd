\name{simFunctional}
\alias{simFunctional}
\alias{simFunctional,yuima-method}
\alias{Fnorm}
\alias{Fnorm,yuima-method}
\alias{F0}
\alias{F0,yuima-method}

\title{Calculate the value of functional}
\description{Calculate the value of functional associated with sde by Euler scheme.}
\usage{
simFunctional(yuima, expand.var="e")
Fnorm(yuima, expand.var="e")
F0(yuima, expand.var="e")
}
\arguments{
  \item{yuima}{a \code{yuima} object containing model, functional and data.}
  \item{expand.var}{default expand.var="e".}
}
\details{
Calculate the value of functional of interest. Fnorm returns normalized one, and F0 returns the value for the case small parameter epsilon = 0. In simFunctional and Fnorm, yuima MUST contains the 'data' slot (X in legacy version)
}
\value{
  \item{Fe}{a real value}
}
\author{YUIMA Project Team}
\note{
we need to fix this routine.
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
samp <- setSampling(Terminal = Terminal, n = division)
yuima <- setYuima(model = model,sampling = samp) 
yuima <- setFunctional( yuima, xinit=xinit, f=f,F=F,e=e)
# evaluate the functional value

yuima <- simulate(yuima,xinit=xinit,true.par=e)
Fe <- simFunctional(yuima)
Fe
Fenorm <- Fnorm(yuima)
Fenorm
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
