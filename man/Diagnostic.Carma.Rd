\name{Diagnostic.Carma}
\alias{Diagnostic.Carma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostic Carma model}
\description{
This function verifies if the condition of stationarity is satisfied.}
\usage{
Diagnostic.Carma(carma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{carma}{ An object of class \code{yuima.qmle-class} where the slot \code{model} is a carma process.
}
}
\value{Logical variable. If \code{TRUE}, Carma is stationary.}
\author{
YUIMA TEAM
}


\examples{

mod1 <- setCarma(p = 2, q = 1, scale.par = "sig",
          Carma.var = "y")

param1 <- list(a1 = 1.39631, a2 = 0.05029, b0 = 1,
            b1 = 1, sig = 1)
samp1 <- setSampling(Terminal = 100, n = 200)

set.seed(123)

sim1 <- simulate(mod1, true.parameter = param1,
          sampling = samp1)

est1 <- qmle(sim1, start = param1)

Diagnostic.Carma(est1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Carma}% use one of  RShowDoc("KEYWORDS")
\keyword{Noise}% __ONLY ONE__ keyword per line
