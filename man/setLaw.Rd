\name{setLaw}
\alias{setLaw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random Variable Constructor }
\description{'\code{setLaw}' constructs an object of class \code{\link{yuima.law}} that contains user-defined random number generator, density, cdf, quantile function and characteristic function of a driving noise for a stochastic model.}
\usage{
setLaw(rng = function(n, ...) {
    NULL
}, density = function(x, ...) {
    NULL
}, cdf = function(q, ...) {
    NULL
}, quant = function(p, ...) {
    NULL
}, characteristic = function(u, ...) {
    NULL
}, time.var = "t", dim = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rng}{A user defined
function for the generation of the noise sample.}
  \item{density}{A user defined function for the computation of the noise density.}
  \item{cdf}{A user defined function for the computation of the noise cumulative distribution function.}
  \item{characteristic}{A user defined function for the computation of the characteristi function.}
  \item{quant}{A user defined function for the computation of the quantile.}
  \item{time.var}{A label of the time variable.}
  \item{dim}{Dimension of the noise.}
}
\details{The constructor produces an \code{\link{yuima.law}}-object  entirely specified  by the user. This object can be employed in all stochastic processes available in \code{yuima} as a driving noise. 

The minimal requirement for simulating the stochastic model is the specification of the \code{rng} function. 

The \code{density} function is used internally in the  \code{\link{qmleLevy}} function for the estimation of the Levy parameters.}
\value{
The constructor returns an object of class \code{\link{yuima.law}}}
\author{
YUIMA TEAM}
\note{
There may be missing information in the description. 
Please contribute with suggestions and fixings.

Contacts: Lorenzo Mercuri \email{lorenzo.mercuri@unimi.it}}

\references{
               Masuda, H., Mercuri, L. and Uehara, Y. (2022),
Noise Inference For Ergodic Levy Driven SDE, \emph{Electronic Journal of Statistics}, \bold{16}(1), 2432-2474. \doi{10.1214/22-EJS2006}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%
\examples{
\dontrun{
### The following example reports all steps 
### for the construction of an yuima.law
### object using the external R package 
### VarianceGamma available on CRAN

library(VarianceGamma)

## Definition of a yuima.law object 

# User defined rng
myrng <- function(n, eta, t){ 
          rvg(n, vgC = 0, sigma = sqrt(t), theta = 0, nu = 1/(eta*t))
}

# User defined density
mydens <- function(x, eta, t){
           dvg(x, vgC = 0, sigma = sqrt(t), theta = 0, nu = 1/(eta*t))
}

mylaw <- setLaw(rng = myrng, density = mydens, dim = 1)

}

}
