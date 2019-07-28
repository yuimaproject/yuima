\name{setPPR}
\alias{setPPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Point Process}
\description{Constructor of a Point Process Regression Model}
\usage{
setPPR(yuima, counting.var = "N", gFun, Kernel,
  var.dx = "s", var.dt = "s", lambda.var = "lambda",
  lower.var = "0", upper.var = "t", nrow = 1, ncol = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{an object of \code{\link{yuima.model-class}} that describes the mathematical features of counting and covariates processes \code{Y[t]=(X[t],N[t])}.}
  \item{counting.var}{a label denoting the name of the counting process.}
  \item{gFun}{a vector string that is the mathematical expression of the vector function \code{g(t,Y[t-],theta)} in the intensity process.}
  \item{Kernel}{a matrix string that is the kernel \code{kappa(t-s,Y[s],theta)} in the definition of the intensity process.}
  \item{var.dx}{a string denoting the integration variable in the intensity process.}
  \item{var.dt}{a string denoting the integration time variable in the intensity process.}
  \item{lambda.var}{name of the intensity process.}
  \item{lower.var}{Lower bound of the support for the integral in the definition of the intensity process.}
  \item{upper.var}{Upper bound of the support for the integral in the definition of the intensity process.}
  \item{nrow}{number of rows in the kernel.}
  \item{ncol}{number of columns in the kernel.}
}
\value{An object of \code{\link{yuima.PPR}}}
\references{Insert Here References}
\author{The YUIMA Project Team

Contacts: Lorenzo Mercuri \email{lorenzo.mercuri@unimi.it}}
\note{There may be missing information in the model description. Please contribute with suggestions and fixings.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
## Hawkes process with power law kernel

# I. Law Definition:
my.rHwk2 <- function(n){
  as.matrix(rep(1,n))
  }
Law.Hwk2 <- setLaw(rng = my.rHwk2, dim = 1)

# II. Definition of the counting process N_t
mod.Hwk2 <- setModel(drift = c("0"), diffusion = matrix("0",1,1),
  jump.coeff = matrix(c("1"),1,1), measure = list(df = Law.Hwk2),
  measure.type = "code", solve.variable = c("N"),
  xinit=c("0"))

# III. Definition of g() and kappa()
g.Hwk2 <- "mu"
Kern.Hwk2 <- "alpha/(1+(t-s))^beta"

# IV. Construction of an yuima.PPR object
PPR.Hwk2 <- setPPR(yuima = mod.Hwk2, gFun=g.Hwk2,
  Kernel = as.matrix(Kern.Hwk2),var.dx = "N")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{PPR model}% use one of  RShowDoc("KEYWORDS")
\keyword{PPR Constructor}% __ONLY ONE__ keyword per line