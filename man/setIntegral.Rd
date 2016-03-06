\name{setIntegral}
\alias{setIntegral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Integral of Stochastic Differential Equation}
\description{
'\code{setIntegral}' is the constructor of an object of class \link{yuima.Integral}
}
\usage{
setIntegral(yuima, integrand, var.dx, lower.var, upper.var,
 out.var = "", nrow = 1, ncol = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{
an object of class \code{\link{yuima.model}} that is the SDE.
}
  \item{integrand}{
A matrix or a vector of strings that describe each component of the integrand.
}
  \item{var.dx}{ A label that indicates the variable of integration }
  \item{lower.var}{
 A label that indicates the lower variable in the support of integration, by default \code{lower.var = 0}.}
  \item{upper.var}{
A label that indicates the upper variable in the support of integration, by default \code{upper.var = t}.}
  \item{out.var}{
Label for the output
}
  \item{nrow}{
Dimension of output if \code{integrand} is a vector of string.
}
  \item{ncol}{Dimension of output if \code{integrand} is a vector of string.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{The constructor returns an object of class \code{\link{yuima.Integral}}.}
\references{Yuima Documentation}
\author{The YUIMA Project Team}

\examples{
\dontrun{
# Definition Model

Mod1<-setModel(drift=c("a1"), diffusion = matrix(c("s1"),1,1),
  solve.variable = c("X"), time.variable = "s")

# In this example we define an integral of SDE such as
# \[
# I=\int^{t}_{0} b*exp(-a*(t-s))*(X_s-a1*s)dX_s
# \]

integ <- matrix("b*exp(-a*(t-s))*(X-a1*s)",1,1)

Integral <- setIntegral(yuima = Mod1,integrand = integ,
  var.dx = "X", lower.var = "0", upper.var = "t",
  out.var = "", nrow =1 ,ncol=1)

# Structure of slots

is(Integral)
# Function h in the above definition
Integral@Integral@Integrand@IntegrandList
# Dimension of Intgrand
Integral@Integral@Integrand@dimIntegrand

# all parameters are $\left(b,a,a1,s1\right)$
Integral@Integral@param.Integral@allparam

# the parameters in the integrand are $\left(b,a,a1\right)$ \newline
Integral@Integral@param.Integral@Integrandparam

# common parameters are $a1$
Integral@Integral@param.Integral@common

# integral variable dX_s
Integral@Integral@variable.Integral@var.dx
Integral@Integral@variable.Integral@var.time

# lower and upper vars
Integral@Integral@variable.Integral@lower.var
Integral@Integral@variable.Integral@upper.var

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
