\name{setMap}
\alias{setMap}
\alias{Map of SDE}
\alias{Map of yuima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Map of a Stochastic Differential Equation
}
\description{
'\code{setMap}' is the constructor of an object of class \code{\link{yuima.Map}} that describes a map of a SDE
}
\usage{
setMap(func, yuima, out.var = "", nrow = 1, ncol = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{a matrix or a vector of strings that describe each component of the map.}
  \item{yuima}{ an object of class \code{\link{yuima.model}} that is the SDE.}
  \item{out.var}{ label for the output}
  \item{nrow}{ dimension of Map if \code{func} is a vector of string. }
  \item{ncol}{ dimension of output if \code{func} is a vector of string.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{The constructor returns an object of class \code{\link{yuima.Map}}.}
\references{Yuima Documentation}
\author{The YUIMA Project Team}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

\dontrun{
# Definition of a yuima model
mod <- setModel(drift=c("a1", "a2"),
  diffusion = matrix(c("s1","0","0","s2"),2,2),
  solve.variable = c("X","Y"))

# Definition of a map
my.Map <- matrix(c("(X+Y)","-X-Y",
  "a*exp(X-a1*t)","b*exp(Y-a2*t)"),
  nrow=2,ncol=2)

# Construction of yuima.Map

yuimaMap <- setMap(func = my.Map, yuima = mod,
  out.var = c("f11","f21","f12","f22"))

# Simulation of a Map

set.seed(123)
samp <- setSampling(0, 100,n = 1000)
mypar <- list(a=1, b=1, s1=0.1, s2=0.2, a1=0.1, a2=0.1)
sim1 <- simulate(object = yuimaMap, true.parameter = mypar,
  sampling = samp)

# plot

plot(sim1, ylab = yuimaMap@Output@param@out.var,
  main = "simulation Map", cex.main = 0.8)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
