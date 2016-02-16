\name{setMultiModel}
\alias{setMultiModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multidimensional Jump Diffusion Model
}
\description{
'\code{setMultiModel}' generalizes \code{\link{setModel}} since it allows the user to construct a multidimensional Jump Diffusion Model where the jump noise is a multidimensional Levy process.
}
\usage{
setMultiModel(drift = NULL, diffusion = NULL, hurst = 0.5, jump.coeff = NULL,
  intensity = NULL, df = NULL, measure.type = character(),
  state.variable = "x", jump.variable = "z", time.variable = "t",
  solve.variable, xinit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drift}{ see \code{\link{setModel}}.}
  \item{diffusion}{ see \code{\link{setModel}}.}
  \item{hurst}{ see \code{\link{setModel}}.}
  \item{jump.coeff}{see \code{\link{setModel}}.}
  \item{intensity}{ a string that indicates the label of intensity parameter. If \code{intensity = NULL} the Levy jumps are not CP type.}
  \item{df}{ a list containing a density defined by user}
  \item{measure.type}{ see \code{\link{setModel}}.}
  \item{state.variable}{see \code{\link{setModel}}.}
  \item{jump.variable}{see \code{\link{setModel}}.}
  \item{time.variable}{see \code{\link{setModel}}.}
  \item{solve.variable}{see \code{\link{setModel}}.}
  \item{xinit}{see \code{\link{setModel}}.}
}
\details{The user has to specify its proper function for multivariate density and random number generator}
\value{The function returns an object of class \code{\link{yuima.multimodel}}.}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{YUIMA Team}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
# We define the density function of the underlying Levy

dmyexp <- function(z, sig1, sig2, sig3){
  rep(0,3)
}

# We define the random number generator

rmyexp <- function(z, sig1, sig2, sig3){
  cbind(rnorm(z,0,sig1), rgamma(z,1,sig2), rnorm(z,0,sig3))
}

# Model Definition: in this case we consider only a multi
# compound poisson process with a common intensity as underlying
# noise

mod <- setMultiModel(drift = matrix(c("0","0","0"),3,1), diffusion = NULL,
  jump.coeff = matrix(c("1","0","0","0","1","-1","1","0","0"),3,3),
  intensity = c("lambda1"), df = list("dmyexp(z,sig1,sig2,sig3)"),
  jump.variable = c("z"), measure.type=c("CP"),
  solve.variable=c("X1","X2","X3"))

# Sample scheme

samp<-setSampling(0,100,n=1000)
param <- list(lambda1 = 1, sig1 = 0.1, sig2 = 0.1, sig3 = 0.1)

# Simulation

traj <- simulate(object = mod, sampling = samp,
  true.parameter = param)

# Plot

plot(traj, main = " driven noise. Multidimensional CP",
  cex.main = 0.8)

# We construct a multidimensional SDE driven by a multivariate
# levy process without CP components.

# Definition multivariate density

dmyexp1 <- function(z, sig1, sig2, sig3){
  rep(0,3)
}

# Definition of random number generator
# In this case user must define the delta parameter in order to
# control the effect of time interval in the simulation.

rmyexp1 <- function(z, sig1, sig2, sig3, delta){
  cbind(rexp(z,sig1*delta), rgamma(z,1*delta,sig2), rexp(z,sig3*delta))
}

# Model defintion

mod1 <- setMultiModel(drift=matrix(c("0.1*(0.01-X1)",
  "0.05*(1-X2)","0.1*(0.1-X3)"),3,1), diffusion=NULL,
   jump.coeff = matrix(c("0.01","0","0","0","0.01",
                          "0","0","0","0.01"),3,3),
   df=list("dmyexp1(z,sig1,sig2,sig3)"),
   jump.variable = c("z"), measure.type=c("code"),
   solve.variable=c("X1","X2","X3"),xinit=c("10","1.2","10"))

# Simulation sample paths

samp<-setSampling(0,100,n=1000)
param <- list(sig1 = 1, sig2 = 1, sig3 = 1)

# Simulation

set.seed(1)
traj1 <- simulate(object = mod1, sampling = samp,
  true.parameter = param)

# Plot

plot(traj1, main = "driven noise: multi Levy without CP",
  cex.main = 0.8)

# We construct a multidimensional SDE driven by a multivariate
# levy process.

# We consider a mixed situation where some
# noise are driven by a multivariate Compuond Poisson that
# shares a common intensity parameters.

### Multi Levy model

rmyexample2<-function(z,sig1,sig2,sig3, delta){
    if(missing(delta)){
      delta<-1
    }
    cbind(rexp(z,sig1*delta), rgamma(z,1*delta,sig2),
        rexp(z,sig3*delta), rep(1,z),
        rep(1,z))
}

dmyexample2<-function(z,sig1,sig2,sig3){
  rep(0,5)
}

# Definition Model

mod2 <- setMultiModel(drift=matrix(c("0.1*(0.01-X1)",
  "0.05*(1-X2)","0.1*(0.1-X3)", "0", "0"),5,1), diffusion=NULL,
  jump.coeff = matrix(c("0.01","0","0","0","0",
                        "0","0.01","0","0","0",
                        "0","0","0.01","0","0",
                        "0","0","0","0.01","0",
                        "0","0","0","0","0.01"),5,5),
  df = list("dmyexample2(z,sig1,sig2,sig3)"),
  intensity = c("lambda1"), jump.variable = c("z"),
  measure.type=c("code","code","code","CP","CP"),
  solve.variable=c("X1","X2","X3","X4","X5"),
  xinit=c("10","1.2","10","0","0"))

# Simulation scheme
samp <- setSampling(0, 100, n = 1000)
param <- list(sig1 = 1, sig2 = 1, sig3 = 1, lambda1 = 1)

# Simulation

set.seed(1)
traj2 <- simulate(object = mod2, sampling = samp,
  true.parameter = param)

plot(traj2, main = "driven noise: general multi Levy", cex.main = 0.8)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
