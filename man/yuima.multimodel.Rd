\name{yuima.multimodel-class}
\docType{class}
\alias{yuima.multimodel-class}
\alias{yuima.multimodel}
\alias{initialize,yuima.multimodel-method}
\alias{simulate,yuima.multimodel-method}

\title{Class for the mathematical description of Multi dimensional Jump Diffusion processes}
\description{
  The \code{yuima.multimodel} class is a class of the  \pkg{yuima} package that extends the \code{\link{yuima.model-class}}.
}
\section{Slots}{
  \describe{
    \item{\code{drift}:}{always \code{expression((0))}.}
    \item{\code{diffusion}:}{a list of \code{expression((0))}.}
    \item{\code{hurst}:}{always \code{h=0.5}, but ignored for this model.}
    \item{\code{jump.coeff}:}{set according to \code{scale} in \code{\link{setPoisson}}.}
    \item{\code{measure}:}{a list containting the intensity measure and the jump distribution.}
    \item{\code{measure.type}:}{always \code{"CP"}.}
    \item{state.variable}{a vector of names identifying the names used to
     denote the state variable in the drift and diffusion specifications.}
    \item{\code{parameter}:}{which is a short name for ``parameters'', is an
	 object of class \code{\link{model.parameter-class}}. For more details see
	 \code{\link{model.parameter-class}} documentation page.}
    \item{\code{state.variable}:}{identifies the state variables in the \R{}
     expression.}
    \item{\code{jump.variable}:}{identifies the variable for the jump
	 coefficient.}
    \item{\code{time.variable}:}{the time variable.}
    \item{\code{noise.number}:}{denotes the number of sources of noise.}
    \item{\code{equation.number}:}{denotes the dimension of the stochastic
     differential equation.}
    \item{\code{dimension}:}{the dimensions of the parameter given in the
     \code{parameter} slot.}
    \item{\code{solve.variable}:}{identifies the variable with respect to which
	 the stochastic differential equation has to be solved.}
    \item{\code{xinit}:}{contains the initial value of the stochastic
	 differential equation.}
    \item{\code{J.flag}:}{wheather jump.coeff include jump.variable.}
  }
}
\section{Methods}{
  \describe{
    \item{simulate}{simulation method. For more information see
	  \code{\link{simulate}}.}
%    \item{toLatex}{This method converts an object of \code{yuima.poisson-class} to character vectors with LaTeX markup.}
    \item{qmle}{Quasi maximum likelihood estimation procedure. For more information see \code{\link{qmle}}. }
  }
}
\author{The YUIMA Project Team}
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

mod <- setModel(drift = matrix(c("0","0","0"),3,1), diffusion = NULL,
  jump.coeff = matrix(c("1","0","0","0","1","-1","1","0","0"),3,3),
  measure = list( intensity = "lambda1", df = "dmyexp(z,sig1,sig2,sig3)"),
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

mod1 <- setModel(drift=matrix(c("0.1*(0.01-X1)",
  "0.05*(1-X2)","0.1*(0.1-X3)"),3,1), diffusion=NULL,
   jump.coeff = matrix(c("0.01","0","0","0","0.01",
                          "0","0","0","0.01"),3,3),
   measure = list(df="dmyexp1(z,sig1,sig2,sig3)"),
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

mod2 <- setModel(drift=matrix(c("0.1*(0.01-X1)",
  "0.05*(1-X2)","0.1*(0.1-X3)", "0", "0"),5,1), diffusion=NULL,
  jump.coeff = matrix(c("0.01","0","0","0","0",
                        "0","0.01","0","0","0",
                        "0","0","0.01","0","0",
                        "0","0","0","0.01","0",
                        "0","0","0","0","0.01"),5,5),
  measure = list(df = "dmyexample2(z,sig1,sig2,sig3)",
            intensity = "lambda1"),
  jump.variable = c("z"),
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
\keyword{classes}



