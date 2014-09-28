\name{setPoisson}
\alias{setPoisson}

\title{
  Basic constructor for Compound Poisson processes
}

\description{
  'setPoisson' construct a Compound Poisson model specification for a process of the form:

  \code{Mt = m0+sum_{i=0}^Nt c*Y_{tau_i},  M0=m0}

  where Nt is a homogeneous or time-inhomogeneous Poisson process, tau_i is the
  sequence of random times of Nt and Y is a sequence of i.i.d. random jumps.
}

\usage{
setPoisson(intensity = 1, df = NULL, scale = 1, dimension=1, ...)
}

\arguments{
  \item{intensity}{either and expression or a numerical value representing the intensity function of the Poisson process Nt.}
  \item{df}{is the density of jump random variables Y.}
  \item{scale}{this is the scaling factor \code{c}.}
  \item{dimension}{this is the dimension of the jump component.}
  \item{...}{passed to \code{\link{setModel}}}
 }

\details{

An object of \code{\link{yuima.model-class}} where the \code{model} slot is of class \code{yuima.poisson-class}.

}

\value{
  \item{model}{an object of  \code{\link{yuima.model-class}}.}
}

\author{The YUIMA Project Team}

\examples{
Terminal <- 10
samp <- setSampling(T=Terminal,n=1000)

# Ex 1. (Simple homogeneous Poisson process)
mod1 <- setPoisson(intensity="lambda", df=list("dconst(z,1)"))
set.seed(123)
y1 <- simulate(mod1, true.par=list(lambda=1),sampling=samp)
plot(y1)

# scaling the jumps
mod2 <- setPoisson(intensity="lambda", df=list("dconst(z,1)"),scale=5)
set.seed(123)
y2 <- simulate(mod2, true.par=list(lambda=1),sampling=samp)
plot(y2)

# scaling the jumps through the constant distribution
mod3 <- setPoisson(intensity="lambda", df=list("dconst(z,5)"))
set.seed(123)
y3 <- simulate(mod3, true.par=list(lambda=1),sampling=samp)
plot(y3)

# Ex 2. (Time inhomogeneous Poisson process)
mod4 <- setPoisson(intensity="beta*(1+sin(lambda*t))", df=list("dconst(z,1)"))
set.seed(123)
lambda <- 3
beta <- 5
y4 <- simulate(mod4, true.par=list(lambda=lambda,beta=beta),sampling=samp)
par(mfrow=c(2,1))
par(mar=c(3,3,1,1))
plot(y4)
f <- function(t) beta*(1+sin(lambda*t))
curve(f, 0, Terminal, col="red")

# Ex 2. (Time inhomogeneous Compound Poisson process with Gaussian Jumps)
mod5 <- setPoisson(intensity="beta*(1+sin(lambda*t))", df=list("dnorm(z,mu,sigma)"))
set.seed(123)
y5 <- simulate(mod5, true.par=list(lambda=lambda,beta=beta,mu=0, sigma=2),sampling=samp)
plot(y5)
f <- function(t)  beta*(1+sin(lambda*t))
curve(f, 0, Terminal, col="red")

}
\keyword{ts}


