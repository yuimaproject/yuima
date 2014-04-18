\name{CarmaNoise}
\alias{Recovering.Noise}
\alias{Carma.Recovering}
\alias{CarmaRecovNoise}
\alias{CarmaNoise}
\alias{Levy.Carma}
\title{Estimation for the underlying Levy in a carma model}
\description{Retrieve the increment of the underlying Levy for the carma(p,q) process using the approach developed in Brockwell et al.(2011)}
\usage{
CarmaNoise(yuima, param, data=NULL, NoNeg.Noise=FALSE)
}
\arguments{
  \item{yuima}{a yuima object or an object of \code{\link{yuima.carma-class}}.}
  \item{param}{\code{list} of parameters for the carma.}
  \item{data}{an object of class \code{\link{yuima.data-class}} contains the observations available at uniformly spaced time. If \code{data=NULL}, the default, the 'CarmaNoise' uses the data in an object of \code{\link{yuima.data-class}}.}
  \item{NoNeg.Noise}{Estimate a non-negative Levy-Driven Carma process. By default \code{NoNeg.Noise=FALSE}.}
}
\value{
  %\item{QL}{a real value.}
  \item{incr.Levy}{a numeric object contains the estimated increments.}
  }
\author{The YUIMA Project Team}
\note{
  %The function ml.ql uses the function optim internally.
  The function \code{qmle} uses the function \code{CarmaNoise} for estimation of underlying Levy in the carma model.
}

 \references{
               Brockwell, P., Davis, A. R. and Yang. Y. (2011)
Estimation for Non-Negative Levy-Driven CARMA Process, \emph{Journal of Business And Economic Statistics}, \bold{29} - 2, 250-259.
}


\examples{
\dontrun{
#Ex.1: Carma(p=3, q=0) process driven by a brownian motion.

mod0<-setCarma(p=3,q=0)

# We fix the autoregressive and moving average parameters
# to ensure the existence of a second order stationary solution for the process.

true.parm0 <-list(a1=4,a2=4.75,a3=1.5,b0=1)

# We simulate a trajectory of the Carma model.

numb.sim<-1000
samp0<-setSampling(Terminal=100,n=numb.sim)
set.seed(100)
incr.W<-matrix(rnorm(n=numb.sim,mean=0,sd=sqrt(100/numb.sim)),1,numb.sim)

sim0<-simulate(mod0,
               true.parameter=true.parm0,
               sampling=samp0, increment.W=incr.W)

#Applying the CarmaNoise

system.time(
  inc.Levy0<-CarmaNoise(sim0,true.parm0)
)

# We compare the orginal with the estimated noise increments 

par(mfrow=c(1,2))
plot(t(incr.W)[1:998],type="l", ylab="",xlab="time")
title(main="True Brownian Motion",font.main="1")
plot(inc.Levy0,type="l", main="Filtered Brownian Motion",font.main="1",ylab="",xlab="time")

# Ex.2: carma(2,1) driven  by a compound poisson
# where jump size is normally distributed and
# the lambda is equal to 1.

mod1<-setCarma(p=2,               
               q=1,
               measure=list(intensity="Lamb",df=list("dnorm(z, 0, 1)")),
               measure.type="CP") 

true.parm1 <-list(a1=1.39631, a2=0.05029,
                  b0=1,b1=2,
                  Lamb=1)

# We generate a sample path.

samp1<-setSampling(Terminal=100,n=200)
set.seed(123)
sim1<-simulate(mod1,
               true.parameter=true.parm1,
               sampling=samp1)

# We estimate the parameter using qmle.
carmaopt1 <- qmle(sim1, start=true.parm1)
summary(carmaopt1)
# Internally qmle uses CarmaNoise. The result is in 
plot(carmaopt1)

# Ex.3: Carma(p=2,q=1) with scale and location parameters 
# driven by a Compound Poisson
# with jump size normally distributed.
mod2<-setCarma(p=2,                
               q=1,
               loc.par="mu",
               scale.par="sig",
               measure=list(intensity="Lamb",df=list("dnorm(z, 0, 1)")),
               measure.type="CP") 

true.parm2 <-list(a1=1.39631,
                  a2=0.05029,
                  b0=1,
                  b1=2,
                  Lamb=1,
                  mu=0.5,
                  sig=0.23)
# We simulate the sample path 
set.seed(123)
sim2<-simulate(mod2,
               true.parameter=true.parm2,
               sampling=samp1)

# We estimate the Carma and we plot the underlying noise.

carmaopt2 <- qmle(sim2, start=true.parm2)
summary(carmaopt2)

# Increments estimated by CarmaNoise
plot(carmaopt2)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
