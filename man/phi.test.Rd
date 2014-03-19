\name{phi.test}
\alias{phi.test}
\title{Phi-divergence test statistic for stochastic differential equations}
\description{Phi-divergence test statistic for stochastic differential equations.}
\usage{
phi.test(yuima, H0, H1, phi, print=FALSE,...)
}
\arguments{
  \item{yuima}{a yuima object.}
  \item{H0}{a named list of parameter under H0.}
  \item{H1}{a named list of parameter under H1.}
  \item{phi}{the phi function to be used in the test. See Details.}
  \item{print}{you can see a progress of the estimation when print is \code{TRUE}.}
  \item{...}{passed to \code{\link{qmle}} function.}
}
\details{
  
  \code{phi.test} executes a Phi-divergence test. If \code{H1} is not specified
  this hypothesis is filled with the QMLE estimates.
  
  If \code{phi} is missing, then \code{phi(x)=1-x+x*log(x)} and the
  Phi-divergence statistic corresponds to the likelihood ratio test statistic.
}
\value{
  \item{ans}{an obkect of class \code{phitest}.}
}
\author{The YUIMA Project Team}
\examples{
model<- setModel(drift="t1*(t2-x)",diffusion="t3")
T<-10
n<-1000
sampling <- setSampling(Terminal=T,n=n)
yuima<-setYuima(model=model, sampling=sampling)

h0 <- list(t1=0.3, t2=1, t3=0.25)
X <- simulate(yuima, xinit=1, true=h0)
h1 <- list(t1=0.3, t2=0.2, t3=0.1)

phi1 <- function(x) 1-x+x*log(x)

phi.test(X, H0=h0, H1=h1,phi=phi1)
phi.test(X, H0=h0, phi=phi1, start=h0, lower=list(t1=0.1, t2=0.1, t3=0.1), 
   upper=list(t1=2,t2=2,t3=2),method="L-BFGS-B")
phi.test(X, H0=h1, phi=phi1, start=h0, lower=list(t1=0.1, t2=0.1, t3=0.1), 
  upper=list(t1=2,t2=2,t3=2),method="L-BFGS-B")
}
 
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
