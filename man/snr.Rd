\encoding{UTF-8}
\name{snr}
\alias{snr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating self-normalized residuals for SDEs.
}
\description{
Calculate self-normalized residuals based on the Gaussian quasi-likelihood estimator.
}
\usage{
snr(yuima, start, lower, upper, withdrift)
}
\arguments{
  \item{yuima}{a yuima object.}
  \item{lower}{a named list for specifying lower bounds of parameters.}
  \item{upper}{a named list for specifying upper bounds of parameters.}
  \item{start}{initial values to be passed to the optimizer.}
  \item{withdrift}{use drift information for constructing self-normalized residuals. by default, withdrift = FALSE}
}
\details{
This function calculates the Gaussian quasi maximum likelihood estimator and associated self-normalized residuals.
}
\value{
\item{estimator}{Gaussian quasi maximum likelihood estimator}
\item{snr}{self-normalized residuals based on the Gaussian quasi maximum likelihood estimator}
}

\references{
Masuda, H. (2013). Asymptotics for functionals of self-normalized residuals of discretely observed stochastic processes. 
                   Stochastic Processes and their Applications 123 (2013), 2752--2778
}
\author{
The YUIMA Project Team

Contacts: Yuma Uehara \email{y-uehara@ism.ac.jp}
}



%\seealso{
%}
\examples{
# Test code (1. diffusion case)
yuima.mod <- setModel(drift="-theta*x",diffusion="theta1/sqrt(1+x^2)")
n <- 10000
ysamp <- setSampling(Terminal=n^(1/3),n=n)
yuima <- setYuima(model=yuima.mod, sampling=ysamp)
set.seed(123)
yuima <- simulate(yuima, xinit=0, true.parameter = list(theta=2,theta1=3))
start=list(theta=3,theta1=0.5)
lower=list(theta=1,theta1=0.3)
upper=list(theta=5,theta1=3)
res <- snr(yuima,start,lower,upper)
str(res)

# Test code (2.jump diffusion case)
a<-3
b<-5
mod <- setModel(drift="10-theta*x", #drift="10-3*x/(1+x^2)",
                diffusion="theta1*(2+x^2)/(1+x^2)",
                jump.coeff="1",
                # measure=list(intensity="10",df=list("dgamma(z, a, b)")),
                measure=list(intensity="10",df=list("dunif(z, a, b)")),
                measure.type="CP")

T <- 100 ## Terminal
n <- 10000 ## generation size
samp <- setSampling(Terminal=T, n=n) ## define sampling scheme
yuima <- setYuima(model = mod, sampling = samp)

yuima <- simulate(yuima, xinit=1,
                  true.parameter=list(theta=2,theta1=sqrt(2),a=a,b=b), 
                  sampling = samp)
start=list(theta=3,theta1=0.5)
lower=list(theta=1,theta1=0.3)
upper=list(theta=5,theta1=3)
res <- snr(yuima,start,lower,upper)
str(res)
}