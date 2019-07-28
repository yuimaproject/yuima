\encoding{UTF-8}
\name{JBtest}
\alias{JBtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove jumps and calculate the Gaussian quasi-likelihood estimator based on the Jarque-Bera normality test 
}
\description{
Remove jumps and calculate the Gaussian quasi-likelihood estimator based on the Jarque-Bera normality test 
}
\usage{
JBtest(yuima,start,lower,upper,alpha,skewness=TRUE,kurtosis=TRUE,withdrift=FALSE)
}
\arguments{
  \item{yuima}{a yuima object (diffusion with compound Poisson jumps).}
  \item{lower}{a named list for specifying lower bounds of parameters.}
  \item{upper}{a named list for specifying upper bounds of parameters.}
  \item{alpha}{Insert Description Here.}
  \item{start}{initial values to be passed to the optimizer.}
  \item{skewness}{use third moment information ? by default, skewness=TRUE}
  \item{kurtosis}{use fourth moment information ? by default, kurtosis=TRUE}
  \item{withdrift}{use drift information for constructing self-normalized residuals or not? by default, withdrift = FALSE}
}
\details{
This function removes large increments which are regarded as jumps based on the iterative Jarque-Bera normality test, and after that, calculates the Gaussian quasi maximum likelihood estimator.
}
\value{
\item{Removed}{Removed jumps and jump times}
\item{OGQMLE}{Gaussian quasi maximum likelihood estimator before jump removal}
\item{JRGQMLE}{Gaussian quasi maximum likelihood estimator after jump removal}
\item{Figures}{For visualization, the jump points are presented. In addition, the histgram of the jump removed self-normalized residuals, transition of the estimators and the logarithm of Jarque-Bera statistics are given as figures}
}

\references{
Masuda, H. (2013). Asymptotics for functionals of self-normalized residuals of discretely observed stochastic processes. 
                   Stochastic Processes and their Applications 123 (2013), 2752--2778

Masuda, H and Uehara, Y. (2018) Estimating Diffusion With Compound Poisson Jumps Based On Self-normalized Residuals, arXiv:1802.03945
}
\author{
The YUIMA Project Team

Contacts: Yuma Uehara \email{y-uehara@ism.ac.jp}
}



%\seealso{
%}
\examples{
set.seed(123)
mod <- setModel(drift="10-3*x",
                diffusion="theta*(2+x^2)/(1+x^2)",
               jump.coeff="1",
               measure=list(intensity="1",df=list("dunif(z, 3, 5)")),
               measure.type="CP")

T <- 10 ## Terminal
n <- 5000 ## generation size
samp <- setSampling(Terminal=T, n=n) ## define sampling scheme
yuima <- setYuima(model = mod, sampling = samp)

yuima <- simulate(yuima, xinit=1,true.parameter=list(theta=sqrt(2)), sampling = samp)

JBtest(yuima,start=list(theta=0.5),upper=c(theta=100)
,lower=c(theta=0),alpha=0.01)
}