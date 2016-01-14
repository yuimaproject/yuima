\name{rng}
\alias{rng}
\alias{rIG}
\alias{dIG}
\alias{rNIG}
\alias{dNIG}
\alias{rbgamma}
\alias{dbgamma}
\alias{rngamma}
\alias{dngamma}
\alias{rstable}

\title{Random numbers and densities}
\description{\code{simulate} function can use the specific random number generators to generate Levy paths.}
\usage{
rIG(x,delta,gamma)
dIG(x,delta,gamma)
rNIG(x,alpha,beta,delta,mu,Lambda)
dNIG(x,alpha,beta,delta,mu,Lambda)
rbgamma(x,delta.plus,gamma.plus,delta.minus,gamma.minus)
dbgamma(x,delta.plus,gamma.plus,delta.minus,gamma.minus)
rngamma(x,lambda,alpha,beta,mu,Lambda)
dngamma(x,lambda,alpha,beta,mu,Lambda)
rstable(x,alpha,beta,sigma,gamma)
%dstable(x,alpha,beta,sigma,gamma)
}
\arguments{
  \item{x}{Number of R.Ns to be geneated.}
  \item{delta}{parameter}
  \item{gamma}{parameter}
  \item{mu}{parameter}
  \item{Lambda}{parameter}
  \item{alpha}{parameter}
  \item{lambda}{parameter}
  \item{sigma}{parameter}
  \item{beta}{parameter}
  \item{delta.plus}{parameter}
  \item{gamma.plus}{parameter}
  \item{delta.minus}{parameter}
  \item{gamma.minus}{parameter}

  % \item{IG (inverse Gaussian)}{delta and gamma are positive valued parameter.}
  % \item{NIG (normal inverse Gaussian)}{alpha and delta are nonnegative number, beta and mu are vector and Lambda is matrix.}
  % \item{bgamma (bilateral gamma)}{All of parameters are positive number.}
  % \item{ngamma (normal gamma)}{lamdba and alpha are positive number, beta and mu are vector and Lambda is matrix.}
  % \item{stable}{Stable index 0<alpha<=2; Skewness -1<=beta<=1; Scale sigma>0; Location gamma being a real number.}
}


\details{
\code{IG} (inverse Gaussian): Delta and gamma are positive (the case of gamma=0 corresponds to the positive half stable, provided by the "rstable").

\code{NIG} (normal inverse Gaussian): It is defined by the normal mean-variance mixuture of inverse Gaussian distribution. The parameters alpha, beta, delta and mu express the heaviness of tails, degree of asymmetry, scale and location, respectively. They satisfy the following conditions: 
Lambda is positive definite with det(Lambda)=1; delta>0; alpha>0 with alpha^2-beta^T Lambda beta >0.

\code{bgamma} (bilateral gamma): Bilateral gamma distribution is defined by the difference of independent gamma distributions Gamma(delta.plus,gamma.plus) and Gamma(delta.minus,gamma.minus).

\code{ngamma} (normal gamma): Normal gamma distribution is defined by the normal mean-variance mixture of gamma distribution. The parameters satisfy the following conditions: 
Lambda is positive definite with det(Lambda)=1; lambda>0; alpha>0 with alpha^2-beta^T Lambda beta >0. Especially in the case of beta=0 it is variance gamma distribution.

\code{stable} (stable): Parameters alpha, beta, sigma and gamma express stability, degree of skewness, scale and location, respectively. They satisfy the following condition: 0<alpha<=2; -1<=beta<=1; scale>0; gamma is a real number.
}

\value{
\item{rXXX}{Collection of of random numbers or vectors}
\item{dXXX}{Density dunction}
}

\author{The YUIMA Project Team}
\note{
  Some density-plot functions ar estill missing: as for the non-Gaussian stable densities, one can use, e.g., stabledist package.
}

\examples{

set.seed(123)

# Ex 1. (One-dimensional standard Cauchy distribution)
# The value of parameters is alpha=1,beta=0,sigma=1,gamma=0.
# Choose the value of x.
x<-10 # the number of r.n
rstable(x,1,0,1,0)

# Ex 2. (One-dimensional Levy distribution)
# Choose the value of sigma, gamma, x.
# alpha = 0.5, beta=1
x<-10 # the number of r.n
beta <- 1
sigma <- 0.1
gamma <- 0.1
rstable(x,0.5,beta,sigma,gamma)

# Ex 3. (Symmetric bilateral gamma)
# delta=delta.plus=delta.minus, gamma=gamma.plus=gamma.minus.
# Choose the value of delta and gamma and x.
x<-10 # the number of r.n
rbgamma(x,1,1,1,1)

# Ex 4. (One-dimensional normal inverse Gaussian distribution)
# Lambda=1.
# Choose the value of parameters and x.
x<-10 # the number of r.n
rNIG(x,1,1,1,1)

# Ex 5. (Multi-dimensional normal inverse Gaussian distribution)
# Choose the value of parameters and x.
beta<-c(.5,.5)
mu<-c(0,0)
Lambda<-matrix(c(1,0,0,1),2,2)
x<-10 # the number of r.n
rNIG(x,1,beta,1,mu,Lambda)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
