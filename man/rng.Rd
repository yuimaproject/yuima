\encoding{UTF-8}
\name{rng}
\alias{rng}
\alias{rIG}
\alias{dIG}
\alias{rNIG}
\alias{dNIG}
\alias{rbgamma}
\alias{dbgamma}
\alias{rvgamma}
\alias{dvgamma}
\alias{rstable}
\alias{rpts}
\alias{rnts}

\title{Random numbers and densities}
\description{\code{simulate} function can use the specific random number generators to generate Levy paths.}
\usage{
rIG(x,delta,gamma)
dIG(x,delta,gamma)
rNIG(x,alpha,beta,delta,mu,Lambda)
dNIG(x,alpha,beta,delta,mu,Lambda)
rbgamma(x,delta.plus,gamma.plus,delta.minus,gamma.minus)
dbgamma(x,delta.plus,gamma.plus,delta.minus,gamma.minus)
rvgamma(x,lambda,alpha,beta,mu,Lambda)
dvgamma(x,lambda,alpha,beta,mu,Lambda)
rstable(x,alpha,beta,sigma,gamma)
rpts(x,alpha,a,b)
rnts(x,alpha,a,b,beta,mu,Lambda)
}
\arguments{
  \item{x}{Number of R.Ns to be geneated.}
  \item{a}{parameter}
  \item{b}{parameter}
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
  % \item{vgamma (variance gamma)}{lamdba and alpha are positive number, beta and mu are vector and Lambda is matrix.}
  % \item{stable}{Stable index 0<alpha<=2; Skewness -1<=beta<=1; Scale sigma>0; Location gamma being a real number.}
}


\details{
\code{IG} (inverse Gaussian): Delta and gamma are positive (the case of gamma=0 corresponds to the positive half stable, provided by the "rstable").

\code{NIG} (normal inverse Gaussian): It is defined by the normal mean-variance mixuture of inverse Gaussian distribution. The parameters alpha, beta, delta and mu express the heaviness of tails, degree of asymmetry, scale and location, respectively. They satisfy the following conditions:
Lambda is symmetric and positive definite with det(Lambda)=1; delta>0; alpha>0 with alpha^2-beta^T Lambda beta >0.

\code{bgamma} (bilateral gamma): Bilateral gamma distribution is defined by the difference of independent gamma distributions Gamma(delta.plus,gamma.plus) and Gamma(delta.minus,gamma.minus).

\code{vgamma} (variance gamma): Variance gamma distribution is defined by the normal mean-variance mixture of gamma distribution. The parameters satisfy the following conditions:
Lambda is symmetric and positive definite with det(Lambda)=1; lambda>0; alpha>0 with alpha^2-beta^T Lambda beta >0. Especially in the case of beta=0 it is variance gamma distribution.

\code{stable} (stable): Parameters alpha, beta, sigma and gamma express stability, degree of skewness, scale and location, respectively. They satisfy the following condition: 0<alpha<=2; -1<=beta<=1; scale>0; gamma is a real number.

\code{pts}(positive tempered stable): Positive tempered stable distribution is defined by the tilting of positive stable distribution. The parameters alpha, a and b express stability, scale and degree of tilting, respectively. They satisfy the following condition: 0<alpha<1; a>0; b>0.

\code{nts}(normal tempered stable): Normal tempered stable distribution is defined by the normal mean-variance mixture of positive tempered stable distribution. The parameters alpha, a, b, beta, mu and Lambda express stability, scale, degree of tilting, degree of asymemtry, location and degree of mixture, respectively. They satisfy the following condition: Lambda is symmetric and positive definite with det(Lambda)=1; 0<alpha<1; a>0; b>0.
}

\value{
\item{rXXX}{Collection of of random numbers or vectors}
\item{dXXX}{Density dunction}
}

\author{The YUIMA Project Team}
\note{
  Some density-plot functions are still missing: as for the non-Gaussian stable densities, one can use, e.g., stabledist package.
  The rejection-acceptance method is used for generating pts and nts. It should be noted that its acceptance rate decreases at exponential order as a and b become larger: specifically, the rate is given by exp( a*gamma(-alpha)*b^(alpha) )
}

\references{
## rstable

Chambers, John M., Colin L. Mallows, and B. W. Stuck.  (1976) A method for simulating stable random variables, Journal of the american statistical association, 71(354), 340-344.

Weron, \enc{Rafał}{Rafal}. (1996) On the Chambers-Mallows-Stuck method for simulating skewed stable random variables, Statistics & probability letters, 28.2, 165-171.

Weron, \enc{Rafał}{Rafal}. (2010) Correction to:" On the Chambers-Mallows-Stuck Method for Simulating Skewed Stable Random Variables", No. 20761, University Library of Munich, Germany.

## rbgamma, dbgamma

\enc{Küchler}{Kuchler}, U., & Tappe, S. (2008). Bilateral Gamma distributions and processes in financial mathematics. Stochastic Processes and their Applications, 118(2), 261-283.

\enc{Küchler}{Kuchler}, U., & Tappe, S. (2008). On the shapes of bilateral Gamma densities. Statistics & Probability Letters, 78(15), 2478-2484.

## rIG, dIG

Chhikara, R. (1988). The Inverse Gaussian Distribution: Theory: Methodology, and Applications (Vol. 95). CRC Press.

Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating random variates using transformations with multiple roots. The American Statistician, 30(2), 88-90.

## rNIG, dNIG, rvgamma, dvgamma (the elements of GH distribution)

Barndorff-Nielsen, O. E. (1997). Processes of normal inverse Gaussian type. Finance and stochastics, 2(1), 41-68.

Eberlein, E. (2001). Application of generalized hyperbolic \enc{Lévy}{Levy} motions to finance. In \enc{Lévy}{Levy} processes (pp. 319-336). \enc{Birkhäuser}{Birkhauser} Boston.

Madan, D. B., Carr, P. P., & Chang, E. C. (1998). The variance gamma process and option pricing. European finance review, 2(1), 79-105.

## rpts

Kawai, R., & Masuda, H. (2011). On simulation of tempered stable random variates. Journal of Computational and Applied Mathematics, 235(8), 2873-2887.

## rnts

Barndorff-Nielsen, O. E., & Shephard, N. (2001). Normal modified stable processes. Aarhus: MaPhySto, Department of Mathematical Sciences, University of Aarhus.

}

\examples{

set.seed(123)

# Ex 1. (One-dimensional standard Cauchy distribution)
# The value of parameters is alpha=1,beta=0,sigma=1,gamma=0.
# Choose the values of x.
x<-10 # the number of r.n
rstable(x,1,0,1,0)

# Ex 2. (One-dimensional Levy distribution)
# Choose the values of sigma, gamma, x.
# alpha = 0.5, beta=1
x<-10 # the number of r.n
beta <- 1
sigma <- 0.1
gamma <- 0.1
rstable(x,0.5,beta,sigma,gamma)

# Ex 3. (Symmetric bilateral gamma)
# delta=delta.plus=delta.minus, gamma=gamma.plus=gamma.minus.
# Choose the values of delta and gamma and x.
x<-10 # the number of r.n
rbgamma(x,1,1,1,1)

# Ex 4. ((Possibly skewed) variance gamma)
# lambda, alpha, beta, mu
# Choose the values of lambda, alpha, beta, mu and x.
x<-10 # the number of r.n
rvgamma(x,2,1,-0.5,0)

# Ex 5. (One-dimensional normal inverse Gaussian distribution)
# Lambda=1.
# Choose the parameter values and x.
x<-10 # the number of r.n
rNIG(x,1,1,1,1)

# Ex 6. (Multi-dimensional normal inverse Gaussian distribution)
# Choose the parameter values and x.
beta<-c(.5,.5)
mu<-c(0,0)
Lambda<-matrix(c(1,0,0,1),2,2)
x<-10 # the number of r.n
rNIG(x,1,beta,1,mu,Lambda)

# Ex 7. (Positive tempered stable)
# Choose the parameter values and x.
alpha<-0.7
a<-0.2
b<-1
x<-10 # the number of r.n
rpts(x,alpha,a,b)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
