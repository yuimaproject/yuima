\name{fitCIR}
\alias{fitCIR}

\title{Calculate Preliminary Estimator and One-step Improvements of a Cox-Ingersoll-Ross Diffusion}

\description{
This is a function to simulate the preliminary estimator and the corresponding one step estimators based on the Newton-Raphson and the scoring method of the Cox-Ingersoll-Ross process given via the SDE 

	\eqn{\mathrm{d} X_t = (\alpha-\beta X_t)\mathrm{d} t + \sqrt{\gamma X_t}\mathrm{d} W_t}
		
with parameters \eqn{\beta>0,} \eqn{2\alpha>5\gamma>0} and a Brownian motion \eqn{(W_t)_{t\geq 0}}. This function uses the Gaussian quasi-likelihood, hence requires that data is sampled at high-frequency.
}

\usage{
fitCIR(data)
}

\arguments{
  \item{data}{
  a numeric matrix 
containing the realization of \eqn{(t_0,X_{t_0}), \dots,(t_n,X_{t_n})} with \eqn{t_j} denoting the \eqn{j}-th sampling times. \code{data[1,]} contains the sampling times \eqn{t_0,\dots, t_n} and \code{data[2,]} the corresponding value of the process \eqn{X_{t_0},\dots,X_{t_n}.} In other words \code{data[,j]=}\eqn{(t_j,X_{t_j})}. The observations should be equidistant. 
  }
}

\value{
  A list with three entries each contain a vector in the following order: The result of the preliminary estimator, Newton-Raphson method and the method of scoring.
 % A list with four entries. The first three entries each contain a vector in the following order: The result of the preliminary estimator, Newton-Raphson method and the method of scoring. The last entry contains the model,	an object of \code{\link{yuima.model-class}}.
  
  If the sampling points are not equidistant the function will return \code{'Please use equidistant sampling points'.}

}

\details{
The estimators calculated by this function can be found in the reference below.
}

\references{
Y. Cheng, N. Hufnagel, H. Masuda. Estimation of ergodic square-root diffusion under high-frequency sampling. Econometrics and Statistics, Article Number: 346 (2022).
}

\author{
Nicole Hufnagel 

Contacts: \email{nicole.hufnagel@math.tu-dortmund.de}
}

\examples{
#You can make use of the function simCIR to generate the data 
data <- simCIR(alpha=3,beta=1,gamma=1, n=5000, h=0.05, equi.dist=TRUE)
results <- fitCIR(data)
}

\keyword{data}






