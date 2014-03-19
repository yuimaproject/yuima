\name{limiting.gamma}
\alias{limiting.gamma}

\title{calculate the value of limiting covariance matrices : Gamma}
\description{To confirm assysmptotic normality of theta estimators.}
\usage{
limiting.gamma(obj,theta,verbose=FALSE)
}
\arguments{
  \item{obj}{an yuima or yuima.model object.}
  \item{theta}{true theta}
  \item{verbose}{an option for display a verbose process.}
}
\details{
  Calculate the value of limiting covariance matrices Gamma.
  The returned values gamma1 and gamma2 are used to confirm assysmptotic normality of theta estimators.
  this program is limitted to 1-dimention-sde model for now.
}
\value{
  \item{gamma1}{a theoretical figure for variance of theta1 estimator}
  \item{gamma2}{a theoretical figure for variance of theta2 estimator}
}
\author{The YUIMA Project Team}
\note{
  we need to fix this routine.
}
\examples{

set.seed(123)

## Yuima
diff.matrix <- matrix(c("theta1"), 1, 1)
myModel <- setModel(drift=c("(-1)*theta2*x"), diffusion=diff.matrix, 
time.variable="t", state.variable="x")
n <- 100
mySampling <- setSampling(Terminal=(n)^(1/3), n=n)
myYuima <- setYuima(model=myModel, sampling=mySampling)
myYuima <- simulate(myYuima, xinit=1, true.parameter=list(theta1=0.6, theta2=0.3))

## theorical figure of theta
theta1 <- 3.5
theta2 <- 1.3

theta <- list(theta1, theta2)
lim.gamma <- limiting.gamma(obj=myYuima, theta=theta, verbose=TRUE)

## return theta1 and theta2 with list
lim.gamma$list

## return theta1 and theta2 with vector
lim.gamma$vec

}
\keyword{ts}
