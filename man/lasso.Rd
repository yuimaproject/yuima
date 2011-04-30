\name{lasso}
\alias{lasso}
\title{Adaptive LASSO estimation for stochastic differential equations}
\description{Adaptive LASSO estimation for stochastic differential equations.}
\usage{
lasso(yuima, lambda0, start, delta=1, ...)
}
\arguments{
  \item{yuima}{a yuima object.}
  \item{lambda0}{a named list with penalty for each parameter.}
  \item{start}{initial values to be passed to the optimizer.}
  \item{delta}{controls the amount of shrinking in the adaptive sequences.}
  \item{...}{passed to \code{\link{optim}} method. See Examples.}
}
\details{
  
  \code{lasso} behaves more likely the standard \code{\link{qmle}} function in  and
  argument \code{method} is one of the methods available in \code{\link{optim}}.

  From initial guess of QML estimates, performs adaptive LASSO estimation using 
  the Least Squares Approximation (LSA) as in Wang and Leng (2007, JASA).
  
  
}
\value{
  \item{ans}{a list with both QMLE and LASSO estimates.}
}
\author{The YUIMA Project Team}
\examples{
##multidimension case
diff.matrix <- matrix(c("theta1.1","theta1.2", "1", "1"), 2, 2)

drift.c <- c("-theta2.1*x1", "-theta2.2*x2", "-theta2.2", "-theta2.1")
drift.matrix <- matrix(drift.c, 2, 2)

ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, time.variable="t",
                   state.variable=c("x1", "x2"), solve.variable=c("x1", "x2"))
n <- 100
ysamp <- setSampling(Terminal=(n)^(1/3), n=n)
yuima <- setYuima(model=ymodel, sampling=ysamp)
set.seed(123)

truep <- list(theta1.1=0.6, theta1.2=0,theta2.1=0.5, theta2.2=0)
yuima <- simulate(yuima, xinit=c(1, 1), 
 true.parameter=truep)


est <- lasso(yuima, start=list(theta2.1=0.8, theta2.2=0.2, theta1.1=0.7, theta1.2=0.1),
 lower=list(theta1.1=1e-10,theta1.2=1e-10,theta2.1=.1,theta2.2=1e-10),
 upper=list(theta1.1=4,theta1.2=4,theta2.1=4,theta2.2=4), method="L-BFGS-B")

# TRUE
unlist(truep)

# QMLE
round(est$mle,3)

# LASSO
round(est$lasso,3) 
}
 
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
