\name{CPoint}
\alias{CPoint}
\alias{qmleL}
\alias{qmleR}
\title{Adaptive LASSO estimation for stochastic differential equations}
\description{Adaptive LASSO estimation for stochastic differential equations.}
\usage{
CPoint(yuima, param1, param2, print=FALSE, plot=FALSE)
qmleL(yuima, start, method="BFGS", t, print=FALSE, lower, upper, ...)
qmleR(yuima, start, method="BFGS", t, print=FALSE, lower, upper, ...)
}
\arguments{
  \item{yuima}{a yuima object.}
  \item{param1}{parameter values before the change point t}
  \item{param2}{parameter values after the change point t}
  \item{print}{produces some output.}
  \item{t}{right end point for \code{qmleL}, left end point for \code{qmleR}. See Details.}
  \item{plot}{plot change point test statistics? Default is \code{FALSE}.}
  \item{...}{passed to \code{\link{optim}} method. See Examples.}
  \item{lower}{a named list for specifying lower bounds of parameters}
  \item{upper}{a named list for specifying upper bounds of parameters}
  \item{start}{initial values to be passed to the optimizer.}
  \item{method}{see Details.}
}
\details{
  
  \code{CPoint}  estimates the change point using quasi-maximum likelihood approach.
  
  Function \code{qmleL} estimates the parameters in the diffusion matrix using
  observations up to time \code{t}.
  
  Function \code{qmleR} estimates the parameters in the diffusion matrix using
  observations from time \code{t} to the end.
  
  Arguments in both \code{qmleL} and \code{qmleR} follow the same rules
  as in \code{\link{qmle}}.
  
  
}
\value{
  \item{ans}{a list with change point instant, and paramters before and after
  the change point.}
}
\author{The YUIMA Project Team}
\examples{
diff.matrix <- matrix(c("theta1.1*x1","0*x2","0*x1","theta1.2*x2"), 2, 2)

drift.c <- c("1-x1", "3-x2")
drift.matrix <- matrix(drift.c, 2, 1)

ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, time.variable="t",
state.variable=c("x1", "x2"), solve.variable=c("x1", "x2"))
n <- 1000

set.seed(123)

t1 <- list(theta1.1=.1, theta1.2=0.2)
t2 <- list(theta1.1=.6, theta1.2=.6)

tau <- 0.4
ysamp1 <- setSampling(n=tau*n, Initial=0, delta=0.01)
yuima1 <- setYuima(model=ymodel, sampling=ysamp1)
yuima1 <- simulate(yuima1, xinit=c(1, 1), true.parameter=t1)

x1 <- yuima1@data@zoo.data[[1]]
x1 <- as.numeric(x1[length(x1)])
x2 <- yuima1@data@zoo.data[[2]]
x2 <- as.numeric(x2[length(x2)])

ysamp2 <- setSampling(Initial=n*tau*0.01, n=n*(1-tau), delta=0.01)
yuima2 <- setYuima(model=ymodel, sampling=ysamp2)

yuima2 <- simulate(yuima2, xinit=c(x1, x2), true.parameter=t2)


yuima <- yuima1
yuima@data@zoo.data[[1]] <- c(yuima1@data@zoo.data[[1]], yuima2@data@zoo.data[[1]][-1])
yuima@data@zoo.data[[2]] <- c(yuima1@data@zoo.data[[2]], yuima2@data@zoo.data[[2]][-1])

plot(yuima)

# estimation of change point for given parameter values
t.est <- CPoint(yuima,param1=t1,param2=t2, plot=TRUE)


low <- list(theta1.1=0, theta1.2=0)

# first state estimate of parameters using small 
# portion of data in the tails
tmp1 <- qmleL(yuima,start=list(theta1.1=0.3,theta1.2=0.5),t=1.5,low=low,method="L-BFGS-B")
tmp1
tmp2 <- qmleR(yuima,start=list(theta1.1=0.3,theta1.2=0.5), t=8.5,low=low,method="L-BFGS-B")
tmp2


# first stage changepoint estimator
t.est2 <- CPoint(yuima,param1=tmp1$par,param2=tmp2$par)
t.est2$tau


# second stage estimation of parameters given first stage
# change point estimator
tmp11 <- qmleL(yuima,start=as.list(tmp1$par), t=t.est2$tau-0.1,method="L-BFGS-B")
tmp11

tmp21 <- qmleR(yuima,start=as.list(tmp2$par), t=t.est2$tau+0.1,method="L-BFGS-B")
tmp21

# second stage estimator of the change point
CPoint(yuima,param1=tmp11$par,param2=tmp21$par)
}
 
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
