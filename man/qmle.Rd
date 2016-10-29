\encoding{UTF-8}
\name{qmle}
%\alias{ql}
%\alias{ml.ql}
\alias{qmle}
\alias{quasilogl}
%\alias{LSE}
%\alias{ml.ql2}
\alias{rql}
\alias{lse}
\alias{pseudologlikelihood}
\alias{pseudologlikelihood.COGARCH}
%\alias{ql,ANY-method}
%\alias{ml.ql,ANY-method}
%\alias{ml.ql2,ANY-method}
%\alias{rql,ANY-method}
\title{Calculate quasi-likelihood and ML estimator of least squares estimator}
\description{Calculate the quasi-likelihood and estimate of the parameters of the
  stochastic differential equation by the maximum likelihood method or least squares estimator
  of the drift parameter.}
\usage{
%ml.ql(yuima,theta2,theta1,h,theta2.lim=matrix(c(0,1),1,2),theta1.lim=matrix(c(0,1),1,2),print=FALSE,method,param,interval)
%ql(yuima,theta2,theta1,h,print=FALSE,param)
%rql(yuima,theta2,theta1,ptheta2,ptheta1,h,print=FALSE,param,prevparam)
qmle(yuima, start, method = "BFGS", fixed = list(), print =
FALSE, lower, upper, joint = FALSE, Est.Incr =
"NoIncr", aggregation = TRUE, threshold = NULL, rcpp =
FALSE, ...)
quasilogl(yuima, param, print = FALSE, rcpp = FALSE)
lse(yuima, start, lower, upper, method = "BFGS", ...)
}
\arguments{
  \item{yuima}{a yuima object.}
%  \item{theta2,theta1}{parameters of the sdeModel.}
%  \item{h}{ time span of observations.}
%  \item{theta2.lim, theta1.lim}{matrixes to specify the domains of the
%    parameters. Vector can be available only if theta is a scalar.}
%  \item{ptheta2,ptheta1}{}
  \item{print}{you can see a progress of the estimation when print is TRUE.}
  \item{method}{see Details.}
  \item{param}{\code{list} of parameters for the  quasi loglikelihood.}
%  \item{interval}{}
%  \item{prevparam}{}
  \item{lower}{a named list for specifying lower bounds of parameters}
  \item{upper}{a named list for specifying upper bounds of parameters}
  \item{start}{initial values to be passed to the optimizer.}
  \item{fixed}{for conditional (quasi)maximum likelihood estimation.}
  \item{joint}{perform joint estimation or two stage estimation? by default \code{joint=FALSE}.}
  \item{Est.Incr}{If the yuima model is an object of \code{\link{yuima.carma-class}} or \code{\link{yuima.cogarch-class}}  the \code{qmle} returns an object of \code{\link{yuima.carma.qmle-class}}, \code{\link{cogarch.est.incr-class}},\code{\link{cogarch.est-class}}  or object of class \code{mle-class}. By default \code{Est.Incr="NoIncr"}, alternative values are \code{IncrPar} and \code{Incr}. }
  \item{aggregation}{If \code{aggregation=TRUE}, before the estimation of the levy parameters we aggregate the increments.}
  \item{threshold}{If the model has Compund Poisson type jumps, the threshold is
  used to perform thresholding of the increments.}
  \item{...}{passed to \code{\link{optim}} method. See Examples.}
  \item{rcpp}{use C++ code?}
}
\details{
%   A function ql calculate the quasi-likelihood of a time series data X with any
%   parameters. A function ml.pl estimates parameters of the sdeModel by
%   maximizing the quasi-likelihood.
  \code{qmle} behaves more likely the standard \code{mle} function in \pkg{stats4} and
  argument \code{method} is one of the methods available in \code{\link{optim}}.

  \code{lse} calculates least squares estimators of the drift parameters. This is
  useful for initial guess of \code{qmle} estimation.
  \code{quasilogl} returns the value of the  quasi loglikelihood for a given
  \code{yuima} object and list of parameters \code{coef}.

}
\value{
  \item{QL}{a real value.}
  \item{opt}{a list with components the same as 'optim' function.}
  \item{carmaopt}{if the model is an object of \code{\link{yuima.carma-class}}, \code{qmle} returns an object \code{\link{yuima.carma.qmle-class}}}
  \item{cogarchopt}{if the model is an object of \code{\link{yuima.cogarch-class}}, \code{qmle} returns an object of class \code{\link{cogarch.est-class}}. The estimates are obtained by  maximizing the pseudo-loglikelihood function as shown in Iacus et al. (2015)}
}
\references{
## Non-ergodic diffucion

Genon-Catalot, V., & Jacod, J. (1993). On the estimation of the diffusion coefficient for multi-dimensional diffusion processes. In Annales de l'IHP \enc{Probabilités}{Probabilit?s} et statistiques, 29(1), 119-151.

Uchida, M., & Yoshida, N. (2013). Quasi likelihood analysis of volatility and nondegeneracy of statistical random field. Stochastic Processes and their Applications, 123(7), 2851-2876.

## Ergodic diffusion

Kessler, M. (1997). Estimation of an ergodic diffusion from discrete observations. Scandinavian Journal of Statistics, 24(2), 211-229.

## Jump diffusion

Shimizu, Y., & Yoshida, N. (2006). Estimation of parameters for diffusion processes with jumps from discrete observations. Statistical Inference for Stochastic Processes, 9(3), 227-277.

Ogihara, T., & Yoshida, N. (2011). Quasi-likelihood analysis for the stochastic differential equation with jumps. Statistical Inference for Stochastic Processes, 14(3), 189-229.

## COGARCH

Iacus S. M., Mercuri L. and Rroji E.(2015) Discrete time approximation of a COGARCH (p, q) model and its estimation. \href{http://arxiv.org/abs/1511.00253}{http://arxiv.org/abs/1511.00253}

## CARMA

Iacus S. M., Mercuri L. (2015) Implementation of Levy CARMA model in Yuima package. Comp. Stat. (30) 1111-1141. \href{http://link.springer.com/article/10.1007/s00180-015-0569-7}{http://link.springer.com/article/10.1007/s00180-015-0569-7}
}

\author{The YUIMA Project Team}
\note{
  %The function ml.ql uses the function optim internally.
  The function qmle uses the function optim internally.

  The function qmle uses the function \code{\link{CarmaNoise}} internally for estimation of underlying Levy if the model is an object of \code{\link{yuima.carma-class}}.
}

\examples{
#dXt^e = -theta2 * Xt^e * dt + theta1 * dWt
diff.matrix <- matrix(c("theta1"), 1, 1)
ymodel <- setModel(drift=c("(-1)*theta2*x"), diffusion=diff.matrix,
  time.variable="t", state.variable="x", solve.variable="x")
n <- 100

ysamp <- setSampling(Terminal=(n)^(1/3), n=n)
yuima <- setYuima(model=ymodel, sampling=ysamp)
set.seed(123)
yuima <- simulate(yuima, xinit=1, true.parameter=list(theta1=0.3,
theta2=0.1))
QL <- quasilogl(yuima, param=list(theta2=0.8, theta1=0.7))
##QL <- ql(yuima, 0.8, 0.7, h=1/((n)^(2/3)))
QL

## another way of parameter specification
##param <- list(theta2=0.8, theta1=0.7)
##QL <- ql(yuima, h=1/((n)^(2/3)), param=param)
##QL


## old code
##system.time(
##opt <- ml.ql(yuima, 0.8, 0.7, h=1/((n)^(2/3)), c(0, 1), c(0, 1))
##)
##cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
##print(coef(opt))


system.time(
opt2 <- qmle(yuima, start=list(theta1=0.8, theta2=0.7), lower=list(theta1=0,theta2=0),
 upper=list(theta1=1,theta2=1), method="L-BFGS-B")
)
cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
print(coef(opt2))

## initial guess for theta2 by least squares estimator
tmp <- lse(yuima, start=list(theta2=0.7), lower=list(theta2=0), upper=list(theta2=1))
tmp

system.time(
opt3 <- qmle(yuima, start=list(theta1=0.8, theta2=tmp), lower=list(theta1=0,theta2=0),
 upper=list(theta1=1,theta2=1), method="L-BFGS-B")
)
cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
print(coef(opt3))


## perform joint estimation? Non-optimal, just for didactic purposes
system.time(
opt4 <- qmle(yuima, start=list(theta1=0.8, theta2=0.7), lower=list(theta1=0,theta2=0),
 upper=list(theta1=1,theta2=1), method="L-BFGS-B", joint=TRUE)
)
cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
print(coef(opt4))

## old code
##system.time(
##opt <- ml.ql(yuima, 0.8, 0.7, h=1/((n)^(2/3)), c(0, 1), c(0, 1), method="Newton")
##)
##cat(sprintf("\nTrue param. theta2 = .3, theta1 = .1\n"))
##print(coef(opt))


\dontrun{
###multidimension case
##dXt^e = - drift.matrix * Xt^e * dt + diff.matrix * dWt
diff.matrix <- matrix(c("theta1.1","theta1.2", "1", "1"), 2, 2)

drift.c <- c("-theta2.1*x1", "-theta2.2*x2", "-theta2.2", "-theta2.1")
drift.matrix <- matrix(drift.c, 2, 2)

ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, time.variable="t",
                   state.variable=c("x1", "x2"), solve.variable=c("x1", "x2"))
n <- 100
ysamp <- setSampling(Terminal=(n)^(1/3), n=n)
yuima <- setYuima(model=ymodel, sampling=ysamp)
set.seed(123)

##xinit=c(x1, x2) #true.parameter=c(theta2.1, theta2.2, theta1.1, theta1.2)
yuima <- simulate(yuima, xinit=c(1, 1),
true.parameter=list(theta2.1=0.5, theta2.2=0.3, theta1.1=0.6, theta1.2=0.2))

## theta2 <- c(0.8, 0.2) #c(theta2.1, theta2.2)
##theta1 <- c(0.7, 0.1) #c(theta1.1, theta1.2)
## QL <- ql(yuima, theta2, theta1, h=1/((n)^(2/3)))
## QL

## ## another way of parameter specification
## #param <- list(theta2=theta2, theta1=theta1)
## #QL <- ql(yuima, h=1/((n)^(2/3)), param=param)
## #QL

## theta2.1.lim <- c(0, 1)
## theta2.2.lim <- c(0, 1)
## theta1.1.lim <- c(0, 1)
## theta1.2.lim <- c(0, 1)
## theta2.lim <- t( matrix( c(theta2.1.lim, theta2.2.lim), 2, 2) )
## theta1.lim <- t( matrix( c(theta1.1.lim, theta1.2.lim), 2, 2) )

## system.time(
## opt <- ml.ql(yuima, theta2, theta1, h=1/((n)^(2/3)), theta2.lim, theta1.lim)
## )
## opt@coef

system.time(
opt2 <- qmle(yuima, start=list(theta2.1=0.8, theta2.2=0.2, theta1.1=0.7, theta1.2=0.1),
 lower=list(theta1.1=.1,theta1.2=.1,theta2.1=.1,theta2.2=.1),
 upper=list(theta1.1=4,theta1.2=4,theta2.1=4,theta2.2=4), method="L-BFGS-B")
)
opt2@coef
summary(opt2)

## unconstrained optimization
system.time(
opt3 <- qmle(yuima, start=list(theta2.1=0.8, theta2.2=0.2, theta1.1=0.7, theta1.2=0.1))
)
opt3@coef
summary(opt3)

quasilogl(yuima, param=list(theta2.1=0.8, theta2.2=0.2, theta1.1=0.7, theta1.2=0.1))

##system.time(
##opt <- ml.ql(yuima, theta2, theta1, h=1/((n)^(2/3)), theta2.lim, theta1.lim, method="Newton")
##)
##opt@coef
##

# carma(p=2,q=0) driven by a brownian motion without location parameter

mod0<-setCarma(p=2,
               q=0,
               scale.par="sigma")

true.parm0 <-list(a1=1.39631,
                 a2=0.05029,
                 b0=1,
                 sigma=0.23)

samp0<-setSampling(Terminal=100,n=250)
set.seed(123)
sim0<-simulate(mod0,
               true.parameter=true.parm0,
               sampling=samp0)

system.time(
carmaopt0 <- qmle(sim0, start=list(a1=1.39631,a2=0.05029,
                              b0=1,
                               sigma=0.23))
)

summary(carmaopt0)

# carma(p=2,q=1) driven by a brownian motion without location parameter

mod1<-setCarma(p=2,
               q=1)

true.parm1 <-list(a1=1.39631,
                  a2=0.05029,
                  b0=1,
                  b1=2)

samp1<-setSampling(Terminal=100,n=250)
set.seed(123)
sim1<-simulate(mod1,
               true.parameter=true.parm1,
               sampling=samp1)

system.time(
  carmaopt1 <- qmle(sim1, start=list(a1=1.39631,a2=0.05029,
                                     b0=1,b1=2),joint=TRUE)
)

summary(carmaopt1)

# carma(p=2,q=1) driven by a compound poisson process where the jump size is normally distributed.

mod2<-setCarma(p=2,
               q=1,
               measure=list(intensity="1",df=list("dnorm(z, 0, 1)")),
               measure.type="CP")

true.parm2 <-list(a1=1.39631,
                  a2=0.05029,
                  b0=1,
                  b1=2)

samp2<-setSampling(Terminal=100,n=250)
set.seed(123)
sim2<-simulate(mod2,
               true.parameter=true.parm2,
               sampling=samp2)

system.time(
  carmaopt2 <- qmle(sim2, start=list(a1=1.39631,a2=0.05029,
                                     b0=1,b1=2),joint=TRUE)
)

summary(carmaopt2)

# carma(p=2,q=1) driven by a normal inverse gaussian process
mod3<-setCarma(p=2,q=1,
               measure=list(df=list("rNIG(z, alpha, beta, delta1, mu)")),
               measure.type="code")
#

# True param
true.param3<-list(a1=1.39631,
                 a2=0.05029,
                 b0=1,
                 b1=2,
                 alpha=1,
                 beta=0,
                 delta1=1,
                 mu=0)

samp3<-setSampling(Terminal=100,n=200)
set.seed(123)

sim3<-simulate(mod3,
               true.parameter=true.param3,
               sampling=samp3)


carmaopt3<-qmle(sim3,start=true.param3)

summary(carmaopt3)

# Simulation and Estimation of COGARCH(1,1) with CP driven noise

# Model parameters
eta<-0.053
b1 <- eta
beta <- 0.04
a0 <- beta/b1
phi<- 0.038
a1 <- phi


# Definition

cog11<-setCogarch(p = 1,q = 1,
  measure = list(intensity = "1",
                 df = list("dnorm(z, 0, 1)")),
  measure.type = "CP",
  XinExpr=TRUE)

# Parameter
paramCP11 <- list(a1 = a1, b1 =  b1,
                  a0 = a0, y01 = 50.31)
# Sampling scheme
samp11 <- setSampling(0, 3200, n=64000)

# Simulation
set.seed(125)

SimTime11 <- system.time(
  sim11 <- simulate(object = cog11,
    true.parameter = paramCP11,
    sampling = samp11,
    method="mixed")
)

plot(sim11)

# Estimation

timeComp11<-system.time(
  res11 <- qmle(yuima = sim11,
    start = paramCP11,
    grideq = TRUE,
    method = "Nelder-Mead")
)

timeComp11

unlist(paramCP11)

coef(res11)

# COGARCH(2,2) model driven by CP

cog22 <- setCogarch(p = 2,q = 2,
  measure = list(intensity = "1",
                 df = list("dnorm(z, 0, 1)")),
  measure.type = "CP",
  XinExpr=TRUE)

# Parameter

paramCP22 <- list(a1 = 0.04, a2 = 0.001,
  b1 =  0.705, b2 = 0.1, a0 = 0.1, y01 = (1 + 2 / 3),
  y02=0)

# Use diagnostic.cog for checking the stat and positivity

check22 <- Diagnostic.Cogarch(cog22, param = paramCP22)

# Sampling scheme

samp22 <- setSampling(0, 3600, n = 64000)

# Simulation

set.seed(125)
SimTime22 <- system.time(
  sim22 <- simulate(object = cog22,
    true.parameter = paramCP22,
    sampling = samp22,
    method = "Mixed")
)

plot(sim22)

timeComp22 <- system.time(
  res22 <- qmle(yuima = sim22,
    start = paramCP22,
    grideq=TRUE,
    method = "Nelder-Mead")
)

timeComp22

unlist(paramCP22)

coef(res22)

}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
