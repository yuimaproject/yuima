\name{adaBayes}
\alias{adaBayes}
\alias{adaBayes,yuima-method}
\title{Adaptive Bayes Estimator for the Parameters in SDE Model}
\description{
  The \code{adabayes.mcmc} class is a class of the  \pkg{yuima} package that extends the \code{mle-class}.}
\usage{
adaBayes(yuima, start, prior, lower, upper, fixed = numeric(0), envir = globalenv(),
method = "mcmc", iteration = NULL,mcmc,rate =1, rcpp = TRUE, 
algorithm = "randomwalk",center=NULL,sd=NULL,rho=NULL,
path = FALSE)
}
\arguments{
  \item{yuima}{a 'yuima' object.}
  \item{start}{initial suggestion for parameter values }
  \item{prior}{a list of prior distributions for the parameters specified by 'code'. Currently, dunif(z, min, max), dnorm(z, mean, sd), dbeta(z, shape1, shape2), dgamma(z, shape, rate) are available. }
  \item{lower}{a named list for specifying lower bounds of parameters}
  \item{upper}{a named list for specifying upper bounds of parameters}
  \item{fixed}{a named list of parameters to be held fixed during optimization.}
  \item{envir}{an environment where the model coefficients are evaluated.}
  \item{method}{\code{"nomcmc"} requires package \code{cubature} }
  \item{iteration}{number of iteration of Markov chain Monte Carlo method}
  \item{mcmc}{number of iteration of Markov chain Monte Carlo method}
  \item{rate}{a thinning parameter. Only the first n^rate observation will be used for inference. }
  \item{rcpp}{Logical value. If \code{rcpp = TRUE} (default), Rcpp code will be performed. Otherwise, usual R code will be performed. }
  \item{algorithm}{If \code{algorithm = "randomwalk"} (default), the random-walk Metropolis algorithm will be performed. If \code{algorithm = "MpCN"}, the Mixed preconditioned Crank-Nicolson algorithm will be performed.}
  \item{center}{A list of parameters used to center MpCN algorithm.}
  \item{sd}{A list for specifying the standard deviation of proposal distributions.}
  \item{path}{Logical value when \code{method = "mcmc"}. If \code{path=TRUE}, then the sample path for each variable will be included in the MCMC object in the output.}
  \item{rho}{A parameter used for MpCN algorithm.}
}
\details{
Calculate the Bayes estimator for stochastic processes by using  the quasi-likelihood function. The calculation is performed by the Markov chain Monte Carlo method. Currently, the Random-walk Metropolis algorithm  and the Mixed preconditioned Crank-Nicolson algorithm is implemented.}
\section{Slots}{
  \describe{
    \item{\code{mcmc}:}{is a list of MCMC objects for all estimated parameters.}
    \item{\code{accept_rate}:}{is a list acceptance rates for diffusion and drift parts.}
    \item{\code{call}:}{is an object of class \code{language}.}
    \item{\code{fullcoef}:}{is an object of class \code{list} that contains estimated parameters.}
    \item{\code{vcov}:}{is an object of class \code{matrix}.}
    \item{\code{coefficients}:}{is an object of class \code{vector} that contains estimated parameters.}
  }
}
\author{Kengo Kamatani with YUIMA project Team}
\note{
\code{algorithm = nomcmc} is unstable.
}
\references{
Yoshida, N. (2011). Polynomial type large deviation inequalities and quasi-likelihood analysis for stochastic differential equations. Annals of the Institute of Statistical Mathematics, 63(3), 431-479.
Uchida, M., & Yoshida, N. (2014). Adaptive Bayes type estimators of ergodic diffusion processes from discrete observations. Statistical Inference for Stochastic Processes, 17(2), 181-219.
Kamatani, K. (2017). Ergodicity of Markov chain Monte Carlo with reversible proposal. Journal of Applied Probability, 54(2).
}
\examples{
\dontrun{
set.seed(123)
b <- c("-theta1*x1+theta2*sin(x2)+50","-theta3*x2+theta4*cos(x1)+25")
a <- matrix(c("4+theta5","1","1","2+theta6"),2,2)
true = list(theta1 = 0.5, theta2 = 5,theta3 = 0.3,
            theta4 = 5, theta5 = 1, theta6 = 1)
lower = list(theta1=0.1,theta2=0.1,theta3=0,
             theta4=0.1,theta5=0.1,theta6=0.1)
upper = list(theta1=1,theta2=10,theta3=0.9,
             theta4=10,theta5=10,theta6=10)
start = list(theta1=runif(1),
             theta2=rnorm(1),
             theta3=rbeta(1,1,1),
             theta4=rnorm(1),
             theta5=rgamma(1,1,1),
             theta6=rexp(1))
yuimamodel <- setModel(drift=b,diffusion=a,state.variable=c("x1", "x2"),solve.variable=c("x1","x2"))
yuimasamp <- setSampling(Terminal=50,n=50*10)
yuima <- setYuima(model = yuimamodel, sampling = yuimasamp)
yuima <- simulate(yuima, xinit = c(100,80),
                  true.parameter = true,sampling = yuimasamp)
prior <-
  list(
    theta1=list(measure.type="code",df="dunif(z,0,1)"),
    theta2=list(measure.type="code",df="dnorm(z,0,1)"),
    theta3=list(measure.type="code",df="dbeta(z,1,1)"),
    theta4=list(measure.type="code",df="dgamma(z,1,1)"),
    theta5=list(measure.type="code",df="dnorm(z,0,1)"),
    theta6=list(measure.type="code",df="dnorm(z,0,1)")
  )
set.seed(123)
mle <- qmle(yuima, start = start, lower = lower, upper = upper, method = "L-BFGS-B",rcpp=TRUE)
print(mle@coef)
center<-list(theta1=0.5,theta2=5,theta3=0.3,theta4=4,theta5=3,theta6=3)
sd<-list(theta1=0.001,theta2=0.001,theta3=0.001,theta4=0.01,theta5=0.5,theta6=0.5)
bayes <- adaBayes(yuima, start=start, prior=prior,lower=lower,upper=upper,
                  method="mcmc",mcmc=1000,rate = 1, rcpp = TRUE,
                   algorithm = "randomwalk",center = center,sd=sd,
                   path=TRUE)
print(bayes@fullcoef)
print(bayes@accept_rate)
print(bayes@mcmc$theta1[1:10])
}
}
\keyword{ts}

