\name{IC}
\alias{IC}

\title{
Information criteria for the stochastic differential equation
}

\description{
Calculate the information criteria BIC, Quasi-BIC (QBIC) and CIC for the stochastic differential equation.
}

\usage{
IC(yuima, data = NULL, start, lower, upper, joint = FALSE, rcpp = FALSE, ...)
}

\arguments{
  \item{yuima}{
  a yuima object.
  }
  \item{data}{
  the data to be used.
  }
  \item{start}{
  a named list of the initial values of the parameters for optimization.
  }
  \item{lower}{
  a named list for specifying lower bounds of the parameters.
  }
  \item{upper}{
  a named list for specifying upper bounds of the parameters.
  }
  \item{joint}{
  perform joint parameter estimation or two stage parameter estimation? (default \code{joint=FALSE})
  If \code{two.step=TRUE}, this argument is not used.
  }
  \item{rcpp}{
  use C++ code? (default \code{rcpp=FALSE})
  }
  \item{\dots}{
  
  }
}

\details{
Please see specifications in \href{https://sites.google.com/site/shoichieguchi90en/specification}{https://sites.google.com/site/shoichieguchi90en/specification}
}

\value{
  \item{par}{
  the estimators of the parameters.
  }
  \item{BIC}{
  a value of BIC.
  }
  \item{QBIC}{
  a value of QBIC.
  }
  \item{CIC}{
  a value of CIC.
  }
}

\references{

## AIC, BIC

Akaike, H. (1973). Information theory and an extension of the maximum likelihood principle. In Second International Symposium on Information Theory (Tsahkadsor, 1971), 267-281. \href{http://link.springer.com/chapter/10.1007/978-1-4612-1694-0_15}{http://link.springer.com/chapter/10.1007/978-1-4612-1694-0_15}

Schwarz, G. (1978). Estimating the dimension of a model. The Annals of Statistics, 6(2), 461-464. \href{http://projecteuclid.org/euclid.aos/1176344136}{http://projecteuclid.org/euclid.aos/1176344136}

## BIC, Quasi-BIC

Eguchi, S. and Masuda, H. (2016). Schwarz type model comparison for LAQ models. \href{https://arxiv.org/abs/1606.01627v2}{arXiv:1606.01627v2}.

## CIC

Uchida, M. (2010). Contrast-based information criterion for ergodic diffusion processes from discrete observations. Annals of the Institute of Statistical Mathematics, 62(1), 161-187. \href{http://link.springer.com/article/10.1007/s10463-009-0245-1}{http://link.springer.com/article/10.1007/s10463-009-0245-1}
}

\author{
The YUIMA Project Team

Contacts: Shoichi Eguchi \email{eguchi@sigmath.es.osaka-u.ac.jp}
}

\note{
The function \code{IC} uses the function \code{\link{qmle}} with \code{method="L-BFGS-B"} internally.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}
\examples{

### Ex.1 
set.seed(123)

N <- 1000   # number of data
h <- N^(-2/3)  # sampling stepsize
Ter <- N*h  # terminal sampling time

## Data generate (dXt = -Xt*dt + exp((-2*cos(Xt) + 1)/2)*dWt)
mod <- setModel(drift="theta21*x", diffusion="exp((theta11*cos(x)+theta12)/2)")
samp <- setSampling(Terminal=Ter, n = N)
yuima <- setYuima(model=mod, sampling=setSampling(Terminal=Ter, n=50*N))
simu.yuima <- simulate(yuima, xinit=1, true.parameter=list(theta11=-2, theta12=1, 
                       theta21=-1), subsampling=samp)
Xt <- NULL
for(i in 1:(N+1)){
  Xt <- c(Xt, simu.yuima@data@original.data[50*(i-1)+1])
}

## Parameter settings
para.init <- list(theta11=runif(1,max=-1.5,min=-2.5), theta12=runif(1,max=1.5,min=0.5), 
                  theta21=runif(1,max=-0.5,min=-1.5))
para.low <- list(theta11=-7, theta12=-4, theta21=-6)
para.upp <- list(theta11=3, theta12=6, theta21=4)

## Ex.1.1 (dXt = (theta21*x)*dt + exp((theta11*cos(x)+theta12)/2)*dWt)
mod1 <- setModel(drift="theta21*x", diffusion="exp((theta11*cos(x)+theta12)/2)")
samp1 <- setSampling(Terminal=Ter, n = N)
yuima1 <- setYuima(model=mod1, sampling=samp1)
ic1 <- IC(yuima1, data=Xt, start=para.init, upper=para.upp, lower=para.low, rcpp=TRUE)
ic1

## Ex.1.2 (dXt = (theta21*x)*dt + exp(theta11*cos(x)/2)*dWt)
mod2 <- setModel(drift="theta21*x", diffusion="exp(theta11*cos(x)/2)")
samp2 <- setSampling(Terminal=Ter, n = N)
yuima2 <- setYuima(model=mod2, sampling=samp2)
ic2 <- IC(yuima2, data=Xt, start=para.init, upper=para.upp, lower=para.low, rcpp=TRUE)
ic2

\dontrun{
### Ex.2 (multidimansion case) 
set.seed(123)

N <- 3000   # number of data
h <- N^(-2/3)  # sampling stepsize
Ter <- N*h  # terminal sampling time

## Data generate
diff.coef.matrix <- matrix(c("beta1+1", "beta3*x1", "-beta2*x1", "beta4+1"), 2, 2)
drif.coef.vec <- c("alpha1*x1", "alpha2*x2")
mod <- setModel(drift = drif.coef.vec, diffusion = diff.coef.matrix, 
                state.variable = c("x1", "x2"), solve.variable = c("x1", "x2"))
samp <- setSampling(Terminal = Ter, n = N)
yuima <- setYuima(model = mod, sampling = setSampling(Terminal = Ter, n = 50*N))
simu.yuima <- simulate(yuima, xinit = c(1,1), true.parameter = list(alpha1=-2, alpha2=-1, 
                       beta1=1, beta2=1, beta3=1, beta4=2), subsampling = samp)
Xt <- matrix(0,(N+1),2)
for(i in 1:(N+1)){
  Xt[i,] <- simu.yuima@data@original.data[50*(i-1)+1,]
}

## Parameter settings
para.init <- list(alpha1 = runif(1,min=-3,max=-1), alpha2 = runif(1,min=-2,max=0), 
                  alpha3 = runif(1,min=-1,max=1), beta1 = runif(1,min=0,max=2), 
                  beta2 = runif(1,min=0,max=2), beta3 = runif(1,min=0,max=2), 
                  beta4 = runif(1,min=1,max=3))
para.low <- list(alpha1 = -5, alpha2 = -5, alpha3 = -5, beta1 = 0, beta2 = 0, beta3 = 0, beta4 = 0)
para.upp <- list(alpha1 = 5, alpha2 = 5, alpha3 = 5, beta1 = 5, beta2 = 5, beta3 = 5, beta4 = 5)

## Ex.2.1 
diff.coef.matrix1 <- matrix(c("beta1+1", "beta3*x1", "-beta2*x1", "beta4+1"), 2, 2)
drif.coef.vec1 <- c("alpha1*x1", "alpha2*x2+alpha3")
mod1 <- setModel(drift = drif.coef.vec1, diffusion = diff.coef.matrix1, 
                 state.variable = c("x1", "x2"), solve.variable = c("x1", "x2"))
samp1 <- setSampling(Terminal=Ter, n = N)
yuima1 <- setYuima(model=mod1, sampling=samp1)
ic1 <- IC(yuima1, data=Xt, start=para.init, upper=para.upp, lower=para.low, rcpp=TRUE)
ic1

## Ex.2.2 
diff.coef.matrix2 <- matrix(c("beta1+1", "beta3*x1", "-beta2*x1", "beta4+1"), 2, 2)
drif.coef.vec2 <- c("alpha1*x1", "alpha2*x2")
mod2 <- setModel(drift = drif.coef.vec2, diffusion = diff.coef.matrix2, 
                 state.variable = c("x1", "x2"), solve.variable = c("x1", "x2"))
samp2 <- setSampling(Terminal=Ter, n = N)
yuima2 <- setYuima(model=mod2, sampling=samp2)
ic2 <- IC(yuima2, data=Xt, start=para.init, upper=para.upp, lower=para.low, rcpp=TRUE)
ic2
}

}

\keyword{Information criteria}






