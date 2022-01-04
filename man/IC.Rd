\name{IC}
\alias{IC}

\title{
Information criteria for the stochastic differential equation
}

\description{
Information criteria BIC, Quasi-BIC (QBIC) and CIC for the stochastic differential equation.
}

\usage{
IC(drif = NULL, diff = NULL, data = NULL, Terminal = 1, 
   add.settings = list(), start, lower, upper, ergodic = TRUE,
   stepwise = FALSE, weight = FALSE, rcpp = FALSE, ...)
}

\arguments{
  \item{drif}{
  a character vector that each element presents the candidate drift coefficient.
  }
  \item{diff}{
  a character vector that each element presents the candidate diffusion coefficient.
  }
  \item{data}{
  the data to be used.
  }
  \item{Terminal}{
  terminal time of the grid.
  }
  \item{add.settings}{
  details of model settings(see \code{\link{setModel}}).
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
  \item{ergodic}{
  whether the candidate models are ergodic SDEs or not(default \code{ergodic=TRUE}).
  }
  \item{stepwise}{
  specifies joint procedure or stepwise procedure(default \code{stepwise=FALSE}).
  }
  \item{weight}{
  calculate model weight? (default \code{weight=FALSE})
  }
  \item{rcpp}{
  use C++ code? (default \code{rcpp=FALSE})
  }
  \item{\dots}{
  
  }
}

\details{
Calculate the information criteria BIC, QBIC, and CIC for stochastic processes.
The calculation and model selection are performed by joint procedure or stepwise procedure.
}

\value{
  \item{BIC}{
  values of BIC for all candidates.
  }
  \item{QBIC}{
  values of QBIC for all candidates.
  }
  \item{CIC}{
  values of CIC for all candidates.
  }
  \item{model}{
  information of all candidate models.
  }
  \item{par}{
  quasi-maximum likelihood estimator for each candidate.
  }
  \item{weight}{
  model weights for all candidates.
  }
  \item{selected}{
  selected model number and selected drift and diffusion coefficients
  }

}

\references{

## AIC, BIC

Akaike, H. (1973). Information theory and an extension of the maximum likelihood principle. In Second International Symposium on Information Theory (Tsahkadsor, 1971), 267-281. \doi{10.1007/978-1-4612-1694-0_15}

Schwarz, G. (1978). Estimating the dimension of a model. The Annals of Statistics, 6(2), 461-464. \doi{10.1214/aos/1176344136}

## BIC, Quasi-BIC

Eguchi, S. and Masuda, H. (2018). Schwarz type model comparison for LAQ models. Bernoulli, 24(3), 2278-2327. \doi{10.3150/17-BEJ928}.

## CIC

Uchida, M. (2010). Contrast-based information criterion for ergodic diffusion processes from discrete observations. Annals of the Institute of Statistical Mathematics, 62(1), 161-187. \doi{10.1007/s10463-009-0245-1}

## Model weight

Burnham, K. P. and Anderson, D. R. (2002). Model Selection and Multimodel Inference. Springer-Verlag, New York.
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
\dontrun{
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

## Candidate coefficients
diffusion <- c("exp((theta11*cos(x)+theta12*sin(x)+theta13)/2)", 
               "exp((theta11*cos(x)+theta12*sin(x))/2)", 
               "exp((theta11*cos(x)+theta13)/2)", "exp((theta12*sin(x)+theta13)/2)")
drift <- c("theta21*x + theta22", "theta21*x")

## Parameter settings
para.init <- list(theta11=runif(1,max=5,min=-5), theta12=runif(1,max=5,min=-5), 
                  theta13=runif(1,max=5,min=-5), theta21=runif(1,max=-0.5,min=-1.5),
                  theta22=runif(1,max=-0.5,min=-1.5))
para.low <- list(theta11=-10, theta12=-10, theta13=-10, theta21=-5, theta22=-5)
para.upp <- list(theta11=10, theta12=10, theta13=10, theta21=-0.001, theta22=-0.001)

## Ex.1.1 Joint
ic1 <- IC(drif=drift, diff=diffusion, data=Xt, Terminal=Ter, start=para.init, lower=para.low, 
          upper=para.upp, stepwise = FALSE, weight = FALSE, rcpp = TRUE)
ic1

## Ex.1.2 Stepwise
ic2 <- IC(drif=drift, diff=diffusion, data=Xt, Terminal=Ter, 
          start=para.init, lower=para.low, upper=para.upp,
          stepwise = TRUE, weight = FALSE, rcpp = TRUE)
ic2


### Ex.2 (multidimansion case) 
set.seed(123)

N <- 3000   # number of data
h <- N^(-2/3)  # sampling stepsize
Ter <- N*h  # terminal sampling time

## Data generate
diff.coef.matrix <- matrix(c("beta1*x1+beta3", "1", "-1", "beta1*x1+beta3"), 2, 2)
drif.coef.vec <- c("alpha1*x1", "alpha2*x2")
mod <- setModel(drift = drif.coef.vec, diffusion = diff.coef.matrix, 
                state.variable = c("x1", "x2"), solve.variable = c("x1", "x2"))
samp <- setSampling(Terminal = Ter, n = N)
yuima <- setYuima(model = mod, sampling = setSampling(Terminal = N^(1/3), n = 50*N))
simu.yuima <- simulate(yuima, xinit = c(1,1), true.parameter = list(alpha1=-2, alpha2=-1, 
                       beta1=-1, beta3=2), subsampling = samp)
Xt <- matrix(0,(N+1),2)
for(i in 1:(N+1)){
  Xt[i,] <- simu.yuima@data@original.data[50*(i-1)+1,]
}

## Candidate coefficients
diffusion <- list(matrix(c("beta1*x1+beta2*x2+beta3", "1", "-1", "beta1*x1+beta2*x2+beta3"), 2, 2),
                  matrix(c("beta1*x1+beta2*x2", "1", "-1", "beta1*x1+beta2*x2"), 2, 2),
                  matrix(c("beta1*x1+beta3", "1", "-1", "beta1*x1+beta3"), 2, 2),
                  matrix(c("beta2*x2+beta3", "1", "-1", "beta2*x2+beta3"), 2, 2),
                  matrix(c("beta1*x1", "1", "-1", "beta1*x1"), 2, 2),
                  matrix(c("beta2*x2", "1", "-1", "beta2*x2"), 2, 2),
                  matrix(c("beta3", "1", "-1", "beta3"), 2, 2))
drift <- list(c("alpha1*x1", "alpha2*x2"), c("alpha1*x2", "alpha2*x1"))
modsettings <- list(state.variable = c("x1", "x2"), solve.variable = c("x1", "x2"))

## Parameter settings
para.init <- list(alpha1 = runif(1,min=-3,max=-1), alpha2 = runif(1,min=-2,max=0),
                  beta1 = runif(1,min=-2,max=0), beta2 = runif(1,min=0,max=2), 
                  beta3 = runif(1,min=1,max=3))
para.low <- list(alpha1 = -5, alpha2 = -5, beta1 = -5, beta2 = -5, beta3 = 1)
para.upp <- list(alpha1 = 0.01, alpha2 = -0.01, beta1 = 5, beta2 = 5, beta3 = 10)

## Ex.2.1 Joint
ic3 <- IC(drif=drift, diff=diffusion, data=Xt, Terminal=Ter, add.settings=modsettings, 
          start=para.init, lower=para.low, upper=para.upp, 
          weight=FALSE, rcpp=FALSE)
ic3

## Ex.2.2 Stepwise
ic4 <- IC(drif=drift, diff=diffusion, data=Xt, Terminal=Ter, add.settings=modsettings, 
             start=para.init, lower=para.low, upper=para.upp,
             stepwise = TRUE, weight=FALSE, rcpp=FALSE)
ic4

}
}

\keyword{Information criteria}






