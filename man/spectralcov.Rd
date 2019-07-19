\name{spectralcov}
\alias{spectralcov}
\alias{lmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Spectral Method for Cumulative Covariance Estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function implements the local method of moments proposed in Bibinger et al. (2014) to estimate the cummulative covariance matrix of a non-synchronously observed multi-dimensional Ito process with noise.
}
\usage{
lmm(x, block = 20, freq = 50, freq.p = 10, K = 4, interval = c(0, 1), 
    Sigma.p = NULL, noise.var = "AMZ", samp.adj = "direct", psd = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.
}
  \item{block}{
%%     ~~Describe \code{block} here~~
a positive integer indicating the number of the blocks which the observation interval is split into.
}
  \item{freq}{
%%     ~~Describe \code{freq} here~~
a positive integer indicating the number of the frequencies used to compute the final estimator.
}
  \item{freq.p}{
%%     ~~Describe \code{freq.p} here~~
a positive integer indicating the number of the frequencies used to compute the pilot estimator for the spot covariance matrix (corresponding to the number \eqn{J_n} in Eq.(29) from Altmeyer and Bibinger (2015)).
}
  \item{K}{
%%     ~~Describe \code{K} here~~
a positive integer indicating the number of the blocks used to compute the pilot estimator for the spot covariance matrix (corresponding to the number \eqn{K_n} in Eq.(29) from Altmeyer and Bibinger (2015)).
}
  \item{interval}{
%%     ~~Describe \code{interval} here~~
a vector indicating the observation interval. The first component represents the initial value and the second component represents the terminal value.
}
  \item{Sigma.p}{
%%     ~~Describe \code{Sigma.p} here~~
a \code{block} by \code{dim(x)} matrix giving the pilot estimates of the spot covariance matrix plugged into the optimal weight matrices. If \code{NULL} (the default), it is computed by using formula (29) from Altmeyer and Bibinger (2015). 
}
  \item{noise.var}{
%%     ~~Describe \code{noise.var} here~~
character string giving the method to estimate the noise variances. There are several options: \code{"AMZ"} (the default) uses equation (3.7) from Gatheral and Oomen (2010), i.e. the quasi-maximum likelihood estimator proposed by Ait-Sahalia et al. (2005) (see also Xiu (2010)). \code{"BR"} uses equation (3.9) from Gatheral and Oomen (2010), i.e. the sample average of the squared returns divided by 2, the estimator proposed by Bandi and Russel (2006). \code{"O"} uses equation (3.8) from Gatheral and Oomen (2010), i.e. another method-of-moments estimator proposed by Oomen (2006). It is also possible to directly specify the noise variances by setting this argument to a numeric vector. In this case the \code{i}-th component of \code{noise.var} must indicates the variance of the noise for the \code{i}-th component of the observation process. 
}
  \item{samp.adj}{
%%     ~~Describe \code{samp.adj} here~~
character string giving the method to adjust the effect of the sampling times on the variances of the spectral statistics for the noise part. The default method \code{"direct"} uses the local sums of the squares of the one-skip differences of the sampling times divided by 2, which directly appears in the representation of the variances of the spectral statistics for the noise part. Another choice is \code{"QVT"}, which uses the local quadratic variations of time as in Altmeyer and Bibinger (2015) and Bibinger et al. (2014).
}
  \item{psd}{
%%     ~~Describe \code{psd} here~~
logical. If \code{TRUE} (the default), the estimated covariance matrix and variance-covariance matrix are converted to their spectral absolute values to ensure their positive semi-definiteness. This procedure does not matter in terms of the asymptotic theory.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The default implementation is the adaptive version of the local method of moments estimator, which is only based on observation data. It is possible to implement oracle versions of the estimator by setting user-specified \code{Sigma.p} and/or \code{noise.var}. An example is given below.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{"yuima.specv"}, which is a list with the following elements:
\item{covmat}{the estimated covariance matrix}
\item{vcov}{the estimated variance-covariance matrix of \code{as.vector(covmat)}}
\item{Sigma.p}{the pilot estimates of the spot covariance matrix}
}

\references{
%% ~put references to the literature/web site here ~
Ait-Sahalia, Y., Mykland, P. A. and Zhang, L. (2005)
  How often to sample a continuous-time process in the presence of market microstructure noise,
  \emph{The Review of Financial Studies}, \bold{18}, 351--416.

Altmeyer, R. and Bibinger, M. (2015)
  Functional stable limit theorems for quasi-efficient spectral covolatility estimators,
    to appear in \emph{Stochastic processes and their applications}, doi:10.1016/j.spa.2015.07.009.
    
Bandi, F. M. and Russell, J. R. (2006)
  Separating microstructure noise from volatility,
  \emph{Journal of Financial Economics}, \bold{79}, 655--692.
    
Bibinger, M., Hautsch, N., Malec, P. and Reiss, M. (2014)
  Estimating the quadratic covariation matrix from noisy observations: local method of moments and efficiency,
  \emph{Annals of Statistics}, \bold{42}, 80--114.
  
Gatheral J. and Oomen, R. C. A. (2010)
  Zero-intelligence realized variance estimation,
  \emph{Finance Stochastics}, \bold{14}, 249--283.
  
Oomen, R. C. A. (2006)
  Properties of realized variance under alternative sampling schemes,
  \emph{Journal of Business and Economic Statistics}, \bold{24}, 219--237.
  
Reiss, M. (2011)
  Asymptotic equivalence for inference on the volatility from noisy observations,
  \emph{Annals of Statistics}, \bold{39}, 772--802.
  
Xiu, D. (2010)
  Quasi-maximum likelihood estimation of volatility with high frequency data,
  \emph{Journal of Econometrics}, \bold{159}, 235--250.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA Project Team
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{cce}}, \code{\link{setData}}
}
\examples{
# Example. One-dimensional and regular sampling case
# Here the simulated model is taken from Reiss (2011)

## Set a model
sigma <- function(t) sqrt(0.02 + 0.2 * (t - 0.5)^4)
modI <- setModel(drift = 0, diffusion = "sigma(t)")

## Generate a path of the process
set.seed(117)
n <- 12000
yuima.samp <- setSampling(Terminal = 1, n = n) 
yuima <- setYuima(model = modI, sampling = yuima.samp) 
yuima <- simulate(yuima, xinit = 0)

delta <- 0.01 # standard deviation of microstructure noise
yuima <- noisy.sampling(yuima, var.adj = delta^2) # generate noisy observations
plot(yuima)

## Estimation of the integrated volatility
est <- lmm(yuima)
est

## True integrated volatility and theoretical standard error
disc <- seq(0, 1, by = 1/n)
cat("true integrated volatility\n")
print(mean(sigma(disc[-1])^2))
cat("theoretical standard error\n")
print(sqrt(8*delta*mean(sigma(disc[-1])^3))/n^(1/4))

# Plotting the pilot estimate of the spot variance path
block <- 20
G <- seq(0,1,by=1/block)[1:block]
Sigma.p <- sigma(G)^2 # true spot variance
plot(zoo(Sigma.p, G), col = "blue",, xlab = "time", 
     ylab = expression(sigma(t)^2))
lines(zoo(est$Sigma.p, G))

## "Oracle" implementation
lmm(yuima, block = block, Sigma.p = Sigma.p, noise.var = delta^2)

# Example. Multi-dimensional case
# We simulate noisy observations of a correlated bivariate Brownian motion
# First we examine the regular sampling case since in this situsation the theoretical standard 
# error can easily be computed via the formulae given in p.88 of Bibinger et al. (2014)

## Set a model
drift <- c(0,0)

rho <- 0.5 # correlation

diffusion <- matrix(c(1,rho,0,sqrt(1-rho^2)),2,2)

modII <- setModel(drift=drift,diffusion=diffusion,
                  state.variable=c("x1","x2"),solve.variable=c("x1","x2"))
                  
## Generate a path of the latent process
set.seed(123)

## We regard the unit interval as 6.5 hours and generate the path on it 
## with the step size equal to 1 seconds

n <- 8000
yuima.samp <- setSampling(Terminal = 1, n = n) 
yuima <- setYuima(model = modII, sampling = yuima.samp) 
yuima <- simulate(yuima)

## Generate noisy observations
eta <- 0.05
yuima <- noisy.sampling(yuima, var.adj = diag(eta^2, 2))
plot(yuima)

## Estimation of the integrated covariance matrix
est <- lmm(yuima)
est

## Theoretical standard error
a <- sqrt(4 * eta * (sqrt(1 + rho) + sqrt(1 - rho)))
b <- sqrt(2 * eta * ((1 + rho)^(3/2) + (1 - rho)^(3/2)))
cat("theoretical standard error\n")
print(matrix(c(a,b,b,a),2,2)/n^(1/4))

## "Oracle" implementation
block <- 20
Sigma.p <- matrix(c(1,rho,rho,1),block,4,byrow=TRUE) # true spot covariance matrix
lmm(yuima, block = block, Sigma.p = Sigma.p, noise.var = rep(eta^2,2))

# Next we extract nonsynchronous observations from 
# the path generated above by Poisson random sampling
psample <- poisson.random.sampling(yuima, rate = c(1/2,1/2), n = n)

## Estimation of the integrated covariance matrix
lmm(psample)

## "Oracle" implementation
lmm(psample, block = block, Sigma.p = Sigma.p, noise.var = rep(eta^2,2))

## Other choices of tuning parameters (estimated values are not varied so much)
lmm(psample, block = 25)
lmm(psample, freq = 100)
lmm(psample, freq.p = 15)
lmm(psample, K = 8)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
