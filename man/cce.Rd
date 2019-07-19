\name{cce}
\alias{cce}
%\alias{cce,list-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Nonsynchronous Cumulative Covariance Estimator
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function estimates the covariance between two Ito processes when they are observed at discrete times possibly nonsynchronously. It can apply to 
irregularly sampled one-dimensional data as a special case.
}
\usage{
cce(x, method="HY", theta, kn, g=function(x)min(x,1-x), refreshing = TRUE,
    cwise = TRUE, delta = 0, adj = TRUE, K, c.two, J = 1, c.multi, kernel, H,
    c.RK, eta = 3/5, m = 2, ftregion = 0, vol.init = NA,
    covol.init = NA, nvar.init = NA, ncov.init = NA, mn, alpha = 0.4,
    frequency = 300, avg = TRUE, threshold, utime, psd = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
the method to be used. See `Details'.
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
a numeric vector or matrix. If it is a matrix, each of its components   indicates the tuning parameter which determines the pre-averaging window lengths \code{kn} to be used for estimating the corresponding component. If it is a numeric vector, it is converted to a matrix as \code{(C+t(C))/2}, where \code{C=matrix(theta,d,d)} and \code{d=dim(x)}. The default value is \code{0.15} for the method \code{"PHY"} or \code{"PTHY"} following Christensen et al. (2013), while it is \code{1} for the method \code{"MRC"} following Christensen et al. (2010). 
}
  \item{kn}{
%%     ~~Describe \code{kn} here~~
an integer-valued vector or matrix indicating the pre-averaging window length(s). For the methods \code{"PHY"} or \code{"PTHY"}, see `Details' for the default value. For the method \code{"MRC"}, the default value is \code{ceiling(theta*n^(1+delta))}, where \code{n} is the number of the refresh times associated with the data minus 1.  
}
  \item{g}{
%%     ~~Describe \code{g} here~~
a function indicating the weight function to be used. The default value is the Bartlett window: \code{function(x)min(x,1-x)}.
}
  \item{refreshing}{
%%     ~~Describe \code{refreshing} here~~
logical. If \code{TRUE}, the data is pre-synchronized by the next-tick interpolation in the refresh times.
}
  \item{cwise}{
%%     ~~Describe \code{cwise} here~~
logical. If \code{TRUE}, the estimator is calculated componentwise.
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
a non-negative number indicating the order of the pre-averaging window length(s) \code{kn}.
}
  \item{adj}{
%%     ~~Describe \code{adj} here~~
logical. If \code{TRUE}, a finite-sample adjustment is performed. For the method \code{"MRC"}, see Christensen et al. (2010) for details. For the method \code{"TSCV"}, see Zhang (2011) and Zhang et al. (2005) for details.
}
  \item{K}{
%%     ~~Describe \code{K} here~~
a positive integer indicating the large time-scale parameter. The default value is \code{ceiling(c.two*n^(2/3))}, where \code{n} is the number of the refresh times associated with the data minus 1.
}
  \item{c.two}{
%%     ~~Describe \code{c.two} here~~
a positive number indicating the tuning parameter which determines the scale of the large time-scale parameter \code{K}. The default value is the average of the  numeric vector each of whose components is the roughly estimated optimal value in the sense of the minimizer of the theoretical asymptotic variance of the estimator of the corresponding diagonal component. The theoretical asymptotic variance is considered in the standard case and given by Eq.(63) of Zhang et al. (2005).
}
  \item{J}{
%%     ~~Describe \code{J} here~~
a positive integer indicating the small time-scale parameter.
}
  \item{c.multi}{
%%     ~~Describe \code{c.multi} here~~
a numeric vector or matrix. If it is a matrix, each of its components   indicates the tuning parameter which determines (the scale of) the number of the time scales to be used for estimating the corresponding component. If it is a numeric vector, it is converted to a matrix as \code{(C+t(C))/2}, where \code{C=matrix(c.multi,d,d)} and \code{d=dim(x)}. The default value is the numeric vector each of whose components is the roughly estimated optimal value in the sense of minimizing the theoretical asymptotic variance of the estimator of the corresponding diagonal component. The theoretical asymptotic variance is considered in the standard case and given by Eq.(37) of Zhang (2006).
}
  \item{kernel}{
%%     ~~Describe \code{kernel} here~~
a function indicating the kernel function to be used. The default value is the Parzan kernel, which is recommended in Barndorff-Nielsen et al. (2009, 2011).
}
  \item{H}{
%%     ~~Describe \code{H} here~~
a positive number indicating the bandwidth parameter. The default value is \code{c.RK*n^eta}, where \code{n} is the number of the refresh times associated with the data minus 1.
}
  \item{c.RK}{
%%     ~~Describe \code{c.RK} here~~
a positive number indicating the tuning parameter which determines the scale of the bandwidth parameter \code{H}. The default value is the average of the numeric vector each of whose components is the roughly estimated optimal value in the sense of minimizing the theoretical asymptotic variance of the estimator of the corresponding diagonal component. The theoretical asymptotic variance is considered in the standard case and given in Barndorff-Nielsen et al. (2009, 2011).
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
a positive number indicating the tuning parameter which determines the order of the bandwidth parameter \code{H}.
}
  \item{m}{
%%     ~~Describe \code{m} here~~
a positive integer indicating the number of the end points to be jittered.
}
  \item{ftregion}{
%%     ~~Describe \code{ftregion} here~~
a non-negative number indicating the length of the flat-top region. \code{ftregion=0} (the default) means that a non-flat-top realized kernel studied in Barndorff-Nielsen et al. (2011) is used. \code{ftregion=1/H} means that a flat-top realized kernel studied in Barndorff-Nielsen et al. (2008) is used. See Varneskov (2015) for other values. 
}
  \item{vol.init}{
%%     ~~Describe \code{vol.init} here~~
a numeric vector each of whose components indicates the initial value to be used to estimate the integrated volatility of the corresponding component, which is passed to the optimizer. %Components involving \code{NA} are filled with the rough estimated true values. 
}
  \item{covol.init}{
%%     ~~Describe \code{covol.init} here~~
a numeric matrix each of whose columns indicates the initial value to be used to estimate the integrated covariance of the corresponding component, which is passed to the optimizer. %Components involving \code{NA} are filled with the rough estimated true values.
}
  \item{nvar.init}{
%%     ~~Describe \code{nvar.init} here~~
a numeric vector each of whose components indicates the initial value to be used to estimate the variance of noise of the corresponding component, which is passed to the optimizer. %Components involving \code{NA} are filled with the rough estimated true values.
}
  \item{ncov.init}{
%%     ~~Describe \code{ncov.init} here~~
a numeric matrix each of whose columns indicates the initial value to be used to estimate the covariance of noise of the corresponding component, which is passed to the optimizer. %Components involving \code{NA} are filled with the rough estimated true values.
}
  \item{mn}{
%%     ~~Describe \code{mn} here~~
a positive integer indicating the number of terms to be used for calculating the SIML estimator. The default value is \code{ceiling(n^alpha)}, where \code{n} is the number of the refresh times associated with the data minus 1.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
a postive number indicating the order of \code{mn}.
}
  \item{frequency}{
%%     ~~Describe \code{frequency} here~~
a positive integer indicating the frequency (seconds) of the calendar time sampling to be used.
}
  \item{avg}{
%%     ~~Describe \code{avg} here~~
logical. If \code{TRUE}, the averaged subsampling estimator is calculated. Otherwise the simple sparsely subsampled estimator is calculated.
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
a numeric vector or list indicating the threshold parameter(s). Each of its components indicates the threshold parameter or process to be used for estimating the corresponding component. If it is a numeric vector, the elements in \code{threshold} are recycled if there are two few elements in \code{threshold}. The default value is determined following Koike (2014) (for the method \code{"THY"}) and Koike (2015) (for the method \code{"PTHY"}).
}
  \item{utime}{
%%     ~~Describe \code{utime} here~~
a positive number indicating what seconds the interval [0,1] corresponds to. 
%The default value is 23400 (this means 6.5 hours) if the time indices of the (first component) of \code{x} is numeric, otherwise 1.
The default value is the difference between the maximum and the minimum of the sampling times, multiplied by 23,400. Here, 23,400 seconds correspond to 6.5 hours, hence if the data is sampled on the interval [0,1], then the sampling interval is regarded as 6.5 hours. 
}
\item{psd}{
logical. If \code{TRUE}, the estimated covariance matrix \code{C} is converted to \code{(C\%*\%C)^(1/2)} for ensuring the positive semi-definiteness. In this case the absolute values of the estimated correlations are always ensured to be less than or equal to 1. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function is a method for objects of \code{\link{yuima.data-class}} and \code{\link{yuima-class}}.
It extracts the \code{data} slot when applied to a an object of \code{\link{yuima-class}}.

Typical usages are
\preformatted{
cce(x,psd=FALSE)
cce(x,method="PHY",theta,kn,g,refreshing=TRUE,cwise=TRUE,psd=FALSE)
cce(x,method="MRC",theta,kn,g,delta=0,avg=TRUE,psd=FALSE)
cce(x,method="TSCV",K,c.two,J=1,adj=TRUE,utime,psd=FALSE)
cce(x,method="GME",c.multi,utime,psd=FALSE)
cce(x,method="RK",kernel,H,c.RK,eta=3/5,m=2,ftregion=0,utime,psd=FALSE)
cce(x,method="QMLE",vol.init=NULL,covol.init=NULL,
    nvar.init=NULL,ncov.init=NULL,psd=FALSE)
cce(x,method="SIML",mn,alpha=0.4,psd=FALSE)
cce(x,method="THY",threshold,psd=FALSE)
cce(x,method="PTHY",theta,kn,g,threshold,refreshing=TRUE,cwise=TRUE,psd=FALSE)
cce(x,method="SRC",frequency=300,avg=TRUE,utime,psd=FALSE)
cce(x,method="SBPC",frequency=300,avg=TRUE,utime,psd=FALSE)
}

The default method is method \code{"HY"}, which is an implementation of the Hayashi-Yoshida estimator proposed in Hayashi and Yoshida (2005).
\cr
\cr
Method \code{"PHY"} is an implementation of the Pre-averaged Hayashi-Yoshida estimator proposed in Christensen et al. (2010).
\cr
\cr
Method \code{"MRC"} is an implementation of the Modulated Realized Covariance based on refresh time sampling proposed in Christensen et al. (2010).
\cr
\cr
Method \code{"TSCV"} is an implementation of the previous tick Two Scales realized CoVariance based on refresh time sampling proposed in Zhang (2011).
\cr
\cr
Method \code{"GME"} is an implementation of the Generalized Multiscale Estimator proposed in Bibinger (2011).
\cr
\cr
Method \code{"RK"} is an implementation of the multivariate Realized Kernel based on refresh time sampling proposed in Barndorff-Nielsen et al. (2011).
\cr
\cr
Method \code{"QMLE"} is an implementation of the nonparametric Quasi Maximum Likelihood Estimator proposed in Ait-Sahalia et al. (2010).
\cr
\cr
Method \code{"SIML"} is an implementation of the Separating Information Maximum Likelihood estimator proposed in Kunitomo and Sato (2013) with the basis of refresh time sampling.
\cr
\cr
Method \code{"THY"} is an implementation of the Truncated Hayashi-Yoshida estimator proposed in Mancini and Gobbi (2012).
\cr
\cr
Method \code{"PTHY"} is an implementation of the Pre-averaged Truncated Hayashi-Yoshida estimator, which is a thresholding version of the pre-averaged Hayashi-Yoshida estimator.
\cr
\cr
Method \code{"SRC"} is an implementation of the calendar time Subsampled Realized Covariance.
\cr
\cr
Method \code{"SBPC"} is an implementation of the calendar time Subsampled realized BiPower Covariation.
\cr
\cr
The rough estimation procedures for selecting the default values of the tuning parameters are based on those in Barndorff-Nielsen et al. (2009).
\cr
\cr
For the methods \code{"PHY"} or \code{"PTHY"}, the default value of \code{kn} changes depending on the values of \code{refreshing} and \code{cwise}. If both \code{refreshing} and \code{cwise} are \code{TRUE} (the default), the default value of \code{kn} is given by the matrix \code{ceiling(theta*N)}, where \code{N} is a matrix whose diagonal components are identical with the vector \code{length(x)-1} and whose \eqn{(i,j)}-th component is identical with the number of the refresh times associated with \eqn{i}-th and \eqn{j}-th components of \code{x} minus 1. If \code{refreshing} is \code{TRUE} while \code{cwise} is \code{FALSE}, the default value of \code{kn} is given by \code{ceiling(mean(theta)*sqrt(n))}, where \code{n} is the number of the refresh times associated with the data minus 1. If \code{refreshing} is \code{FALSE} while \code{cwise} is \code{TRUE}, the default value of \code{kn} is given by the matrix \code{ceiling(theta*N0)}, where \code{N0} is a matrix whose diagonal components are identical with the vector \code{length(x)-1} and whose \eqn{(i,j)}-th component is identical with \code{(length(x)[i]-1)+(length(x)[j]-1)}. If both \code{refreshing} and \code{cwise} are \code{FALSE}, the default value of \code{kn} is given by \code{ceiling(mean(theta)*sqrt(sum(length(x)-1)))} (following Christensen et al. (2013)).
\cr
\cr
For the method \code{"QMLE"}, the optimization of the quasi-likelihood function is implemented via \link{arima0} using the fact that it can be seen as an MA(1) model's one: See Hansen et al. (2008) for details.
}


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list with components:
\item{covmat}{the estimated covariance matrix}
\item{cormat}{the estimated correlation matrix}
}
\references{
%% ~put references to the literature/web site here ~
Ait-Sahalia, Y., Fan, J. and Xiu, D. (2010)
  High-frequency covariance estimates with noisy and asynchronous financial data,
  \emph{Journal of the American Statistical Association}, \bold{105}, no. 492, 1504--1517.
  
Barndorff-Nielsen, O. E., Hansen, P. R., Lunde, A. and Shephard, N. (2008)
  Designing realised kernels to measure the ex-post variation of equity prices in the presence of noise,
  \emph{Econometrica}, \bold{76}, no. 6, 1481--1536.
  
Barndorff-Nielsen, O. E., Hansen, P. R., Lunde, A. and Shephard, N. (2009)
  Realized kernels in practice: trades and quotes,
  \emph{Econometrics Journal}, \bold{12}, C1--C32.

Barndorff-Nielsen, O. E., Hansen, P. R., Lunde, A. and Shephard, N. (2011)
  Multivariate realised kernels: Consistent positive semi-definite estimators of the covariation of equity prices with noise and non-synchronous trading,
  \emph{Journal of Econometrics}, \bold{162}, 149--169.

Bibinger, M. (2011)
  Efficient covariance estimation for asynchronous noisy high-frequency data,
  \emph{Scandinavian Journal of Statistics}, \bold{38}, 23--45.
  
Bibinger, M. (2012)
  An estimator for the quadratic covariation of asynchronously observed Ito processes with noise: asymptotic distribution theory,
  \emph{Stochastic processes and their applications}, \bold{122}, 2411--2453.

Christensen, K., Kinnebrock, S. and Podolskij, M. (2010)
  Pre-averaging estimators of the ex-post covariance matrix in noisy diffusion models with non-synchronous data,
  \emph{Journal of Econometrics}, \bold{159}, 116--133.
  
Christensen, K., Podolskij, M. and Vetter, M. (2013)
  On covariation estimation for multivariate continuous Ito semimartingales with noise in non-synchronous observation schemes,
  \emph{Journal of Multivariate Analysis} \bold{120} 59--84.

Hansen, P. R., Large, J. and Lunde, A. (2008)
  Moving average-based estimators of integrated variance,
  \emph{Econometric Reviews}, \bold{27}, 79--111.
  
Hayashi, T. and Yoshida, N. (2005)
  On covariance estimation of non-synchronously observed diffusion processes,
  \emph{Bernoulli}, \bold{11}, no. 2, 359--379.

Hayashi, T. and Yoshida, N. (2008)
  Asymptotic normality of a covariance estimator for nonsynchronously observed diffusion processes, 
  \emph{Annals of the Institute of Statistical Mathematics}, \bold{60}, no. 2, 367--406.

Koike, Y. (2016)
  Estimation of integrated covariances in the simultaneous presence of nonsynchronicity, microstructure noise and jumps,
  \emph{Econometric Theory}, \bold{32}, 533--611.
  
Koike, Y. (2014)
  An estimator for the cumulative co-volatility of asynchronously observed semimartingales with jumps,
  \emph{Scandinavian Journal of Statistics}, \bold{41}, 460--481.
  
Kunitomo, N. and Sato, S. (2013)
  Separating information maximum likelihood estimation of realized volatility and covariance with micro-market noise,
  \emph{North American Journal of Economics and Finance}, \bold{26}, 282--309.
  
Mancini, C. and Gobbi, F. (2012)
  Identifying the Brownian covariation from the co-jumps given discrete observations,
  \emph{Econometric Theory}, \bold{28}, 249--273.
  
Varneskov, R. T. (2016)
  Flat-top realized kernel estimation of quadratic covariation with non-synchronous and noisy asset prices,
  \emph{Journal of Business & Economic Statistics}, \bold{34}, no.1, 1--22.
  
Zhang, L. (2006)
  Efficient estimation of stochastic volatility using noisy observations: a multi-scale approach,
  \emph{Bernoulli}, \bold{12}, no.6, 1019--1043.

Zhang, L. (2011)
  Estimating covariation: Epps effect, microstructure noise,
  \emph{Journal of Econometrics}, \bold{160}, 33--47.
  
Zhang, L., Mykland, P. A. and Ait-Sahalia, Y. (2005)
  A tale of two time scales: Determining integrated volatility with noisy high-frequency data,
  \emph{Journal of the American Statistical Association}, \bold{100}, no. 472, 1394--1411.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA Project Team
}
\note{
%%  ~~further notes~~
The example shows the central limit theorem for the nonsynchronous 
covariance estimator.  
%Estimation of the asymptotic variance with p-values and
Estimation of the asymptotic variance can be implemented by \code{\link{hyavar}}. 
The second-order correction will be provided in a future version of the package.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{setModel}}, \code{\link{setData}}, \code{\link{hyavar}}, \code{\link{lmm}}, \code{\link{cce.factor}}
}
\examples{
## Set a model
diff.coef.1 <- function(t, x1 = 0, x2 = 0) sqrt(1+t)
diff.coef.2 <- function(t, x1 = 0, x2 = 0) sqrt(1+t^2)
cor.rho <- function(t, x1 = 0, x2 = 0) sqrt(1/2)
diff.coef.matrix <- matrix(c("diff.coef.1(t,x1,x2)", 
"diff.coef.2(t,x1,x2) * cor.rho(t,x1,x2)", 
"", "diff.coef.2(t,x1,x2) * sqrt(1-cor.rho(t,x1,x2)^2)"), 2, 2) 
cor.mod <- setModel(drift = c("", ""), 
diffusion = diff.coef.matrix,solve.variable = c("x1", "x2")) 

set.seed(111) 

## We use a function poisson.random.sampling to get observation by Poisson sampling.
yuima.samp <- setSampling(Terminal = 1, n = 1200) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima) 
psample<- poisson.random.sampling(yuima, rate = c(0.2,0.3), n = 1000) 

## cce takes the psample and returns an estimate of the quadratic covariation. 
cce(psample)$covmat[1, 2]
##cce(psample)[1, 2]

## True value of the quadratic covariation.
cc.theta <- function(T, sigma1, sigma2, rho) { 
  tmp <- function(t) return(sigma1(t) * sigma2(t) * rho(t)) 
	integrate(tmp, 0, T) 
}

theta <- cc.theta(T = 1, diff.coef.1, diff.coef.2, cor.rho)$value 
cat(sprintf("theta =\%.5f\n", theta))

names(psample@zoo.data)






# Example. A stochastic differential equation with nonlinear feedback. 

## Set a model
drift.coef.1 <- function(x1,x2) x2
drift.coef.2 <- function(x1,x2) -x1
drift.coef.vector <- c("drift.coef.1","drift.coef.2")
diff.coef.1 <- function(t,x1,x2) sqrt(abs(x1))*sqrt(1+t)
diff.coef.2 <- function(t,x1,x2) sqrt(abs(x2))
cor.rho <- function(t,x1,x2) 1/(1+x1^2)
diff.coef.matrix <- matrix(c("diff.coef.1(t,x1,x2)", 
"diff.coef.2(t,x1,x2) * cor.rho(t,x1,x2)","", 
"diff.coef.2(t,x1,x2) * sqrt(1-cor.rho(t,x1,x2)^2)"), 2, 2) 
cor.mod <- setModel(drift = drift.coef.vector,
 diffusion = diff.coef.matrix,solve.variable = c("x1", "x2"))

## Generate a path of the process
set.seed(111) 
yuima.samp <- setSampling(Terminal = 1, n = 10000) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima, xinit=c(2,3)) 
plot(yuima)


## The "true" value of the quadratic covariation.
cce(yuima)

## We use the function poisson.random.sampling to generate nonsynchronous 
## observations by Poisson sampling.
psample<- poisson.random.sampling(yuima, rate = c(0.2,0.3), n = 3000) 

## cce takes the psample to return an estimated value  of the quadratic covariation. 
## The off-diagonal elements are the value of the Hayashi-Yoshida estimator. 
cce(psample)




# Example. Epps effect for the realized covariance estimator

## Set a model
drift <- c(0,0)

sigma1 <- 1
sigma2 <- 1
rho <- 0.5

diffusion <- matrix(c(sigma1,sigma2*rho,0,sigma2*sqrt(1-rho^2)),2,2)

model <- setModel(drift=drift,diffusion=diffusion,
                  state.variable=c("x1","x2"),solve.variable=c("x1","x2"))
                  
## Generate a path of the latent process
set.seed(116)

## We regard the unit interval as 6.5 hours and generate the path on it 
## with the step size equal to 2 seconds

yuima.samp <- setSampling(Terminal = 1, n = 11700) 
yuima <- setYuima(model = model, sampling = yuima.samp) 
yuima <- simulate(yuima)

## We extract nonsynchronous observations from the path generated above 
## by Poisson random sampling with the average duration equal to 10 seconds

psample <- poisson.random.sampling(yuima, rate = c(1/5,1/5), n = 11700)

## Hayashi-Yoshida estimator consistetly estimates the true correlation
cce(psample)$cormat[1,2]

## If we synchronize the observation data on some regular grid 
## by previous-tick interpolations and compute the correlation 
## by therealized covariance based on such synchronized observations, 
## we underestimate the true correlation (known as the Epps effect). 
## This is illustrated by the following examples.

## Synchronization on the grid with 5 seconds steps
suppressWarnings(s1 <- cce(subsampling(psample, sampling = setSampling(n = 4680)))$cormat[1,2])
s1

## Synchronization on the grid with 10 seconds steps
suppressWarnings(s2 <- cce(subsampling(psample, sampling = setSampling(n = 2340)))$cormat[1,2])
s2

## Synchronization on the grid with 20 seconds steps
suppressWarnings(s3 <- cce(subsampling(psample, sampling = setSampling(n = 1170)))$cormat[1,2])
s3

## Synchronization on the grid with 30 seconds steps
suppressWarnings(s4 <- cce(subsampling(psample, sampling = setSampling(n = 780)))$cormat[1,2])
s4

## Synchronization on the grid with 1 minute steps
suppressWarnings(s5 <- cce(subsampling(psample, sampling = setSampling(n = 390)))$cormat[1,2])
s5

plot(zoo(c(s1,s2,s3,s4,s5),c(5,10,20,30,60)),type="b",xlab="seconds",ylab="correlation",
main = "Epps effect for the realized covariance")



# Example. Non-synchronous and noisy observations of a correlated bivariate Brownian motion

## Generate noisy observations from the model used in the previous example
Omega <- 0.005*matrix(c(1,rho,rho,1),2,2) # covariance matrix of noise
noisy.psample <- noisy.sampling(psample,var.adj=Omega)
plot(noisy.psample)

## Hayashi-Yoshida estimator: inconsistent
cce(noisy.psample)$covmat

## Pre-averaged Hayashi-Yoshida estimator: consistent
cce(noisy.psample,method="PHY")$covmat

## Generalized multiscale estimator: consistent
cce(noisy.psample,method="GME")$covmat

## Multivariate realized kernel: consistent
cce(noisy.psample,method="RK")$covmat

## Nonparametric QMLE: consistent
cce(noisy.psample,method="QMLE")$covmat

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
