\name{hyavar}
\alias{hyavar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Asymptotic Variance Estimator for the Hayashi-Yoshida estimator
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function estimates the asymptotic variances of covariance and correlation estimates by the Hayashi-Yoshida estimator.
}
\usage{
hyavar(yuima, bw, nonneg = TRUE, psd = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{
%%     ~~Describe \code{yuima} here~~
an object of yuima-class or yuima.data-class.
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
a positive number or a numeric matrix. If it is a matrix, each component indicate the bandwidth parameter for the kernel estimators used to estimate the asymptotic variance of the corresponding component (necessary only for off-diagonal components). If it is a number, it is converted to a matrix as \code{matrix(bw,d,d)}, where \code{d=dim(x)}. The default value is the matrix whose \eqn{(i,j)}-th component is given by \eqn{min(n_i,n_j)^{0.45}}, where \eqn{n_i} denotes the number of the observations for the \eqn{i}-th component of the data.
}
  \item{nonneg}{
%%     ~~Describe \code{psd} here~~
logical. If \code{TRUE}, the asymptotic variance estimates for correlations are always ensured to be non-negative. See `Details'.
}

\item{psd}{
%%     ~~Describe \code{psd} here~~
passed to \code{\link{cce}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The precise description of the method used to estimate the asymptotic variances is as follows. 
For diagonal components, they are estimated by the realized quarticity multiplied by \code{2/3}. Its theoretical validity is ensured by Hayashi et al. (2011), for example. 
For off-diagonal components, they are estimated by the naive kernel approach descrived in Section 8.2 of Hayashi and Yoshida (2011). Note that the asymptotic covariance between a diagonal component and another component, which is necessary to evaluate the asymptotic variances of correlation estimates, is not provided in Hayashi and Yoshida (2011), but it can be derived in a similar manner to that paper.
\cr
\cr
If \code{nonneg} is \code{TRUE}, negative values of the asymptotic variances of correlations are avoided in the following way. The computed asymptotic varaince-covariance matrix of the vector \eqn{(HYii,HYij,HYjj)} is converted to its spectral absolute value. Here, \eqn{HYij} denotes the Hayashi-Yohida estimator for the \eqn{(i,j)}-th component.
\cr
\cr
The function also returns the covariance and correlation matrices calculated by the Hayashi-Yoshida estimator (using \code{\link{cce}}).
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
\item{avar.cov}{the estimated asymptotic variances for covariances}
\item{avar.cor}{the estimated asymptotic variances for correlations}
}
\references{
%% ~put references to the literature/web site here ~
Barndorff-Nilesen, O. E. and Shephard, N. (2004)
  Econometric analysis of realized covariation: High frequency based covariance, regression, and correlation in financial economics, 
  \emph{Econometrica}, \bold{72}, no. 3, 885--925.

Bibinger, M. (2011)
  Asymptotics of Asynchronicity, 
  technical report, Available at \href{https://arxiv.org/abs/1106.4222}{arXiv:1106.4222}.

Hayashi, T., Jacod, J. and Yoshida, N. (2011)
  Irregular sampling and central limit theorems for power variations: The continuous case,
  \emph{Annales de l'Institut Henri Poincare - Probabilites et Statistiques}, \bold{47}, no. 4, 1197--1218.

Hayashi, T. and Yoshida, N. (2011)
  Nonsynchronous covariation process and limit theorems,
  \emph{Stochastic processes and their applications}, \bold{121}, 2416--2454.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA Project Team
}
\note{
%%  ~~further notes~~
Construction of kernel-type estimators for off-diagonal components is implemented after pseudo-aggregation described in Bibinger (2011).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{setData}}, \code{\link{cce}}
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

## Constructing a 95\% confidence interval for the quadratic covariation from psample
result <- hyavar(psample)
thetahat <- result$covmat[1,2] # estimate of the quadratic covariation
se <- sqrt(result$avar.cov[1,2]) # estimated standard error
c(lower = thetahat + qnorm(0.025) * se, upper = thetahat + qnorm(0.975) * se)

## True value of the quadratic covariation.
cc.theta <- function(T, sigma1, sigma2, rho) { 
  tmp <- function(t) return(sigma1(t) * sigma2(t) * rho(t)) 
  integrate(tmp, 0, T) 
}

# contained in the constructed confidence interval
cc.theta(T = 1, diff.coef.1, diff.coef.2, cor.rho)$value

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


## The "true" values of the covariance and correlation.
result.full <- cce(yuima)
(cov.true <- result.full$covmat[1,2]) # covariance
(cor.true <- result.full$cormat[1,2]) # correlation

## We use the function poisson.random.sampling to generate nonsynchronous 
## observations by Poisson sampling.
psample<- poisson.random.sampling(yuima, rate = c(0.2,0.3), n = 3000) 

## Constructing 95\% confidence intervals for the covariation from psample
result <- hyavar(psample)
cov.est <- result$covmat[1,2] # estimate of covariance
cor.est <- result$cormat[1,2] # estimate of correlation
se.cov <- sqrt(result$avar.cov[1,2]) # estimated standard error of covariance
se.cor <- sqrt(result$avar.cor[1,2]) # estimated standard error of correlation

## 95\% confidence interval for covariance
c(lower = cov.est + qnorm(0.025) * se.cov,
 upper = cov.est + qnorm(0.975) * se.cov) # contains cov.true

## 95\% confidence interval for correlation
c(lower = cor.est + qnorm(0.025) * se.cor,
 upper = cor.est + qnorm(0.975) * se.cor) # contains cor.true

## We can also use the Fisher z transformation to construct a
## 95\% confidence interval for correlation
## It often improves the finite sample behavior of the asymptotic
## theory (cf. Section 4.2.3 of Barndorff-Nielsen and Shephard (2004))
z <- atanh(cor.est) # the Fisher z transformation of the estimated correlation
se.z <- se.cor/(1 - cor.est^2) # standard error for z (calculated by the delta method)
## 95\% confidence interval for correlation via the Fisher z transformation
c(lower = tanh(z + qnorm(0.025) * se.z), upper = tanh(z + qnorm(0.975) * se.z)) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
