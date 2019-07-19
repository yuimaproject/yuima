\name{cce.factor}
\alias{cce.factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
High-Dimensional Cumulative Covariance Estimator by Factor Modeling and Regularization
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function estimates the covariance and precision matrices of a high-dimesnional Ito process by factor modeling and regularization when it is observed at discrete times possibly nonsynchronously with noise. 
}
\usage{
cce.factor(yuima, method = "HY", factor = NULL, PCA = FALSE, 
           nfactor = "interactive", regularize = "glasso", taper, 
           group = 1:(dim(yuima) - length(factor)), lambda = "bic", 
           weight = TRUE, nlambda = 10, ratio, N, thr.type = "soft", 
           thr = NULL, tau = NULL, par.alasso = 1, par.scad = 3.7, 
           thr.delta = 0.01, frequency = 300, utime, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{
%%     ~~Describe \code{yuima} here~~
an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.
}
  \item{method}{
%%     ~~Describe \code{method} here~~
the method to be used in \code{\link{cce}}.
}
  \item{factor}{
%%     ~~Describe \code{factor} here~~
an integer or character vector indicating which components of \code{yuima} are factors. 
If \code{NULL}, no factor structure is taken account of.
}
  \item{PCA}{
%%     ~~Describe \code{PCA} here~~
logical. If \code{TRUE}, a principal component analysis is performed to construct factors.
}
  \item{nfactor}{
%%     ~~Describe \code{nfactor} here~~
the number of factors constructed when \code{PCA} is \code{TRUE}. 
If \code{nfactor = "interactive"}, the scree plot of the principal component analysis is depicted and the user can set this argument interactively. 
}
  \item{regularize}{
%%     ~~Describe \code{regularize} here~~
the regularizaton method to be used. 
Possible choices are \code{"glasso"} (the default), \code{"tapering"}, \code{"thresholding"} and \code{"eigen.cleaning"}. 
See `Details'.
}
  \item{taper}{
%%     ~~Describe \code{taper} here~~
the tapering matrix used when \code{regularize = "tapering"}. 
If missing, the tapering matrix is constructed according to \code{group}. See `Details'.
}
  \item{group}{
%%     ~~Describe \code{group} here~~
an integer vector having the length equal to \code{dim(yuima)-length(factor)}. 
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
the penalty parameter used when \code{regularize = "glasso"}. If it is \code{"aic"} (resp. \code{"bic"}), it is selected by minimizing the formally defined AIC (resp. BIC). 
See `Details'. 
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
logical. If \code{TRUE}. A weighted version is used for \code{regularize = "glasso"} as in Koike (2019).
}
  \item{nlambda}{
%%     ~~Describe \code{nlambda} here~~
a positive integer indicating the number of candidate penalty parameters for which AIC or BIC is evaluated when \code{lambda} is \code{"aic"} or \code{"bic"}. %See `Details'. 
}
  \item{ratio}{
%%     ~~Describe \code{ratio} here~~
a positive number indicating the ratio of the largest and smallest values in candidate penalty parameters for which AIC or BIC is evaluated when \code{lambda} is \code{"aic"} or \code{"bic"}. See `Details'. 
The default value is \code{sqrt(log(d)/N)}, where \code{d} is the dimension of \code{yuima}. 
}
  \item{N}{
%%     ~~Describe \code{N} here~~
a positive integer indicating the "effective" sampling size, which is necessary to evealuate AIC and BIC when \code{lambda} is \code{"aic"} or \code{"bic"}. 
In a standard situation, it is equal to the sample size \eqn{-} 1, but it might be different when the data are observed nonsynchronously and/or with noise. 
If missing, it is automatically determined according to \code{method}. %See `Details'.
}
  \item{thr.type}{
%%     ~~Describe \code{thr.type} here~~
a character string indicating the type of the thresholding method used when \code{regularize = "thresholding"}. 
Possible choices are \code{"hard"}, \code{"soft"}, \code{"alasso"} and \code{"scad"}. See Section 2.3 of Dai et al. (2019) for the definition of each method. 
}
  \item{thr}{
%%     ~~Describe \code{thr} here~~
a numeric matrix indicating the threshold levels used when \code{regularize = "thresholding"}. Its entries indicate the threshold levels for the corresponding entries of the covariance matrix (values for \eqn{\lambda} in the notation of Dai et al. (2019)). 
A single number is converted to the matrix with common entries equal to that number. 
If \code{NULL}, it is determined according to \code{tau}. 
See `Details'.
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
a number between 0 and 1 used to determine the threshold levels used when \code{regularize = "thresholding"} and \code{thr=NULL} (a value for \eqn{\tau} in the notation of Dai et al. (2019)). 
If \code{NULL}, it is determined by a grid search procedure as suggested in Section 4.3 of Dai et al. (2019). 
See `Details'.
}
  \item{par.alasso}{
%%     ~~Describe \code{par.alasso} here~~
the tuning parameter for \code{thr.type = "alasso"} (a value for \eqn{\eta} in the notation of Dai et al. (2019)). 
}
  \item{par.scad}{
%%     ~~Describe \code{par.scad} here~~
the tuning parameter for \code{thr.type = "scad"} (a value for \eqn{a} in the notation of Dai et al. (2019)). 
}
  \item{thr.delta}{
%%     ~~Describe \code{thr.delta} here~~
a positive number indicating the step size used in the grid serach procedure to determine \code{tau}. 
%See `Details'.
}
  \item{frequency}{
%%     ~~Describe \code{frequency} here~~
passed to \code{\link{cce}}.
}
  \item{utime}{
%%     ~~Describe \code{utime} here~~
passed to \code{\link{cce}}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
passed to \code{\link{cce}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
One basic approach to estimate the covariance matrix of high-dimensional time series is to take account of the factor structure and perform regularization for the residual covariance matrix. 
This function implements such an estimation procedure for high-frequency data modeled as a discretely observed semimartingale. 
Specifically, let \eqn{Y} be a \eqn{d}-dimensional semimartingale which describes the dynamics of the observation data. We consider the following continuous-time factor model:
\deqn{
Y_t = \beta X_t + Z_t,     0\le t\le T,
}
where \eqn{X} is an \eqn{r}-dimensional semimartingale (the factor process), \eqn{Z} is a \eqn{d}-dimensional semimartingale (the residual process), and \eqn{\beta} is a constant \eqn{d\times r} matrix (the factor loading matrix). We assume that \eqn{X} and \eqn{Z} are orthogonal in the sense that \eqn{[X,Z]_T=0}. Then, the quadratic covariation matrix of \eqn{Y} is given by
\deqn{
[Y,Y]_T=\beta[X,X]_T\beta^\top+[Z,Z]_T.
}
Also, \eqn{\beta} can be written as \eqn{\beta=[Y,X]_T[X,X]_T^{-1}}. Thus, if we have observation data both for \eqn{Y} and \eqn{X}, we can construct estimators for \eqn{[Y,Y]_T}, \eqn{[X,X]_T} and \eqn{\beta} by \code{\link{cce}}. Moreover, plugging these estimators into the above equation, we can also construct an estimator for \eqn{[Z,Z]_T}.  Since this estimator is often poor due to the high-dimensionality, we regularize it by some method. Then, by plugging the regularized estimator for \eqn{[Z,Z]_T} into the above equation, we obtain the final  estimator for \eqn{[Y,Y]_T}. 

Even if we do not have observation data for \eqn{X}, we can (at least formally) construct a pseudo factor process by performing principal component analysis for the initial estimator of \eqn{[Y,Y]_T}. See Ait-Sahalia and Xiu (2017) and Dai et al. (2019) for details. 

Currently, the following four options are available for the regularization method applied to the residual covariance matrix estimate:
\enumerate{

\item \code{regularize = "glasso"} (the default). 

This performs the glaphical Lasso. When \code{weight=TRUE} (the default), a weighted version of the graphical Lasso is performed as in Koike (2019). Otherwise, the standard graphical Lasso is performed as in Brownlees et al. (2018). 

If \code{lambda="aic"} (resp.~\code{lambda="bic"}), the penalty parameter for the graphical Lasso is selected by minimizing the formally defined AIC (resp.~BIC). 
The minimization is carried out by grid search, where the grid is determined as in Section 5.1 of Koike (2019).

The optimization problem in the graphical Lasso is solved by the GLASSOFAST algorithm of Sustik and Calderhead (2012), which is available from the package \pkg{glassoFast}. 

\item \code{regularize = "tapering"}. 

This performs tapering, i.e. taking the entry-wise product of the residual covariance matrix estimate and a tapering matrix specified by \code{taper}. 
See Section 3.5.1 of Pourahmadi (2011) for an overview of this method.

If \code{taper} is missing, it is constructed according to \code{group} as follows: \code{taper} is a 0-1 matrix and the \eqn{(i,j)}-th entry is equal to 1 if and only if \code{group[i]==group[j]}. Thus, by default it makes the residual covariance matrix diagonal.

\item \code{regularize = "thresholding"}. 

This performs thresholding, i.e. entries of the residual covariance matrix are shrunk toward 0 according to a thresholding rule (specified by \code{thr.type}) and a threshold level (spencified by \code{thr}). 

If \code{thr=NULL}, the \eqn{(i,j)}-th entry of \code{thr} is given by \eqn{\tau\sqrt{\hat{[Z^i,Z^i]}_T\hat{[Z^j,Z^j]}_T}}, where \eqn{\hat{[Z^i,Z^i]}_T} (resp. \eqn{\hat{[Z^j,Z^j]}_T}) denotes the \eqn{i}-th (resp. \eqn{j}-th) diagonal entry of the non-regularized estimator for the residual covariance matrix \eqn{[Z,Z]_T}, and \eqn{\tau} is a tuning parameter specified by \code{tau}. 

When \code{tau=NULL}, the value of \eqn{\tau} is set to the smallest value in the grid with step size \code{thr.delta} such that the regularized estimate of the residual covariance matrix becomes positive definite.

\item \code{regularize = "eigen.cleaning"}. 

This performs the eigenvalue cleaning algorithm described in Hautsch et al. (2012). 

}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list with components:
\item{covmat.y}{the estimated covariance matrix}
\item{premat.y}{the estimated precision matrix}
\item{beta.hat}{the estimated factor loading matrix}
\item{covmat.x}{the estimated factor covariance matrix}
\item{covmat.z}{the estimated residual covariance matrix}
\item{premat.z}{the estimated residual precision matrix}
\item{sigma.z}{the estimated residual covariance matrix before regularization}
\item{pc}{the variances of the principal components (it is \code{NULL} if \code{PCA = FALSE})}
}
\references{
%% ~put references to the literature/web site here ~
Ait-Sahalia, Y. and Xiu, D. (2017).
  Using principal component analysis to estimate a high dimensional factor model with high-frequency data,
  \emph{Journal of Econometrics}, \bold{201}, 384--399.
  
Brownlees, C., Nualart, E. and Sun, Y. (2018). 
  Realized networks,
  \emph{Journal of Applied Econometrics}, \bold{33}, 986--1006.

Dai, C., Lu, K. and Xiu, D. (2019).
  Knowing factors or factor loadings, or neither? Evaluating estimators of large covariance matrices with noisy and asynchronous data,
  \emph{Journal of Econometrics}, \bold{208}, 43--79.
  
Hautsch, N., Kyj, L. M. and Oomen, R. C. (2012).
  A blocking and regularization approach to high-dimensional realized covariance estimation,
  \emph{Journal of Applied Econometrics}, \bold{27}, 625--645.
  
Koike, Y. (2019). 
De-biased graphical Lasso for high-frequency data. 
\href{https://arxiv.org/abs/1905.01494}{arXiv:1905.01494}.

Pourahmadi, M. (2011). 
  Covariance estimation: The GLM and regularization perspectives. 
  \emph{Statistical Science}, \bold{26}, 369--387.

Sustik, M. A. and Calderhead, B. (2012).
  GLASSOFAST: An efficient GLASSO implementation,
  UTCSTechnical Report TR-12-29, The University of Texas at Austin.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA project Team
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{cce}}, \code{\link{lmm}}, \code{\link{glassoFast}}
}
\examples{
set.seed(123)

## Simulating a factor process (Heston model)
drift <- c("mu*S", "-theta*(V-v)")
diffusion <- matrix(c("sqrt(max(V,0))*S", "gamma*sqrt(max(V,0))*rho",
                       0, "gamma*sqrt(max(V,0))*sqrt(1-rho^2)"),
                    2,2)
mod <- setModel(drift = drift, diffusion = diffusion,
                state.variable = c("S", "V")) 
n <- 2340
samp <- setSampling(n = n)
heston <- setYuima(model = mod, sampling = samp)
param <- list(mu = 0.03, theta = 3, v = 0.09, 
              gamma = 0.3, rho = -0.6) 
result <- simulate(heston, xinit = c(1, 0.1), 
                   true.parameter = param)

zdata <- get.zoo.data(result) # extract the zoo data
X <- log(zdata[[1]]) # log-price process
V <- zdata[[2]] # squared volatility process


## Simulating a residual process (correlated BM)
d <- 100 # dimension
Q <- 0.1 * toeplitz(0.7^(1:d-1)) # residual covariance matrix
dZ <- matrix(rnorm(n*d),n,d) \%*\% chol(Q)/sqrt(n)
Z <- zoo(apply(dZ, 2, "diffinv"), samp@grid[[1]])


## Constructing observation data
b <- runif(d, 0.25, 2.25) # factor loadings
Y <- X \%o\% b + Z
yuima <- setData(cbind(X, Y))

# We subsample yuima to construct observation data
yuima <- subsampling(yuima, setSampling(n = 78))


## Estimating the covariance matrix (factor is known)
cmat <- tcrossprod(b) * mean(V[-1]) + Q # true covariance matrix
pmat <- solve(cmat) # true precision matrix

# (1) Regularization method is glasso (the default)
est <- cce.factor(yuima, factor = 1) 
norm(est$covmat.y - cmat, type = "2") 
norm(est$premat.y - pmat, type = "2")

# (2) Regularization method is tapering
est <- cce.factor(yuima, factor = 1, regularize = "tapering") 
norm(est$covmat.y - cmat, type = "2") 
norm(est$premat.y - pmat, type = "2")

# (3) Regularization method is thresholding
est <- cce.factor(yuima, factor = 1, regularize = "thresholding") 
norm(est$covmat.y - cmat, type = "2") 
norm(est$premat.y - pmat, type = "2")

# (4) Regularization method is eigen.cleaning
est <- cce.factor(yuima, factor = 1, regularize = "eigen.cleaning") 
norm(est$covmat.y - cmat, type = "2") 
norm(est$premat.y - pmat, type = "2")


## Estimating the covariance matrix (factor is unknown)
yuima2 <- setData(Y)

# We subsample yuima to construct observation data
yuima2 <- subsampling(yuima2, setSampling(n = 78))

# (A) Ignoring the factor structure (regularize = "glasso")
est <- cce.factor(yuima2) 
norm(est$covmat.y - cmat, type = "2") 
norm(est$premat.y - pmat, type = "2")

# (B) Estimating the factor by PCA (regularize = "glasso")
est <- cce.factor(yuima2, PCA = TRUE, nfactor = 1) # use 1 factor 
norm(est$covmat.y - cmat, type = "2") 
norm(est$premat.y - pmat, type = "2")

# One can interactively select the number of factors
# after implementing PCA (the scree plot is depicted)
# Try: est <- cce.factor(yuima2, PCA = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}% use one of  RShowDoc("KEYWORDS")
