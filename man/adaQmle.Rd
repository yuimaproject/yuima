\name{adaQmle}
\alias{adaQmle}
\title{Adaptive Quasi-Maximum Likelihood Estimation for Diffusion Models}
\description{
Computes adaptive estimators for drift and diffusion parameters using initial
Gaussian (V) or least-squares (W) contrasts and higher-order coupled (U) or
plug-in (V) contrasts.
}
\usage{
adaQmle(yuima, start, p, lower, upper, fixed = list(),
  initial.contrast = c(drift = "gaussian", diffusion = "gaussian"),
  refinement = "plugin",
  expansion = c("progressive", "terminal"),
  estimator = "optim", prior = NULL,
  method = "L-BFGS-B", envir = globalenv(),
  mcmc = 1000L, iteration = NULL, rate = 1,
  algorithm = c("randomwalk", "mpcn"), center = NULL,
  sd = NULL, rho = 0.8, path = FALSE, ...)
}
\arguments{
  \item{yuima}{A \code{yuima} object containing an observed diffusion.}
  \item{start}{A named list of starting parameter values.}
  \item{p}{An integer at least two: the exponent in the sampling condition
  \eqn{n h^p \to 0}. The internal orders are derived as
  \eqn{k_0=\lfloor p/2\rfloor} and
  \eqn{l_0=\lfloor(p-1)/2\rfloor}.}
  \item{lower, upper}{Optional named lists or vectors of parameter bounds.}
  \item{fixed}{A named list of parameters held fixed.}
  \item{initial.contrast}{One value, or drift and diffusion values (in that
  order), selected
  from \code{"gaussian"} (alias \code{"V"}) and
  \code{"least_squares"} (alias \code{"W"}).}
  \item{refinement}{\code{"plugin"} (alias \code{"V"}) uses V contrasts whose correction terms
  are evaluated at preceding estimates. \code{"coupled"} uses U contrasts
  whose correction terms depend on the current candidate (alias \code{"U"}).}
  \item{expansion}{\code{"progressive"} raises the contrast order at each
  step. \code{"terminal"} (alias \code{"full"}) uses the final derived order
  in every refinement step.}
  \item{estimator}{A per-step vector containing \code{"m"}/\code{"optim"}
  or \code{"b"}/\code{"bayes"}. A scalar is repeated; a short vector is
  extended by repeating its final entry.}
  \item{prior}{Named prior density specifications in the same code format as
  \code{adaBayes}. A missing prior gives a flat prior.}
  \item{method}{Optimization method passed to \code{optim}.}
  \item{envir}{Parent environment used to evaluate model coefficients.}
  \item{mcmc}{Number of MCMC draws in a Bayes step.}
  \item{iteration}{Optional alias for \code{mcmc}, retained for
  \code{adaBayes} compatibility.}
  \item{rate}{Initial Bayes sample exponent \eqn{G}, as in
  \code{adaBayes}. The first \eqn{n_0=\lfloor n^G\rfloor} observations are
  used in an initial Bayes step.}
  \item{algorithm}{\code{"randomwalk"} or \code{"mpcn"}.}
  \item{center}{Optional named centers for MpCN proposals.}
  \item{sd}{Optional named proposal standard deviations for random walk. For
  MpCN these are the diagonal preconditioning values used by
  \code{adaBayes}.}
  \item{rho}{MpCN persistence parameter in \eqn{(0,1)}.}
  \item{path}{Whether MCMC paths are retained in the result details.}
  \item{...}{Additional named arguments passed to \code{optim} at each
  maximum-likelihood step, for example
  \code{control = list(maxit = 300)}.}
}
\details{
Write \eqn{\Delta_i X=X_{t_i}-X_{t_{i-1}}}, let \eqn{h} be the sampling
interval, and put
\eqn{a(x,\alpha)=\sigma(x,\alpha)\sigma(x,\alpha)^{\mathsf T}}. Up to
parameter-independent constants, the initial Gaussian contrasts are
\deqn{V_{\alpha,0}=\sum_i\{\log\det(a_i)+h^{-1}(\Delta_iX)^{\mathsf T}
a_i^{-1}(\Delta_iX)\}}
and
\deqn{V_{\beta,0}=\sum_i h^{-1}(\Delta_iX-hb_i)^{\mathsf T}a_i^{-1}
(\Delta_iX-hb_i),}
where the diffusion estimate is plugged into \eqn{a_i} in the drift
contrast. The corresponding least-squares contrasts are
\deqn{W_{\alpha,0}=\sum_i\|h^{-1}\Delta_iX(\Delta_iX)^{\mathsf T}-a_i\|_F^2}
and
\deqn{W_{\beta,0}=\sum_i h^{-1}\|\Delta_iX-hb_i\|^2.}

For an initial Bayes step, put \eqn{G=\code{rate}},
\eqn{n_0=\lfloor n^G\rfloor}, and \eqn{q=\max(p,2/G)}. Since the implemented
contrast \eqn{Q} is twice the negative quasi-log likelihood, the diffusion
and drift log targets are respectively
\deqn{-\frac{1}{2}n_0^{-(1-2/(qG))}Q_\alpha+\log\pi_1(\alpha)}
and
\deqn{-\frac{1}{2}(n_0h)^{-(1-2/(qG))}Q_\beta+\log\pi_2(\beta).}
Thus the temperature multiplies the quasi-likelihood contribution only; the
prior density is not tempered.

With \code{algorithm="mpcn"}, the proposal follows the mixed preconditioned
Crank--Nicolson kernel used by \code{adaBayes}. Conditional on the current
value, a Gamma mixing precision is drawn from the quadratic radius about
\code{center}, followed by a Crank--Nicolson Gaussian proposal. The
Metropolis ratio includes the corresponding radial Hastings correction.

For higher-order steps, let \eqn{\mu_i^{(r)}} and \eqn{C_i^{(r)}} denote the
order-\eqn{r} expansions of the conditional mean of \eqn{\Delta_iX} and its
conditional covariance, and write \eqn{\bar C_i^{(r)}=C_i^{(r)}/h}. The
coupled U contrast has the Gaussian form
\deqn{U^{(r)}(\theta)=\sum_i\{\log\det(\bar C_i^{(r)}(\theta))+
h^{-1}(\Delta_iX-\mu_i^{(r)}(\theta))^{\mathsf T}
\bar C_i^{(r)}(\theta)^{-1}(\Delta_iX-\mu_i^{(r)}(\theta))\},}
interpreted as the required truncated power series in \eqn{h}.
With \code{refinement = "coupled"}, all correction terms are recomputed at
the current candidate. With \code{refinement = "plugin"}, the leading
coefficient for the parameter block currently being optimized remains
variable, while higher-order mean or second-moment corrections are fixed at
the preceding adaptive estimates. The implemented contrast is twice the
negative quasi-log likelihood, up to parameter-independent constants, and is
therefore minimized.

Diffusion parameters are estimated first. If a parameter occurs in both the
diffusion and drift coefficients, it is estimated in the diffusion block and
then plugged into every drift block; it is not re-estimated jointly.

Likelihood and least-squares sums are evaluated by compiled C++ kernels.
The conditional mean and covariance coefficients required by the high-order
contrasts are generated symbolically with the \pkg{calculus} package.
}
\value{
An object of class \code{yuima.qmle}. Its \code{details} slot contains the
input \code{p}, Bayes values \code{q}, \code{rate}, and \code{bayes.n0}, the
derived diffusion and drift orders, the actual contrast order of every step,
\eqn{nh^p}, the normalized options, and stage results.
}
\references{
Uchida, M. and Yoshida, N. (2012). Adaptive estimation of an ergodic
diffusion process based on sampled data. \emph{Stochastic Processes and their
Applications}, 122, 2885--2924.

Kamatani, K. and Uchida, M. (2015). Hybrid multi-step estimators for
stochastic differential equations based on sampled data. \emph{Statistical
Inference for Stochastic Processes}, 18, 177--204.

Kaino, Y., Uchida, M. and Yoshida, Y. (2017). Hybrid estimation for an
ergodic diffusion process based on reduced data. \emph{Bulletin of Informatics
and Cybernetics}, 49, 89--118.
}
\examples{
\dontrun{
model <- setModel(drift = "-b*x", diffusion = "a")
sampling <- setSampling(Terminal = 10, n = 1000)
yuima <- setYuima(model = model, sampling = sampling)
yuima <- simulate(yuima, xinit = 0,
  true.parameter = list(a = 1, b = 1))
fit <- adaQmle(yuima, start = list(a = 0.8, b = 0.8), p = 3,
  lower = list(a = 0.01, b = 0.01),
  upper = list(a = 3, b = 3),
  initial.contrast = c(drift = "gaussian", diffusion = "gaussian"),
  refinement = "plugin", expansion = "progressive",
  estimator = c("b", "m"), mcmc = 1000)
}
}
\keyword{models}
\keyword{optimize}
