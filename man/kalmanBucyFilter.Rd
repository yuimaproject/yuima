\name{kalmanBucyFilter}
\alias{kalmanBucyFilter}
\alias{kalmanBucyFilter,yuima-method}
\alias{yuima.kalmanBucyFilter-class}
\alias{mean,yuima.kalmanBucyFilter-method}
\alias{vcov,yuima.kalmanBucyFilter-method}
\title{
Kalman-Bucy Filter
}
\description{
Estimates values of unobserved variables from observed variables in a Linear State Space Model.
}
\usage{
kalmanBucyFilter(
  yuima, params, mean_init, vcov_init = NULL, delta.vcov.solve = 0.001,
  are = FALSE, explicit = FALSE, time_homogeneous = FALSE,
  env = globalenv()
)
}
\arguments{
  \item{yuima}{A \code{yuima} object. The class of \code{yuima@model} must be \code{yuima.linear_state_space_model}.}
  \item{params}{A list of numeric values specifying the names and values of the model parameters.}
  \item{mean_init}{A numeric value specifying the initial value of the unobserved variables.}
  \item{vcov_init}{A matrix specifying the initial variance-covariance of the unobserved variables.
    This argument is required if \code{are = FALSE}. If there is only one observed variable, it can be a single numeric value of length 1.}
  \item{delta.vcov.solve}{A numeric value specifying the step size used to solve the mean squared error of the estimator.}
  \item{are}{A logical value. If \code{TRUE}, the algebraic Riccati equation is solved.}
  \item{explicit}{A logical value. If \code{TRUE}, the filtering equation is solved using the explicit formula.
    The formula is available only if the coefficients are time-independent.}
  \item{time_homogeneous}{A logical value. If \code{TRUE}, calculations are performed under the assumption that the model coefficients are independent of time,
    reducing computational effort. This flag is used only if \code{are = FALSE}.}
  \item{env}{An environment object specifying the environment in which the model coefficients are evaluated.}
}
\value{
A \code{yuima.kalmanBucyFilter} object.
\item{model}{A \code{\link{yuima.linear_state_space_model}} object.}
\item{mean}{A \code{\link{ts}} object containing the estimated values of unobserved variables.}
\item{vcov}{An array object containing the estimated mean squared error of the estimator of unobserved variables.}
\item{mean.init}{A numeric value representing the initial value of unobserved variables.}
\item{vcov.init}{A matrix representing the initial mean squared error of the estimator of unobserved variables.}
\item{delta}{A numeric value representing the time step of observations.}
\item{data}{A \code{yuima.data} object.}
}
\references{
Liptser, R. S., & Shiryaev, A. N. (2001). \emph{Statistics of Random Processes: General Theory}. Springer.
}
\author{YUIMA TEAM}
\examples{
vcov_init <- matrix(0.1)
mean_init <- 0
a <- 1.5
b <- 0.3
c <- 1
sigma <- 0.02

n <- 10^4
h <- 0.001

trueparam <- list(a = a, b = b, c = c, sigma = sigma)
mod <- setModel(drift = c("-a*X", "c*X"),
                diffusion = matrix(c("b", "0", "0", "sigma"), 2, 2),
                solve.variable = c("X", "Y"), state.variable = c("X", "Y"),
                observed.variable = "Y")

samp <- setSampling(delta = h, n = n)
yuima <- simulate(mod, sampling = samp, true.parameter = trueparam)

res <- kalmanBucyFilter(
  yuima, trueparam, mean_init, vcov_init, 0.001,
  are = FALSE, env = globalenv()
)

# vcov and mean slots are accesible by mean and vcov method
mean(res)
vcov(res)
}
