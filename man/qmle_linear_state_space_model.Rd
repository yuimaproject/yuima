\name{qmle.linear_state_space_model}
\alias{qmle.linear_state_space_model}
\alias{qmle.linear_state_space_model,yuima-method}
\title{
Calculate quasi-likelihood and ML estimator of the least squares estimator for Linear State Space Model
}
\description{
Linear State Space Model version of \code{\link{qmle}}.
}
\usage{
qmle.linear_state_space_model(yuima, start, lower, upper, method = "L-BFGS-B", fixed = list(), envir = globalenv(), filter_mean_init, explicit = FALSE, rcpp = FALSE, drop_terms = 0, ...)
}
\arguments{
  \item{yuima}{A yuima object. The class of \code{yuima@model} must be `yuima.linear_state_space_model`.}
  \item{start}{Initial values to be passed to the optimizer.}
  \item{lower}{A named list specifying lower bounds of parameters.}
  \item{upper}{A named list specifying upper bounds of parameters.}
  \item{method}{Optimization method. See \code{\link{optim}}.
  \item{fixed}{For conditional (quasi)maximum likelihood estimation.}
  \item{envir}{An environment where the model coefficients are evaluated.}
  \item{filter_mean_init}{Initial values of unobserved variables in filter calculation.}
  \item{explicit}{A logical value. If \code{TRUE}, use the algebraic Riccati equation and calculate filter explicitly.}
  \item{rcpp}{A logical value. If \code{TRUE}, the algebraic Riccati equation is solved by Rcpp.}
  \item{drop_terms}{Losses for the first given number of terms will be ignored during optimization.}
  \item{...}{Passed to \code{\link{optim}} method. See Examples.}
}
\value{
A \code{\link{yuima.linear_state_space_qmle}} object, extending the \code{\link{yuima.qmle}} class.
\item{model}{A \code{\link{yuima.linear_state_space_model}} object.}
\item{drop_terms}{Numeric.}
\item{...}{Slots in \code{yuima.qmle}.}
}
\references{
## Estimation of linear state space model

Kurisaki, M. Parameter estimation for ergodic linear SDEs from partial and discrete observations. Stat Inference Stoch Process 26, 279–330 (2023).

## Algebraic Riccati equation

Bini, D., Iannazzo, B., &amp; Meini, B. (2012). Numerical solution of algebraic Riccati equations. Society for Industrial and Applied Mathematics. 
}
\author{YUIMA TEAM}

\examples{
### Set model
drift = c("a*X", "X")
diffusion = matrix(c("b", "0", "0", "sigma"), nrow = 2)
# Decide to use `model.class="linearStateSpaceModel"` implicitly.
ymodel = setModel(
    drift = drift, 
    diffusion = diffusion, 
    solve.variable = c("X", "Y"),
    state.variable = c("X", "Y"),
    observed.variable = "Y",
    model.class = "linearStateSpaceModel"
)

### Set data
T <- 100
N <- 50000
n <- N
h <- T/N

true.par = list(
    a = -1.5,
    b = 0.3,
    sigma = 0.053
)
tmp.yuima <- simulate(ymodel, true.parameter = true.par, sampling = setSampling(n = N, Terminal = T))
ydata <- tmp.yuima@data
rm(tmp.yuima)

### Set yuima
variable_data_mapping <- list(
    "X" = NA,
    "Y" = 2
)
yuima <- setYuima(model = ymodel, data = ydata, variable_data_mapping = variable_data_mapping)

# Estimate
upper.par = list(
    a = 1,
    b = 5,
    sigma = 1
)

lower.par = list(
    a = -10,
    b = 0.01,
    sigma = 0.001
)

start.par = list(
    a = 0.5,
    b = 4,
    sigma = 0.9
)

filter_mean_init = c(0)

res <- qmle.linear_state_space_model(yuima, start = start.par, upper = upper.par, lower = lower.par, filter_mean_init = filter_mean_init, explicit = FALSE, rcpp = TRUE)

}
