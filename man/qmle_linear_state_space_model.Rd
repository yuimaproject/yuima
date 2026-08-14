\name{qmle.linear_state_space_model}
\alias{qmle.linear_state_space_model}
\alias{qmle.linear_state_space_model,yuima-method}
\alias{yuima.linear_state_space_qmle-class}
\title{
Calculate Quasi-Likelihood and Maximum Likelihood Estimator for Linear State Space Model
}
\description{
A function for the quasi-maximum likelihood estimation of linear state space models, extending the \code{\link{qmle}} function.
}
\usage{
qmle.linear_state_space_model(
  yuima, start, lower, upper, method = "L-BFGS-B", fixed = list(),
  envir = globalenv(), filter_mean_init, explicit = FALSE, drop_terms = 0,
  ...
)
}
\arguments{
  \item{yuima}{A yuima object. The class of \code{yuima@model} must be \code{yuima.linear_state_space_model}.}
  \item{start}{Initial values for the parameters to be passed to the optimizer.}
  \item{lower}{A named list specifying the lower bounds of the parameters.}
  \item{upper}{A named list specifying the upper bounds of the parameters.}
  \item{method}{The optimization method to be used. See \code{\link{optim}}.}
  \item{fixed}{A named list of parameters to be held fixed during optimization.}
  \item{envir}{An environment in which the model coefficients are evaluated.}
  \item{filter_mean_init}{Initial values of unobserved variables for the filter calculation.}
  \item{explicit}{A logical value. If \code{TRUE}, use the explicit formula to solve the filitering equation.}
  \item{drop_terms}{A numeric value. The specified number of initial terms in the filtering result will be excluded from the quasi-likelihood function calculation.}
  \item{...}{Additional arguments to be passed to the \code{\link{optim}} method. See Examples.}
}
\value{
A \code{yuima.linear_state_space_qmle-class} object, extending the \code{\link{yuima.qmle-class}} class.
\item{model}{A \code{\link{yuima.linear_state_space_model}} object.}
\item{drop_terms}{A numeric value representing the number of terms ignored during optimization.}
\item{explicit}{A logical value provided as an argument.}
\item{mean_init}{A numeric value provided as the argument \code{filter_mean_init}.}
\item{...}{Additional slots inherited from the \code{yuima.qmle} class.}
}
\references{
Kurisaki, M. (2023). Parameter estimation for ergodic linear SDEs from partial and discrete observations. Stat Inference Stoch Process, 26, 279-330.

Bini, D., Iannazzo, B., & Meini, B. (2012). Numerical Solution of Algebraic Riccati Equations. Society for Industrial and Applied Mathematics.
}
\author{YUIMA TEAM}
\examples{
### Set model
drift <- c("a*X", "X")
diffusion <- c("b", "sigma")
# Use `model.class="linearStateSpaceModel"` implicitly.
ymodel <- setModel(
    drift = drift, 
    diffusion = diffusion, 
    solve.variable = c("X", "Y"),
    state.variable = c("X", "Y"),
    observed.variable = "Y"
)

### Set data
T <- 100
N <- 50000
n <- N
h <- T / N

true.par <- list(
    a = -1.5,
    b = 0.3,
    sigma = 0.053
)
tmp.yuima <- simulate(
  ymodel, true.parameter = true.par, sampling = setSampling(n = N, Terminal = T)
)
ydata <- tmp.yuima@data
rm(tmp.yuima)

### Set yuima
variable_data_mapping <- list(
    #"X" = NA,
    "Y" = 2
)
yuima <- setYuima(model = ymodel, data = ydata, variable_data_mapping = variable_data_mapping)

# Estimate
upper.par <- list(
    a = 1,
    b = 5,
    sigma = 1
)

lower.par <- list(
    a = -10,
    b = 0.01,
    sigma = 0.001
)

start.par <- list(
    a = 0.5,
    b = 4,
    sigma = 0.9
)

filter_mean_init <- 0

res <- qmle.linear_state_space_model(
  yuima, start = start.par, upper = upper.par, lower = lower.par, 
  filter_mean_init = filter_mean_init
)
}
