\name{poest}
\alias{poest}
\alias{poest,yuima-method}
\title{P-O Estimator}
\description{P-O(penalized method to ordinary method) estimator using QMLE or adaBayes and \eqn{L^q} penalty.}
\usage{
poest(
  yuima, estimator = "qmle", start, prior=NULL, 
  lower, upper, penalty, method = "L-BFGS-B", ...
  )
}
\arguments{
  \item{yuima}{A 'yuima' object.}
  \item{estimator}{Charactor. 'qmle' or 'adaBayes'.}
  \item{start}{A named list. Initial suggestion for parameter values }
  \item{prior}{A list of prior distributions for the parameters specified by 'code'. See \code{\link{adaBayes}} for details.}
  \item{lower}{A named list for specifying lower bounds of parameters.}
  \item{upper}{A named list for specifying upper bounds of parameters.}
  \item{penalty}{Parameters for penalizing. A list of 2 sub-lists named 'drift' and 'diffusion'. Each sub-list may contain \code{q} , \code{gamma} and \code{r}. See example for details.}
  \item{method}{Optimization method for \code{\link{qmle}} or \code{\link{adaBayes}}.}
  \item{...}{Passed to \code{\link{qmle}} or \code{\link{adaBayes}} method.}
}
\details{
  Calculate the P-O estimator for stochastic processes by using  the QMLE or adaBayes function and \eqn{L^q} penalty.
}
\value{
  \item{initial_estimator}{Result of the initial QMLE or adaBayes. }
  \item{final_estimator}{Result of the final QMLE or adaBayes.}
  \item{model}{Given yuima model.}
  \item{estimator}{Name of estimation method.}
  \item{call}{Call}
  \item{coef}{Result of estimation.}
  \item{vcov}{Vcov matrix of coefficients.}
  \item{selected_coef}{Coefficients not regarded to be 0 in the penalization step.}
}
\author{Kosuke Kito and Masahiro Kurisaki with YUIMA project Team}
\references{
Suzuki, T. and Yoshida, N., 2020. Penalized least squares approximation methods and their applications to stochastic processes. Japanese Journal of Statistics and Data Science, 3(2), pp.513-541.

Huang, J., Horowitz, J. and Ma, S., 2008. Asymptotic properties of bridge estimators in sparse high-dimensional regression models. The Annals of Statistics, 36(2).
}
\examples{

  ### example1 normal case
  # model definition
  d <- 5
  p <- d
  sol <- c(paste0("x", 1:d), "y")
  drift <- c(paste0("-0.2*", sol[1:d]), 0)
  logsigma <- paste(paste0("theta", 1:p, "*", sol[1:p]), collapse = "+")
  diffusion_Y <- paste0("pmin(exp(", logsigma, "), 10^5)")
  diffusion <- diag(0.5, d + 1)
  diffusion[d+1,d+1] <- diffusion_Y

  ymodel <- setModel(drift = drift, diffusion = diffusion, 
                    solve.variable = sol, state.variable = sol)
  n <- 1000
  ysamp <- setSampling(Terminal = 1, n = n)
  yuima <- setYuima(model = ymodel, sampling = ysamp)

  true.par <- as.list(c(1, -1, double(p - 2)))
  names(true.par) <- paste0("theta", 1:p)

  set.seed(111)
  yuima <- simulate(yuima, xinit = 1, true.parameter = true.par)

  # other arguments definition
  start.par <- as.list(runif(p,-2,2))
  names(start.par) <- paste0("theta", 1:p)

  lower.par <- as.list(rep(-2, p))
  names(lower.par) <- paste0("theta", 1:p)
  upper.par <- as.list(rep(2, p))
  names(upper.par) <- paste0("theta", 1:p)

  #penalty parameters
  drift.penalty <- list(q = 0.3, gamma = 1,r = 1.2)
  diffusion.penalty <- list(q = 0.3, gamma=1,r = 1.2)
  penalty <- list(drift = drift.penalty, diffusion = diffusion.penalty)

  method = "L-BFGS-B"

  # estimation
  \donttest{
  res.po <- poest(
    yuima = yuima, start = start.par, 
    lower = lower.par, upper = upper.par, 
    penalty = penalty, rcpp = TRUE
    )
  

  # show results
  print(res.po)
  summary(res.po)

  }


  ### example2 not giving penalty
  # model definition
  d <- 5
  p <- d
  sol <- c(paste0("x", 1:d), "y")
  drift <- c(paste0("-0.2*", sol[1:d]), 0)
  logsigma <- paste(paste0("theta", 1:p, "*", sol[1:p]), collapse = "+")
  diffusion_Y <- paste0("pmin(exp(", logsigma, "), 10^5)")
  diffusion <- diag(0.5, d + 1)
  diffusion[d+1,d+1] <- diffusion_Y

  ymodel <- setModel(drift = drift, diffusion = diffusion, 
                    solve.variable = sol, state.variable = sol)
  n <- 1000
  ysamp <- setSampling(Terminal = 1, n = n)
  yuima <- setYuima(model = ymodel, sampling = ysamp)

  true.par <- as.list(c(1, -1, double(p - 2)))
  names(true.par) <- paste0("theta", 1:p)

  set.seed(111)
  yuima <- simulate(yuima, xinit = 1, true.parameter = true.par)

  # other arguments definition
  start.par <- as.list(runif(p,-2,2))
  names(start.par) <- paste0("theta", 1:p)

  lower.par <- as.list(rep(-2, p))
  names(lower.par) <- paste0("theta", 1:p)
  upper.par <- as.list(rep(2, p))
  names(upper.par) <- paste0("theta", 1:p)
  method = "L-BFGS-B"

  # estimation
  \donttest{
  res.po <- poest(
    yuima = yuima, start = start.par, 
    lower = lower.par, upper = upper.par, rcpp = TRUE
    )
  # when penalty or any element of penalty is not given, 
  # the following values will be used for elements not given.
  # list(
  #   drift     = list(q = 1, gamma = 1, r = (3 - q + gamma) / 2),
  #   diffusion = list(q = 1, gamma = 1, r = (3 - q + gamma) / 2),
  # )

  # show results
  print(res.po)
  summary(res.po)
  }

  ### example3 use adaBayes
  # model definition
  d <- 5
  p <- d
  sol <- c(paste0("x", 1:d), "y")
  drift <- c(paste0("-0.2*", sol[1:d]), 0)
  logsigma <- paste(paste0("theta", 1:p, "*", sol[1:p]), collapse = "+")
  diffusion_Y <- paste0("pmin(exp(", logsigma, "), 10^5)")
  diffusion <- diag(0.5, d + 1)
  diffusion[d+1,d+1] <- diffusion_Y

  ymodel <- setModel(drift = drift, diffusion = diffusion, 
                  solve.variable = sol, state.variable = sol)
  n <- 1000
  ysamp <- setSampling(Terminal = 1, n = n)
  yuima <- setYuima(model = ymodel, sampling = ysamp)

  true.par <- as.list(c(1, -1, double(p - 2)))
  names(true.par) <- paste0("theta", 1:p)

  set.seed(111)
  yuima <- simulate(yuima, xinit = 1, true.parameter = true.par)

  # other arguments definition
  estimator <- "adaBayes"
  start.par <- as.list(runif(p,-2,2))
  names(start.par) <- paste0("theta", 1:p)

  prior <- as.list(rep(list(list(measure.type="code",df="dunif(z,-2,2)")), p))
  names(prior) <- paste0("theta", 1:p)

  lower.par <- as.list(rep(-2, p))
  names(lower.par) <- paste0("theta", 1:p)
  upper.par <- as.list(rep(2, p))
  names(upper.par) <- paste0("theta", 1:p)

  #penalty parameters
  drift.penalty <- list(q = 0.3, gamma = 1,r = 1.2)
  diffusion.penalty <- list(q = 0.3, gamma=1,r = 1.2)
  penalty <- list(drift = drift.penalty, diffusion = diffusion.penalty)

  method <- "mcmc"
  mcmc <- 1000
  path <- TRUE

  # estimation
  \donttest{
  res.po <- poest(
    yuima = yuima, estimator = estimator, prior = prior, 
    start = start.par, lower = lower.par, upper = upper.par, 
    penalty = penalty, method = method, 
    mcmc = mcmc, path = path, rcpp = TRUE)

  # show results
  print(res.po)
  summary(res.po)
  }
}
