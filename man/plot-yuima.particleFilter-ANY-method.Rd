\name{plot,yuima.particleFilter,ANY-method}
\alias{plot,yuima.particleFilter,ANY-method}
\title{Plotting Method for Particle Filter}
\usage{
\S4method{plot}{yuima.particleFilter,ANY}(x, plot_truth = FALSE, level = 0)
}
\arguments{
\item{x}{A \code{\link{yuima.particleFilter-class}} object.}

\item{plot_truth}{Logical. If \code{TRUE}, plot true values of state
variables when they are present in the data.}

\item{level}{Numeric. If \code{0 < level < 1}, plot the stored weighted
particle interval of that level. Set \code{level = 0} to omit intervals.}
}
\value{
\code{NULL}, invisibly. A plot is drawn on the active graphics device.
}
\description{
Plots filtered state estimates from a \code{yuima.particleFilter} object.
Its panel layout, colors, line types, labels, truth overlay, and interval band
follow the plotting method for \code{yuima.kalmanBucyFilter}.
}
\details{
Particle intervals are weighted empirical quantiles computed when
\code{particleFilter()} is called. A requested \code{level} must match one of
the values supplied through \code{interval_levels}; otherwise an error lists
the available levels. The method does not substitute a Gaussian interval from
the variance-covariance array.

All writable graphical parameters are restored when the method exits,
including when plotting stops with an error. Thus calling this method does not
leave changes to \code{par()} in the active graphics device.
}
\examples{
\dontrun{
model <- setModel(
  drift = c("-a * X", "X"),
  diffusion = matrix(c("b", "0", "0", "sigma"), 2, 2, byrow = TRUE),
  solve.variable = c("X", "Y"),
  state.variable = c("X", "Y"),
  observed.variable = "Y"
)
sampling <- setSampling(delta = 0.01, n = 100)
parameters <- list(a = 0.7, b = 0.3, sigma = 0.4)
object <- simulate(
  model,
  sampling = sampling,
  xinit = c(0, 0),
  true.parameter = parameters
)
result <- particleFilter(
  object,
  xinits = matrix(rnorm(1000, sd = 0.5), ncol = 1),
  params = parameters,
  interval_levels = 0.95
)
plot(result, plot_truth = TRUE, level = 0.95)
}
}
\author{
The YUIMA Project Team
}
