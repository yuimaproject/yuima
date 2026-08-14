\name{yuima.particleFilter-class}
\alias{yuima.particleFilter-class}
\alias{show,yuima.particleFilter-method}
\alias{summary,yuima.particleFilter-method}
\alias{mean,yuima.particleFilter-method}
\alias{vcov,yuima.particleFilter-method}
\title{Particle Filter Result Class}
\description{
The \code{yuima.particleFilter} class contains filtering summaries,
diagnostics, and optionally particle clouds and ancestry information.
Its \code{show} and \code{summary} methods follow the presentation used by
\code{yuima.kalmanBucyFilter} objects.
}
\section{Slots}{
\describe{
  \item{\code{model}}{The fitted \code{yuima.state_space_model}.}
  \item{\code{data}}{The \code{yuima.data} object used by the filter.}
  \item{\code{mean}}{A time-series object of filtered state means.}
  \item{\code{vcov}}{A state-by-state-by-time array of filtered covariance
  matrices.}
  \item{\code{time}}{Observation times represented in the summaries.}
  \item{\code{particles}}{Stored filtering particle clouds, or an empty
  array.}
  \item{\code{weights}}{Weights corresponding to \code{particles}, or an
  empty matrix.}
  \item{\code{paths}}{Euler-Maruyama particle states when requested, or an
  empty array.}
  \item{\code{ancestors}}{One-based parent indices at observation times, or
  an empty matrix.}
  \item{\code{ess}}{Effective sample size at each observation time.}
  \item{\code{branched}}{Logical branching indicator at each observation
  time.}
  \item{\code{intervals}}{Weighted central particle intervals.}
  \item{\code{interval.levels}}{Levels represented by \code{intervals}.}
  \item{\code{logLik}}{Sum of sequential log normalizing constants.}
  \item{\code{logLik.increment}}{Log normalizing increment at each time.}
  \item{\code{call}}{The matched filtering call.}
  \item{\code{settings}}{Particle count, storage mode, and algorithm
  settings.}
}}
\section{Methods}{
\describe{
  \item{\code{show}}{Displays filtered means and covariance matrices.}
  \item{\code{summary}}{Displays model, time-series, storage, and diagnostic
  metadata.}
  \item{\code{mean}}{Returns the filtered mean time series.}
  \item{\code{vcov}}{Returns the filtered covariance array.}
  \item{\code{plot}}{Draws filtered means using the Kalman-Bucy filter plot
  style, with optional stored particle intervals and true state values.}
}}
\seealso{\code{\link{particleFilter}},
\code{\link[=plot,yuima.particleFilter,ANY-method]{plot}},
\code{\link{kalmanBucyFilter}}}
