\name{setSampling}
\alias{setSampling}
%% \alias{setSampling,ANY-method}
\title{
Set sampling information and create a `sampling' object.
}
\description{
  \code{setSampling} is a constructor for \code{\link{yuima.sampling-class}}.
}
\usage{
  setSampling(Initial = 0, Terminal = 1, division = 100, delta = 0.1, 
   grid = as.numeric(NULL), random = FALSE, scale=1)
}
\arguments{  
  \item{Initial}{Initial time of the grid.}
  \item{Terminal}{Terminal time of the grid.}
  \item{division}{division into this number of intervals.}
  \item{delta}{mesh size in case of regular grid.}
  \item{grid}{a grid of times for the simulation, possibly empty.}
  \item{random}{specify if it is random sampling. See Details.}
  \item{scale}{by default 1. NOT YET IN USE.}
}
\details{
The function creates an object of type 
   \code{\link{yuima.sampling-class}} with several slots.
\describe{
	\item{\code{Initial}:}{initial time of the grid.}
	\item{\code{Terminal}:}{terminal time fo the grid.}
	\item{\code{division}:}{the number of observations - 1.}
	\item{\code{delta}:}{in case of a regular grid is the mesh.}
	\item{\code{grid}:}{the grid of times.}
	\item{\code{random}:}{either \code{FALSE} or the distribution 
	 of the random times.}
	\item{\code{scale}:}{by default 1. NOT YET IN USE.}
	\item{\code{regular}:}{indicator of whether the grid is regular or not.
	 For internal use only.}
}
}
\value{
   An object of type \code{\link{yuima.sampling-class}}.
}
\author{The YUIMA Project Team}
\examples{
samp <- setSampling(Terminal=1, division=1000)
str(samp)
}
\keyword{ts}
