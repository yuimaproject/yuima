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
  setSampling(Initial = 0, Terminal = 1, n = 100, delta, 
   grid, random = FALSE, sdelta=as.numeric(NULL), 
   sgrid=as.numeric(NULL), interpolation="pt" )
}
\arguments{  
  \item{Initial}{Initial time of the grid.}
  \item{Terminal}{Terminal time of the grid.}
  \item{n}{number of time intervals.}
  \item{delta}{mesh size in case of regular time grid.}
  \item{grid}{a grid of times for the simulation, possibly empty.}
  \item{random}{specify if it is random sampling. See Details.}
  \item{sdelta}{mesh size in case of regular space grid.}
  \item{sgrid}{a grid in space for the simulation, possibly empty.}
  \item{interpolation}{a rule of interpolation in case of subsampling. 
   By default, the previous tick interpolation. See Details.}
}
\details{
The function creates an object of type 
   \code{\link{yuima.sampling-class}} with several slots.
\describe{
	\item{\code{Initial}:}{initial time of the grid.}
	\item{\code{Terminal}:}{terminal time fo the grid.}
	\item{\code{n}:}{the number of observations - 1.}
	\item{\code{delta}:}{in case of a regular time grid it is the mesh.}
	\item{\code{grid}:}{the grid of times.}
	\item{\code{random}:}{either \code{FALSE} or the distribution 
	 of the random times.}
	\item{\code{regular}:}{indicator of whether the grid is regular or not.
	 For internal use only.}
	\item{\code{sdelta}:}{in case of a regular space grid it is the mesh.}
	\item{\code{sgrid}:}{the grid in space.}
	\item{\code{oindex}:}{in case of interpolation, a vector of indexes 
	 corresponding to the original observations used for the approximation.}
	\item{\code{interpolation}:}{the name of the interpolation method used.}
}

In case of subsampling, the observations are subsampled on some given 
\code{grid}/\code{sgrid} or according to some \code{random} times. When 
the original observations do not exist at a give point of the grid they are
obtained by some approximation method. Available methods are \code{"pt"} or
\code{"previous tick"} observation method, \code{"nt"} or \code{"next tick"} 
observation method, or by l\code{"linear"} interpolation.
In case of interpolation, the slot \code{oindex} contains the vector of indexes 
corresponding to the original observations used for the approximation. For the
linear method the index corresponds to the left-most observation.

The slot \code{random} is used as information in case a \code{grid} is
already determined (e.g. \code{n} or \code{delta}, etc. ot the \code{grid}
itself are given) or if some subsampling has occurred or if some particular
method which causes a random grid is used in simulation (for example the
space discretized Euler scheme). The slot \code{random} contains a list
of two elements \code{distr} and \code{scale}, where \code{distr} is a
the distribution of independent random times and \code{scale} is either a 
scaling constant or a scaling function.
If the \code{grid} of times is deterministic, then \code{random} is \code{FALSE}.

If not specified and \code{random=FALSE}, the slot \code{grid} is filled 
automatically by the function. It is eventually modified or created
after the call to the function \code{\link{simulate}}.

If \code{delta} is not specified, it is calculated as \code{(Terminal-Initial)/n)}.
If \code{delta} is specified, the \code{Terminal} is adjusted to be equal to
\code{Initial+n*delta}.

The vectors \code{delta}, \code{n}, \code{Initial} and \code{Terminal} may
have different lengths, but then they are extended to the maximal length to
keep consistency. See examples.

If \code{grid} is specified, it takes precedence over all other arguments.
}
\value{
   An object of type \code{\link{yuima.sampling-class}}.
}
\author{The YUIMA Project Team}
\examples{
samp <- setSampling(Terminal=1, n=1000)
str(samp)

samp <- setSampling(Terminal=1, n=1000, delta=0.3)
str(samp)


samp <- setSampling(Terminal=1, n=1000, delta=c(0.1,0.3))
str(samp)

samp <- setSampling(Terminal=1:3, n=1000)
str(samp)
}
\keyword{ts}
