\name{setYuima}
\alias{setYuima}
\title{
Create a 'yuima' object by combining 'model', 'data', 'sampling' and 
'characteristic' slots.
}
\description{
  \code{setYuima} is a constructor for \code{\link{yuima-class}}.
}
\usage{
setYuima(data, model, sampling, characteristic)
}
\arguments{
  \item{data}{an object of class \code{\link{yuima.data-class}}.}
  \item{model}{an object of class \code{\link{yuima.model-class}}.}
  \item{sampling}{an object of class \code{\link{yuima.sampling-class}}.}
  \item{characteristic}{an object of class \code{\link{yuima.characteristic-class}}.}
}
\details{
The \code{yuima-class} object is the main object of the \pkg{yuima} package.
Some of the slots may be missing. 

The \code{data} slot contain the data, either empirical or simulated. 

The \code{model} contains the description of the 
(statistical) model which is used to generate the \code{data} via different 
simulation schemes, to draw inference from the \code{data} or both.

The \code{sampling} slot contains information on how the \code{data} have been
collected or how they should be generated.

The slot \code{characteristic} contains information on \code{PLEASE FINISH THIS}.
}
\value{
  An object of type \code{\link{yuima-class}}.
}
\author{The YUIMA Project Team}

\keyword{ts}
