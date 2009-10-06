\name{setYuima}
\alias{setYuima}
\title{
Create a 'yuima' object by combining 'model', 'data', 'sampling' and 'characteristic' object.
}
\description{
  \code{setYuima} is a constructor for \code{yuima} class.
}
\usage{
setYuima(data,model,sampling,characteristic)
}
\arguments{
  \item{data}{an object of class \code{data}}
  \item{model}{an object of class \code{model}}
  \item{sampling}{an object of class \code{sampling}}
  \item{characteristic}{an object of class \code{characteristic}}
}
\details{
  class \code{yuima} has four slots,
  \code{data} is an object of class \code{data},
  \code{model} is an object of class \code{model},
  \code{sampling} is an object of class \code{sampling}, and
  \code{characteristic} is an object of class \code{characteristic}.
}
\value{
  An object of class \code{yuima} for \code{setYuima}.
}
\author{YUIMA Project Team}

\keyword{ts}
