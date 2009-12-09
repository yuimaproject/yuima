\name{setCharacteristic}
\alias{setCharacteristic}
%% \alias{setSampling,ANY-method}
\title{
  Set characteristic information and create a `characteristic' object.
}
\description{
  \code{setCharacteristic} is a constructor for \code{characteristic} class.
}
\usage{
  setCharacteristic(equation.number,time.scale)
}
\arguments{  
  \item{equation.number}{The number of equations modeled in \code{yuima} object.}
  \item{time.scale}{time.scale assumed in the model.}
}
\details{
  class \code{characteristic} has two slots,
  \code{equation.number} is the number of equations handled in the yuima object, and 
  \code{time.scale} is a hoge of  \code{characteristic}.
}
\value{
   An object of class \code{characteristic}.
}
\author{The YUIMA Project Team}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
