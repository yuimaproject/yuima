\name{setSampling}
\alias{setSampling}
%% \alias{setSampling,ANY-method}
\title{
Set sampling information and create a `sampling' object.
}
\description{
  \code{setSampling} is a constructor for \code{sampling} class.
}
\usage{
  setSampling(Initial = 0, Terminal = 1, division = 100, delta =
                 0.1, grid = as.numeric(NULL), random = FALSE)
}
\arguments{  
  \item{Initial}{Initial}
  \item{Terminal}{Terminal}
  \item{division}{division}
  \item{delta}{mesh size in case of regular grid}
  \item{grid}{a grid of times for the simulation, possibly empty}
  \item{random}{specify if it is random sampling}
}
\details{
  class \code{sampling} has two slots,
  \code{Terminal} is a hoge of \code{sampling}, and
  \code{division} is a hoge of  \code{sampling}	
}
\value{
   An object of class \code{sampling}.
}
\author{YUIMA Project Team}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
