\name{setData}
\alias{setData}
\alias{get.zoo.data}
\alias{dim}
\alias{length}

\title{
Set and access data of an `yuima' or `yuima.data' object.
}
\description{
\code{setData} is a constructor for \code{yuima.data} class.
\code{get.zoo.data} returns the \code{zoo.data} field of an \code{yuima.data} object.(warning: value is LIST of zoo's)
\code{plot} plots the \code{yuima.data} object.
\code{dim} returns \code{dim} of the \code{zoo.data} field.
\code{length} returns the \code{length} of the time series in \code{zoo.data} field.
}
\usage{
  setData(original.data)
  get.zoo.data(x)
}
\arguments{  
  \item{original.data}{Some type of data.}
  \item{x}{An object of class \code{yuima.data} or \code{yuima}.}
}
\details{
    class \code{yuima.data} has two slots,
    \code{zoo.data} is a list of \code{zoo} version of data, and
    \code{original.data} is a copy of the original data
}
\value{
   \item{value}{A list of object(s) of class \code{yuima.data} for \code{setData}. 
   The \code{zoo.data} field for \code{get.zoo.data}}
}
\author{YUIMA Project Team}
\examples{
library(zoo)
X <- ts(matrix(rnorm(200),100,2))
mydata <- setData(X)
str(get.zoo.data(mydata))
#plot(get.zoo.data(mydata))
dim(mydata)
length(mydata)
plot(mydata)

# exactly the same output
mysde <- setYuima(data=setData(X))
str(get.zoo.data(mysde))
plot(mysde)
dim(mysde)
length(mysde)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
