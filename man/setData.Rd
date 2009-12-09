\name{setData}
\alias{setData}
\alias{get.zoo.data}
\alias{dim}
\alias{length}

\title{
Set and access data of an 'yuima-class' or 
'yuima.data-class' object.
}
\description{
\code{setData} is a constructor for \code{\link{yuima.data-class}}.

\code{get.zoo.data} returns the \code{zoo.data} slot of a
\code{\link{yuima.data-class}} object. (Note: value is a \code{list} of 
\code{\link{zoo}} objects).

\code{plot} plots the \code{\link{yuima.data-class}} or 
\code{\link{yuima-class}} objects.

\code{dim} returns the \code{\link{dim}} of the \code{zoo.data} field.

\code{length} returns the \code{\link{length}} of the time series in 
\code{zoo.data} field.
}
\usage{
  setData(original.data)
  get.zoo.data(x)
}
\arguments{  
  \item{original.data}{Some type of data, usually some sort of time series. 
  The function always tries to convert to the input data into an object of
   \code{\link{zoo}}-type.}
  \item{x}{An object of type \code{\link{yuima.data-class}} or 
  \code{\link{yuima-class}}.}
}
\details{
The \code{\link{yuima.data-class}} contains two slots. 

The slot \code{original.data} contains, as the name suggests,
a copy of the original data passed by the user and it is intended for backup 
purposes.

Internally, the \pkg{yuima} package stores and operates on  
\code{\link{zoo}}-type objects.
The \code{\link{zoo}}-fied data are contained in the slot \code{zoo.data}.
}
\value{
   \item{value}{A list of object(s) of class \code{\link{yuima.data-class}} for 
   \code{setData}. 
   The \code{zoo.data} field for \code{get.zoo.data}}
}
\author{The YUIMA Project Team}
\examples{
X <- ts(matrix(rnorm(200),100,2))
mydata <- setData(X)
str(get.zoo.data(mydata))
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
\keyword{ts}
