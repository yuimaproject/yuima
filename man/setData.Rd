\name{setData}
\alias{setData}
\alias{get.zoo.data}
\alias{dim}
\alias{length}
\alias{cbind.yuima}

\title{
Set and access data of an object of type "yuima.data" or "yuima".
}
\description{
\code{setData} constructs an object of \code{\link{yuima.data-class}}.

\code{get.zoo.data} returns the content of the \code{zoo.data} slot of a
\code{\link{yuima.data-class}} object. (Note: value is a \code{list} of 
\code{\link{zoo}} objects).

\code{plot} plot method for object of \code{\link{yuima.data-class}} or 
\code{\link{yuima-class}}.

\code{dim} returns the \code{\link{dim}} of the \code{zoo.data} slot of a
\code{\link{yuima.data-class}} object.

\code{length} returns the \code{\link{length}} of the time series in 
\code{zoo.data} slot of a \code{\link{yuima.data-class}} object.
\code{cbind.yuima} bind yuima.data object.
}
\usage{
  setData(original.data, delta=NULL, t0=0)
  get.zoo.data(x)
}
\arguments{  
  \item{original.data}{some type of data, usually some sort of time series. 
  The function always tries to convert to the input data into an object of
   \code{\link{zoo}}-type. See Details.}
  \item{x}{an object of type \code{\link{yuima.data-class}} or 
  \code{\link{yuima-class}}.}
  \item{delta}{If there is the need to redefine on the fly the \code{delta} increment of the data to make it consistent to statistical theory. See Details.}
\item{t0}{the time origin for the internal \code{zoo.data} slot, defaults to 0.}
}
\details{
Objects in the  \code{\link{yuima.data-class}} contain two slots:
\describe{
  \item{\code{original.data}:}{The slot \code{original.data} contains, as the 
   name suggests, a copy of the original data passed to the function 
   \code{setData}. It is intended for backup purposes.}
  \item{\code{zoo.data}:}{the function \code{setData} tries to convert 
   \code{original.data} into an object of class \code{\link{zoo}}. The 
   coerced \code{\link{zoo}} data are stored in the slot \code{zoo.data}.
   If the conversion fails the function exits with an error.
   Internally, the \pkg{yuima} package stores and operates on  
   \code{\link{zoo}}-type objects.}
}

The function \code{get.zoo.data}
returns the content of the slot \code{zoo.data} of \code{x} if \code{x}
is of \code{\link{yuima.data-class}} or the content of 
\code{x@data@zoo.data} if \code{x} is of \code{\link{yuima-class}}.
}
\value{
   \item{value}{a list of object(s) of \code{\link{yuima.data-class}} for 
   \code{setData}. 
   The content of the \code{zoo.data}  slot for \code{get.zoo.data}}
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

# changing delta on the fly to 1/252
mysde2 <- setYuima(data=setData(X, delta=1/252))
str(get.zoo.data(mysde2))
plot(mysde2)
dim(mysde2)
length(mysde2)

# changing delta on the fly to 1/252 and shifting time to t0=1
mysde2 <- setYuima(data=setData(X, delta=1/252, t0=1))
str(get.zoo.data(mysde2))
plot(mysde2)
dim(mysde2)
length(mysde2)


}
\keyword{ts}
