\name{yuima.data-class}
\docType{class}
\alias{yuima.data-class}

\alias{get.zoo.data,yuima.data-method}
\alias{plot,yuima.data,ANY-method}
%\alias{start,yuima.data-method}
%\alias{end,yuima.data-method}
%\alias{time,yuima.data-method}
\alias{dim,yuima.data-method}
\alias{length,yuima.data-method}
\alias{cce,yuima.data-method}
\alias{llag,yuima.data-method}
\alias{initialize,yuima.data-method}
%\alias{setSampling,yuima.data-method}
\alias{poisson.random.sampling,yuima.data-method}
\alias{subsampling,yuima.data-method}
\alias{cbind.yuima,yuima.data-method}

\title{Class "yuima.data" for the data slot of a "yuima" class object}
\description{
  The \code{yuima.data-class} is a class of the \pkg{yuima} package used to store 
  the data which are hold in the slot \code{data} of an object of
  the \code{yuima-class}.   
  
  Objects from this class contain either true data or simulated data.
}
\section{Slots}{
  \describe{
    \item{\code{original.data}:}{The original data.}
	\item{\code{zoo.data}:}{A list of \code{zoo} format data.}
  }
}
\details{
Objects in this class are created or initialized using the 
methods \code{new} or \code{initialize} or via the function \code{\link{setData}}. 
The preferred way to construct an object in this class is to use the function
\code{\link{setData}}.

Objects in this class are used to store 
the data which are hold in the slot \code{data} of an object of
the \code{yuima-class}.   

Objects in this class contain two slots described here. 
\describe{
  \item{\code{original.data}:}{The slot \code{original.data} contains, as the 
   name suggests, a copy of the original data passed by the user to methods 
   \code{new} or \code{initialize} or to the function \code{\link{setData}}. 
   It is intended for backup purposes.}
  \item{\code{zoo.data}:}{When a new object of this class is created or 
   initialized using the \code{original.data}, the package tries to convert 
   \code{original.data} into an object of class \code{\link{zoo}}. Once 
   coerced to \code{\link{zoo}}, the data are stored in the slot \code{zoo.data}.

   If the conversion fails, the initialization or creation of the object fails.

   Internally, the \pkg{yuima} package stores and operates on  
   \code{\link{zoo}}-type objects.}
}
If data are obtained by simulation, the \code{original.data} slot 
is usually empty.
}
\section{Methods}{
  \describe{
    \item{new}{\code{signature(x = "yuima.data", original.data)}: the function 
	 makes a copy of the prototype object from the class definition of
	 \code{\link{yuima.data-class}}, then
	 calls the \code{initialize} method passing as arguments the newly created 
	 object and the \code{original.data}.} 
    \item{initialize}{\code{signature(x = "yuima.data", original.data)}: makes
	a copy of \code{original.data} into the slot \code{original.data} of 
	\code{x} and tries to coerce \code{original.data} into an object of class
	\code{\link{zoo}}. The result is put in the slot  \code{zoo.data} of 
	\code{x}. If coercion fails, the \code{intialize} method fails as well.}
    \item{get.zoo.data}{\code{signature(x = "yuima.data")}: returns the content
	of the slot \code{zoo.data} of \code{x}.}
    \item{plot}{\code{signature(x = "yuima.data",  \dots)}: calls  
	 \code{\link{plot}} from the \code{\link{zoo}} package with argument
	 \code{x@zoo.data}. Additional arguments \code{\dots} are passed
	 as is to the \code{\link{plot}} function.}
    \item{dim}{\code{signature(x = "yuima.data")}: calls  
	 \code{\link{dim}} from the \code{\link{zoo}} package with argument
	 \code{x@zoo.data}.}
    \item{length}{\code{signature(x = "yuima.data")}: calls  
	 \code{\link{length}} from the \code{\link{zoo}} package with argument
	 \code{x@zoo.data}.}
    \item{cce}{\code{signature(x = "yuima.data")}: calculates asyncronous
	 covariance estimator on the data contained in \code{x@zoo.data}. 
	 For more details see \code{\link{cce}}.}
   \item{llag}{\code{signature(x = "yuima.data")}: calculates lead lag estimate
	  on the data contained in \code{x@zoo.data}. 
	 For more details see \code{\link{llag}}.}   
    \item{cbind.yuima}{\code{signature(x = "yuima.data")}: bind yuima.data object.}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
