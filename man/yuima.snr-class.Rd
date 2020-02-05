\name{yuima.snr-class}
\docType{class}
\alias{yuima.snr-class}

\alias{show,yuima.snr-method}

\title{Class "yuima.snr" for self-normalized residuals of SDE "yuima" class object}
\description{
  The \code{yuima.snr-class} is a class of the \pkg{yuima} package used to store 
  the calculatied self-normalized residuals of an SDEs.   
  
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{The original call.}
	\item{\code{coef}:}{A \code{numeric} vector.}
	\item{\code{snr}:}{A  \code{numeric} vector of residuals.}
	\item{\code{model}:}{A \code{yuima.model}.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{print method}
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
