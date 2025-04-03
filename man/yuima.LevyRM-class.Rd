\name{yuima.LevyRM-class}
\docType{class}
\alias{yuima.LevyRM-class}
\alias{yuima.LevyRM}
\alias{initialize,yuima.LevyRM-method}
\alias{simulate,yuima.LevyRM-method}
%\alias{initialize,yuima.LevyRM-method}
\alias{estimation_RLM,yuima.LevyRM-function}

\title{\code{yuima.LevyRM}: A Class for the Mathematical Description of the t-Student Regression Model}

\description{A yuima class that contains all information on the regression model with t-student Levy process noise. This class extends \code{\link{yuima-class}} and contains information on the regressors used in the definition of the model. The regressors are represented by an object of \code{\link{yuima.model-class}}.

An object of this class can be created by calls of the function \code{\link{setLRM}}.
}

\section{Methods}{
  \describe{
    \item{initialize}{Initialize method. It makes a copy of each argument.}
  	\item{simulate}{simulation method. For more information see
	  \code{\link{simulate}}.}
  }
}


\author{The YUIMA Project Team

Contacts: Lorenzo Mercuri \email{lorenzo.mercuri@unimi.it}
}
\keyword{classes}


