\name{yuima.Map-class}
\docType{class}
\alias{yuima.Map-class}
\alias{yuima.Map}
%\alias{info.Map}
%\alias{param.Map}
\alias{initialize,yuima.Map-method}
\alias{simulate,yuima.Map-method}


\title{Class for the mathematical description of function of a stochastic process}
\description{
  The \code{yuima.Map} class is a class of the  \pkg{yuima} package that extends the \code{\link{yuima-class}} it represents a map of a stochastic process

  \code{ zt = g(theta, Xt, t) : R^{q x d x 1} -> R^{l1 x l2 x ...}}

or an operator between two independent stochasic process:

  \code{ zt = h(theta, Xt, Yt, t)}

where \code{Xt} and \code{Yt}  are object of class \code{\link{yuima.model-class}} or \code{\link{yuima-class}} with the same dimension.
}

\section{Slots}{
Here we report the additional slots of an object of class \code{yuima.Map} with respect to the \code{\link{yuima-class}}:
  \describe{
    \item{\code{Output}:}{It is an object of class \code{info.Map} and it is composed by the following slots:
    \describe{
        \item{\code{formula}:}{It is a \code{vector} that contains the components of map \code{g(theta, Xt, t)} or the operator \code{h(theta, Xt, Yt, t)} }
        \item{\code{dimension}:}{a \code{numeric} object that is the dimensions of the Map.}
        \item{\code{type}:}{If \code{type = "Maps"}, the Map is a map of stochastic process, If \code{type = "Operator"}, the result is an operator between two independent stochastic process}
       \item{\code{param}}{it is an object of class \code{param.Map} and it is composed by the following slots:
       \describe{
          \item{\code{out.var}:}{labels for Map.}
          \item{\code{allparam}:}{labels of all parameters (model and  map/operators).}
          \item{\code{allparamMap}:}{labels of map/operator parameters.}
          \item{\code{common}:}{common parameters.}
          \item{\code{Input.var}:}{labels for inputs.}
          \item{\code{time.var}:}{label for time variable.}
        }
       }
      }
    }
  }
}

\section{Methods}{
  \describe{
    \item{simulate}{simulation method. For more information see
	  \code{\link{simulate}}.}
%    \item{toLatex}{This method converts an object of \code{yuima.poisson-class} to character vectors with LaTeX markup.}
%    \item{qmle}{Quasi maximum likelihood estimation procedure. For more information see \code{\link{qmle}}. }
  }
}
\author{The YUIMA Project Team}
\keyword{classes}
