\name{yuima.Output-class}
\docType{class}
\alias{yuima.Output-class}
\alias{yuima.Output}
%\alias{info.Output}
%\alias{param.Output}
\alias{initialize,yuima.Output-method}
\alias{simulate,yuima.Output-method}


\title{Class for the mathematical description of function of a stochastic process}
\description{
  The \code{yuima.Output} class is a class of the  \pkg{yuima} package that extends the \code{\link{yuima-class}} it represents a map of a stochastic process

  \code{ zt = g(theta, Xt, t) : R^{q x d x 1} -> R^{l1 x l2 x ...}}

or an operator between two independent stochasic process:

  \code{ zt = h(theta, Xt, Yt, t)}

where \code{Xt} and \code{Yt}  areobject of class \code{\link{yuima.model-class}} or \code{\link{yuima-class}} with the same dimension.
}

\section{Slots}{
Here we report the additional slots of an object of class \code{yuima.Output} with respect to the \code{\link{yuima-class}}:
  \describe{
    \item{\code{Output}:}{It is an object of class \code{info.Output} and it is composed by the following slots:
    \describe{
        \item{\code{formula}:}{It is a \code{list} that contains the components of map \code{g(theta, Xt, t)} or the operator \code{h(theta, Xt, Yt, t)} }
        \item{\code{dimension}:}{a \code{numeric} object that is the dimensions of the output.}
        \item{\code{type}:}{If \code{type = "Maps"}, the output is a map of stochastic process, If \code{type = "Operator"}, the result is an operator between two independent stochastic process}
       \item{\code{param}}{it is an object of class \code{param.Output} and it is composed by the following slots:
       \describe{
          \item{\code{out.var}:}{labels for output.}
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
