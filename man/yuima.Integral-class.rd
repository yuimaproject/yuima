\name{yuima.Integral-class}
\docType{class}
\alias{yuima.Integral-class}
\alias{yuima.Integral}
%\alias{info.Output}
%\alias{param.Output}
\alias{initialize,yuima.Integral-method}
\alias{simulate,yuima.Integral-method}


\title{Class for the mathematical description of integral of a stochastic process}

\description{
  The \code{yuima.Integral} class is a class of the  \pkg{yuima} package that extends the \code{\link{yuima-class}} it represents a integral of a stochastic process

  \code{ zt = int^{t}_0 h(theta, Xs, s) dXs}

}

\section{Slots}{
In the following we report the the additional slots of an object of class \code{yuima.Integral} with respect to the \code{\link{yuima-class}}:
  \describe{
    \item{\code{Integral}:}{It is an object of class \code{Integral.sde} and it is composed by the following slots:
    \describe{
        \item{\code{param.Integral}:}{it is an object of class \code{param.Integral} and it is composed by the following slots:
        \describe{
            \item{\code{allparam}:}{labels of all parameters (model and  integral).}
            \item{\code{common}:}{common parameters.}
            \item{\code{Integrandparam}:}{labels of all parameters only in the integral.}
          }
        }
        \item{\code{variable.Integral}:}{it is an object of class \code{variable.Integral} and it is composed by the following slots:
        \describe{
          \item{\code{var.dx}:}{integral variable.}
          \item{\code{lower.var}:}{lower bound of support.}
          \item{\code{upper.var}:}{upper bound of support.}
          \item{\code{out.var}:}{labels of output.}
          \item{\code{var.time}:}{label of time.}
          }
        }
        \item{\code{Integrand}:}{it is an object of class \code{variable.Integral} and it is composed by the following slots:
        \describe{
                \item{\code{IntegrandList}:}{It is a \code{list} that contains the components of integrand \code{h(theta, Xs, s)}.}
                \item{\code{dimIntegrand}:}{a \code{numeric} object that is the dimensions of the output.}
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
  }
}
