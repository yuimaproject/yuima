\name{yuima.cogarch-class}
\docType{class}
\alias{yuima.cogarch-class}
\alias{initialize,yuima.cogarch-method}
\alias{initialize,cogarch.info-method}
\alias{limiting.gamma,yuima.cogarch-method}
\alias{simulate,yuima.cogarch-method}

\title{Class for the mathematical description of CoGarch(p,q) model}
\description{
The \code{yuima.cogarch} class is a class of the  \pkg{yuima} package that extends the \code{\link{yuima.model-class}}. 
}
\section{Objects from the Class}{
Objects can be created by calls of the function \code{\link{setCogarch}}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{info}:}{is an \code{\link{cogarch.info-class}} object that describes the structure of the Cogarch(p,q) model.}
    \item{\code{drift}:}{is an \R{} expression which specifies the drift 
   coefficient (a vector).}
    \item{\code{diffusion}:}{is an \R{} expression which specifies the diffusion 
     coefficient (a matrix).}
    \item{\code{hurst}:}{the Hurst parameter of the gaussian noise. }
    \item{\code{jump.coeff}:}{a vector of \code{"expressions"} for the jump component. }
    \item{\code{measure}:}{Levy measure for the jump component. }
    \item{\code{measure.type}:}{Type of specification for Levy measure }
    \item{\code{parameter}:}{is an object of class \code{\link{model.parameter-class}}. }
    \item{\code{state.variable}:}{the state variable. }
    \item{\code{jump.variable}:}{the jump variable. }
    \item{\code{time.variable}:}{the time variable. }
    \item{\code{noise.number}:}{Object of class \code{"numeric"} }
    \item{\code{equation.number}:}{dimension of the stochastic differential equation.}
    \item{\code{dimension}:}{number of parameters. }
    \item{\code{solve.variable}:}{the solve variable}
    \item{\code{xinit}:}{Object of class \code{"expression"} that contains the starting function for the SDE. }
    \item{\code{J.flag}:}{wheather jump.coeff include jump.variable.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{yuima.model}"}, directly.
}
\section{Methods}{
  \describe{
    \item{simulate}{simulation method. For more information see \code{\link{simulate}} }
    \item{toLatex}{This method converts an object of \code{yuima.cogarch-class} to character vectors with LaTeX markup.
	 }
   \item{qmle}{Quasi maximum likelihood estimation procedure. For more information see \code{\link{qmle}}. }
}
}
\author{The YUIMA Project Team}
\keyword{classes}
