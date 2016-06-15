\name{cogarchNoise}
\alias{Recovering.Noise.cogarch}
\alias{cogarch.Recovering}
\alias{CogarchRecovNoise}
\alias{cogarchNoise}
\alias{Levy.cogarch}
\title{Estimation for the underlying Levy in a COGARCH(p,q) model}
\description{Retrieve the increment of the underlying Levy for the COGARCH(p,q) process}
\usage{
cogarchNoise(yuima, data=NULL, param, mu=1)
}
\arguments{
  \item{yuima}{a yuima object or an object of \code{\link{yuima.cogarch-class}}.}
  \item{data}{an object of class \code{\link{yuima.data-class}} contains the observations available at uniformly spaced time. If \code{data=NULL}, the default, the \code{cogarchNoise} uses the data in an object of \code{\link{yuima.data-class}}.}
  \item{param}{\code{list} of parameters for the COGARCH(p,q).}
  \item{mu}{a numeric object that contains the value of the second moments of the levy measure.}
}
\value{
  %\item{QL}{a real value.}
  \item{incr.Levy}{a numeric object contains the estimated increments.}
  \item{model}{an object of class yuima containing the state, the variance and the cogarch process.}
  }
\author{The YUIMA Project Team}
\note{
  %The function ml.ql uses the function optim internally.
  The function \code{cogarchNoise} assumes the underlying Levy process is centered in zero.

  The function \code{gmm} uses the function \code{cogarchNoise} for estimation of underlying Levy in the COGARCH(p,q) model.
  }

 \references{
               Chadraa. (2009)
Statistical Modelling with COGARCH(P,Q) Processes, \emph{PhD Thesis}.
}


\examples{
# Insert here some examples
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
