\name{qmle.degenerate}
\alias{qmle.degenerate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quasi-likelihood estimation for degenerate diffusion processes
}
\description{
Calculate a quasi maximum likelihood type estimator for the parameters of a degenerate diffusion process based on discrete observation data.
}
\usage{
qmle.degenerate(yuima, start, method = "L-BFGS-B", 
                fixed = list(), print = FALSE, 
                envir = globalenv(), lower, upper, 
                joint = FALSE, ...)
}
\arguments{
  \item{yuima}{a yuima object.}
  \item{start}{initial values to be passed to the optimizer.}
  \item{method}{passed to \code{\link{optim}}.}
  \item{fixed}{a list of fixed parameters.}
  \item{print}{you can see a progress of the estimation when print is \code{TRUE}.}
  \item{envir}{an environment where the model coefficients are evaluated.}
  \item{lower}{a named list for specifying lower bounds of parameters.}
  \item{upper}{a named list for specifying upper bounds of parameters.}
  \item{joint}{perform joint estimation or adaptive estimation with one-step correction? by default \code{joint=FALSE}.}
  \item{\dots}{passed to \code{\link{optim}}.}
}
\details{
When \code{joint=TRUE}, the function calculates the non-adaptive estimator presented in Section 6 of Gloter and Yoshida (2020). 

When \code{joint=FALSE}, the function calculates the adaptive one-step estimator presented in Section 4 of Gloter and Yoshida (2021) with the initial estimator constructed as in Section 3 of Gloter and Yoshida (2021).
}
\value{
An object of class \code{\link{yuima.qmle-class}}.

When \code{joint=FALSE}, the slots \code{min}, \code{details} and \code{minuslogl} are for estimation of the drift parameter of the degenerate components, but the \code{details} slot has an additional component \code{initial} containing the initial estimator. 
}
\references{
Gloter, A. and Yoshida, N. (2020)
  Adaptive and non-adaptive estimation for degenerate diffusion processes, working paper, available at \doi{10.48550/arXiv.2002.10164}.
  
Gloter, A. and Yoshida, N. (2021)
  Adaptive estimation for degenerate diffusion processes, \emph{Electronic Journal of Statistics}, \bold{15}, 1424--1472.
}
\author{
Yuta Koike with YUIMA Project Team
}
\note{
The model coefficients must be symbolically differentiable. 

If "x" is a component of \code{state.variable}, the character \code{"deltax"} should not be used as any variable name of the model as they are internally used for other variable names. 

Time-inhomogeneous diffusions are not supported.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qmle}}
}
\examples{
To be written.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ts}
