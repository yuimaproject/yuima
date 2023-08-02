\name{qmle.degenerate}
\alias{qmle.degenerate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Quasi-likelihood estimation for degenerate diffusion processes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qmle.degenerate(yuima, start, method = "L-BFGS-B", 
                envir = globalenv(), lower, upper, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{a yuima object.}
  \item{start}{initial values to be passed to the optimizer.}
  \item{method}{passed to \code{\link{optim}}.}
  \item{envir}{an environment where the model coefficients are evaluated.}
  \item{lower}{a named list for specifying lower bounds of parameters.}
  \item{upper}{a named list for specifying upper bounds of parameters.}
  \item{\dots}{passed to \code{\link{optim}}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%%  ~~further notes~~
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
