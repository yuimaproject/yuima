\name{setLaw_th}
\alias{setLaw_th}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructior of a t-Levy process.
}
\description{
\code{setLaw_th} constructs an object of class \code{\link{yuima.th-class}}.
}
\usage{
setLaw_th(h = 1, method = "LAG", up = 7, low = -7, N = 180, 
N_grid = 1000, regular_par = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
a numeric object that is the time of the intervals
}
  \item{method}{
Method for the inversion of the characteristic function. Three methods are available: \code{cos}, \code{LAG}, and \code{FFT}. 
}
  \item{up}{
Upper bound for the integration support.
}
  \item{low}{
Lower bound for the integration support.
}
  \item{N}{
Integration grid.
}
  \item{N_grid}{
Number of points in the support.
}
  \item{regular_par}{
A scalar for controlling the Gibbs effect for the inversion of the characteristic function
}
  \item{\dots}{
 Additional arguments. See \code{\link{setLaw}} for more details.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The function returns an object of class \code{\link{yuima.th-class}}.}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
The YUIMA Project Team

Contacts: Lorenzo Mercuri \email{lorenzo.mercuri@unimi.it}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
