\name{FromCF2yuima_law}
\alias{FromCF2yuima_law}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
From a Characteristic Function to an  \code{yuima.law-object}
}
\description{This function returns an object of \code{\link{yuima.law-class}} and requires the characteristic function as the only input. Density, Random Number Generator, Cumulative Distribution Function and quantile function are internally constructed}
\usage{
FromCF2yuima_law(myfun, time.names = "t", var_char = "u", up = 45, 
low = -45, N_grid = 50001, N_Fourier = 2^10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myfun}{ A \code{string} that is the name of the characteristic function defined by Users.}
  \item{time.names}{ Label of time.}
  \item{var_char}{ Argument of the characteristic function.}
  \item{up}{ Upper bound for the internal integration.}
  \item{low}{ Lower bound for the internal integration.}
  \item{N_grid}{ Observation grid. }
  \item{N_Fourier}{ Number of points for the Fourier Inversion.}
}
\details{ The density function is obtained by means of the Fourier Transform.}
\value{
 An object of \code{\link{yuima.law-class}}.
}
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
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
