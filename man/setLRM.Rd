\name{setLRM}
\alias{setLRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A constructor of a t-Student Regression Model.
}
\description{
This function returns an object of \code{\link{yuima.LevyRM-class}}
}
\usage{
setLRM(unit_Levy, yuima_regressors, LevyRM = "Y", coeff = c("mu", "sigma0"), 
data = NULL, sampling = NULL, characteristic = NULL, functional = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{unit_Levy}{An object of \code{\link{yuima.th-class}} that describes the t - noise in the regression model.
}
  \item{yuima_regressors}{An object of \code{\link{yuima.model-class}} that represents the regressors.}
  \item{LevyRM}{The label of the output variable. Default \code{'Y'}.}
  \item{coeff}{Labels for the regressor coefficients and the scale parameter.}
  \item{data}{An object of \code{\link{yuima.data-class}} that contains simulated or real data.}
  \item{sampling}{An object of \code{\link{yuima.sampling-class}}.}
  \item{characteristic}{An object of \code{\link{yuima.characteristic-class}}.}
  \item{functional}{An object of class \code{\link{yuima.functional-class}}.}
  \item{\dots}{Additional arguments. See \code{\link{setYuima}}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An object of \code{\link{yuima.LevyRM-class}.}
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
