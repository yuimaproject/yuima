\name{estimation_LRM}
\alias{estimation_LRM}
\alias{Estimation of t-Levy Regression Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of the t-Levy Regression Model}
\description{The function estimates a t-Levy Regression Model}
\usage{
estimation_LRM(start, model, data, upper, lower, PT = 500, n_obs1 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
Initial values to be passed to the optimizer.
}
  \item{model}{
  A \code{\link{yuima.LevyRM-class}} that contains the mathematical representation of the t-Levy Regression Model. Its slot \code{@data} can contain either real or simulated data.
  }
  \item{data}{	
An object of class \code{\link{yuima.data-class}} contains the observations available at uniformly spaced time. If \code{data=NULL}, the default, the function uses the data in the object \code{model}.
}
  \item{upper}{A named list for specifying upper bounds of parameters.}
  \item{lower}{A named list for specifying lower bounds of parameters.}
  \item{PT}{
 The number of the data for the estimation of the regressor coefficients and the scale parameter. 
}
  \item{n_obs1}{
The number of data used in the estimation of the degree of freedom. As default the number of the whole data is used in this part
}
}
\details{A two-step estimation procedure. Regressor coefficients and scale parameters are obtained by maximizing the quasi-likelihood function based on the Cauchy density. The degree of freedom is estimated used the unitary increment of the t-noise.}
\value{Estimated parameters}
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
