\name{ntv}
\alias{ntv}
\alias{minrv}
\alias{medrv}
\alias{minrv.test}
\alias{medrv.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Volatility Estimation and Jump Test Using Nearest Neighbor Truncation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{minrv} and \code{medrv} respectively compute the MinRV and MedRV estimators introduced in Andersen, Dobrev and Schaumburg (2012).

\code{minrv.test} and \code{medrv.test} respectively perform Haussman type tests for the null hypothesis that the realized path has no jump using the MinRV and MedRV estimators. 
See Section 4.4 in Andersen, Dobrev and Schaumburg (2014) for a concise discussion. 
}
\usage{
minrv(yuima)
medrv(yuima)

minrv.test(yuima, type = "ratio", adj = TRUE)
medrv.test(yuima, type = "ratio", adj = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{
%%     ~~Describe \code{yuima} here~~
an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.
}
\item{type}{
%%     ~~Describe \code{type} here~~
type of the test statistic to use. \code{ratio} is default.
}
  \item{adj}{
%%     ~~Describe \code{adj} here~~
logical; if \code{TRUE}, the maximum adjustment suggested in Barndorff-Nielsen and Shephard (2004) is applied to the test statistic when \code{type} is equal to either \dQuote{\code{log}} or \dQuote{\code{ratio}}. See also Section 2.5 in Dumitru and Urga (2012).
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{minrv} and \code{medrv} return a numeric vector with the same length as \code{dim(yuima)}. Each component of the vector is a volatility estimate for the corresponding component of \code{yuima}.

\code{minrv.test} and \code{medrv.test} return a list with the same length as \code{dim(yuima)}. Each component of the list has class \dQuote{\code{htest}} and contains the following components:
\item{statistic}{the value of the test statistic of the corresponding component of \code{yuima}.}
\item{p.value}{an approximate p-value for the test of the corresponding component.}
\item{method}{the character string \dQuote{\code{Andersen-Dobrev-Schaumburg jump test based on xxx}}, where xxx is either MinRV or MedRV.}
\item{data.name}{the character string \dQuote{\code{xi}}, where \code{i} is the number of the component.}
}
\references{
%% ~put references to the literature/web site here ~
Andersen, T. G., Dobrev D. and Schaumburg, E. (2012)
  Jump-robust volatility estimation using nearest neighbor truncation.
  \emph{Journal of Econometrics}, \bold{169}, 75--93.
  
Andersen, T. G., Dobrev D. and Schaumburg, E. (2014)
  A robust neighborhood truncation approach to estimation of integrated quarticity.
  \emph{Econometric Theory}, \bold{30}, 3--59.
  
Dumitru, A.-M. and Urga, G. (2012)
  Identifying jumps in financial assets: A comparison between nonparametric jump tests.
  \emph{Journal of Business and Economic Statistics}, \bold{30}, 242--255.
  
Maneesoonthorn, W., Martin, G. M. and Forbes, C. S. (2020)
  High-frequency jump tests: Which test should we use?
  \emph{Journal of Econometrics}, \bold{219}, 478--487.
  
Theodosiou, M. and Zikes, F. (2011)
  A comprehensive comparison of alternative tests for jumps in asset prices. 
  Central Bank of Cyprus Working Paper 2011-2.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA Project Team
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{mpv}}, \code{\link{cce}}, \code{\link{bns.test}}, \code{\link{lm.jumptest}}, \code{\link{pz.test}}
}
\examples{
set.seed(123)

# One-dimensional case
## Model: dXt=t*dWt+t*dzt, 
## where zt is a compound Poisson process with intensity 5
## and jump sizes distribution N(0,1).

model <- setModel(drift=0,diffusion="t",jump.coeff="t",measure.type="CP",
                  measure=list(intensity=5,df=list("dnorm(z,0,1)")),
                  time.variable="t")

yuima.samp <- setSampling(Terminal = 1, n = 390) 
yuima <- setYuima(model = model, sampling = yuima.samp) 
yuima <- simulate(yuima)
plot(yuima) # The path evidently has some jumps

## Volatility estimation
minrv(yuima) # minRV (true value = 1/3)
medrv(yuima) # medRV (true value = 1/3)

## Jump test
minrv.test(yuima, type = "standard")
minrv.test(yuima,type="log")
minrv.test(yuima,type="ratio")

medrv.test(yuima, type = "standard")
medrv.test(yuima,type="log")
medrv.test(yuima,type="ratio")


# Multi-dimensional case
## Model: Bivariate standard BM + CP
## Only the first component has jumps

mod <- setModel(drift = c(0, 0), diffusion = diag(2),
                jump.coeff = diag(c(1, 0)),
                measure = list(intensity = 5, 
                               df = "dmvnorm(z,c(0,0),diag(2))"),
                jump.variable = c("z"), measure.type=c("CP"),
                solve.variable=c("x1","x2"))

samp <- setSampling(Terminal = 1, n = 390) 
yuima <- setYuima(model = model, sampling = yuima.samp)
yuima <- simulate(object = mod, sampling = samp)
plot(yuima)

## Volatility estimation
minrv(yuima) # minRV (true value = c(1, 1))
medrv(yuima) # medRV (true value = c(1, 1))

## Jump test
minrv.test(yuima) # test is performed component-wise
medrv.test(yuima) # test is performed component-wise
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ts}
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
