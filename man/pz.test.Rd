\name{pz.test}
\alias{pz.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Podolskij and Ziggel's Test for the Presence of Jumps Using Power Variation with Perturbed Truncation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Performs a test for the null hypothesis that the realized path has no jump following Podolskij and Ziggel (2010).
}
\usage{
pz.test(yuima, p = 4, threshold = "local", tau = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{
%%     ~~Describe \code{yuima} here~~
an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.
}
  \item{p}{
%%     ~~Describe \code{p} here~~
a positive number indicating the exponent of the (truncated) power variation to compute test statistic(s). Theoretically, it must be greater than or equal to 2.
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
a numeric vector or list indicating the threshold parameter(s). Each of its components indicates the threshold parameter or process to be used for estimating the corresponding component. If it is a numeric vector, the elements in \code{threshold} are recycled if there are too few elements in \code{threshold}. 

Alternatively, you can specify either \code{"PZ"} or \code{"local"} to automatically select a (hopefully) appropriate threshold. When \code{threshold="PZ"}, selection is performed following Section 5.1 in Podolskij and Ziggel (2010). When \code{threshold="local"}, selection is performed following Section 5.1 in Koike (2014). The default is \code{threshold="local"}.
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
a probability controlling the strength of perturbation. See Section 2.3 in Podolskij and Ziggel (2010) for details. Podolskij and Ziggel (2010) suggests using a relatively small value for \code{tau}, e.g. \code{tau=0.1} or \code{tau=0.05}.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list with the same length as \code{dim(yuima)}. Each component of the list has class \dQuote{\code{htest}} and contains the following components:
\item{statistic}{the value of the test statistic of the corresponding component of \code{yuima}.}
\item{p.value}{an approximate p-value for the test of the corresponding component.}
\item{method}{the character string \dQuote{\code{Podolskij and Ziggel jump test}}.}
\item{data.name}{the character string \dQuote{\code{xi}}, where \code{i} is the number of the component.}
}
\references{
%% ~put references to the literature/web site here ~
Dumitru, A.-M. and Urga, G. (2012)
  Identifying jumps in financial assets: A comparison between nonparametric jump tests.
  \emph{Journal of Business and Economic Statistics}, \bold{30}, 242--255.
  
Koike, Y. (2014)
  An estimator for the cumulative co-volatility of asynchronously observed semimartingales with jumps,
  \emph{Scandinavian Journal of Statistics}, \bold{41}, 460--481.
  
Maneesoonthorn, W., Martin, G. M. and Forbes, C. S. (2020)
  High-frequency jump tests: Which test should we use?
  \emph{Journal of Econometrics}, \bold{219}, 478--487.

Podolskij, M. and Ziggel, D. (2010)
  New tests for jumps in semimartingale models,
  \emph{Statistical Inference for Stochastic Processes}, \bold{13}, 15--41.
  
Theodosiou, M. and Zikes, F. (2011)
  A comprehensive comparison of alternative tests for jumps in asset prices. 
  Central Bank of Cyprus Working Paper 2011-2.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA Project Team
}
\note{
%%  ~~further notes~~
Podolskij and Ziggel (2010) also introduce a pre-averaged version of the test to deal with noisy observations. Such a test will be implemented in the future version of the package.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{bns.test}}, \code{\link{lm.jumptest}}, \code{\link{minrv.test}}, \code{\link{medrv.test}}
}
\examples{
\dontrun{
set.seed(123)

# One-dimensional case
## Model: dXt=t*dWt+t*dzt, 
## where zt is a compound Poisson process with intensity 5 and jump sizes distribution N(0,1).

model <- setModel(drift=0,diffusion="t",jump.coeff="t",measure.type="CP",
                  measure=list(intensity=5,df=list("dnorm(z,0,sqrt(0.1))")),
                  time.variable="t")

yuima.samp <- setSampling(Terminal = 1, n = 390) 
yuima <- setYuima(model = model, sampling = yuima.samp) 
yuima <- simulate(yuima)
plot(yuima) # The path seems to involve some jumps

#lm.jumptest(yuima) # p-value is very small, so the path would have a jump
#lm.jumptest(yuima, K = floor(sqrt(390))) # different value of K
pz.test(yuima) # p-value is very small, so the path would have a jump
pz.test(yuima, p = 2) # different value of p
pz.test(yuima, tau = 0.1) # different value of tau

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

pz.test(yuima) # test is performed component-wise
}
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
