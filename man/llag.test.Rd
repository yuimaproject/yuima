\name{llag.test}
\alias{llag.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Wild Bootstrap Test for the Absence of Lead-Lag Effects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Tests the absence of lead-lag effects (time-lagged correlations) by the wild bootstrap procedure proposed in Koike (2017) for each pair of components.
}
\usage{
llag.test(x, from = -Inf, to = Inf, division = FALSE, grid, R = 999,
          parallel = "no", ncpus = getOption("boot.ncpus", 1L),
          cl = NULL, tol = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.
}
  \item{from}{
%%     ~~Describe \code{from} here~~
a numeric vector each of whose component(s) indicates the lower end of a finite grid on which the contrast function is evaluated, if \code{grid} is missing.
}
  \item{to}{
%%     ~~Describe \code{to} here~~
a numeric vector each of whose component(s) indicates the upper end of a finite grid on which the contrast function is evaluated, if \code{grid} is missing.
}
  \item{division}{
%%     ~~Describe \code{division} here~~
a numeric vector each of whose component(s) indicates the number of the points of a finite grid on which the contrast function is evaluated, if \code{grid} is missing.
}
  \item{grid}{
%%     ~~Describe \code{grid} here~~
a numeric vector or a list of numeric vectors. See 'Details' of \code{\link{llag}}.
}
  \item{R}{
%%     ~~Describe \code{R} here~~
a single positive integer indicating the number of bootstrap replicates.
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
passed to \code{\link{boot}}.
}
  \item{ncpus}{
%%     ~~Describe \code{ncpus} here~~
passed to \code{\link{boot}}.
}
  \item{cl}{
%%     ~~Describe \code{cl} here~~
passed to \code{\link{boot}}.
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
tolelance parameter to avoid numerical errors in comparison of time stamps. All time stamps are divided by \code{tol} and rounded to integers. Note that the values of \code{grid} are also divided by \code{tol} and rounded to integers. A reasonable choice of \code{tol} is the minimum unit of time stamps. The default value \code{1e-6} supposes that the minimum unit of time stamps is greater or equal to 1 micro-second.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For each pair of components, this function performs the wild bootstrap procedure proposed in Koike (2017) to test whether there is a (possibly) time-lagged correlation. The null hypothesis of the test is that there is no time-lagged correlation and the alternative is its negative. The test regects the null hypothesis if the maximum of the absolute values of cross-covariances is too large. The critical region is constructed by a wild bootstrap procedure with Rademacher variables as the multiplier variables.  
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{p.values}{a matrix whose components indicate the bootstrap p-values for the corresponding pair of components.}
\item{max.cov}{a matrix whose componenets indicate the maxima of the absolute values of cross-covariances for the corresponding pair of components.}
\item{max.corr}{a matrix whose componenets indicate the maxima of the absolute values of cross-correlations for the corresponding pair of components.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Koike, Y. (2019). 
Gaussian approximation of maxima of Wiener functionals and its application to high-frequency data, 
\emph{Annals of Statistics}, \bold{47}, 1663--1687.
\href{https://arxiv.org/abs/1709.00353}{arXiv:1709.00353}.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA Project Team
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{cce}}, \code{\link{hyavar}}, \code{\link{mllag}}, \code{\link{llag}}
}
\examples{
# The following example is taken from mllag

## Set a model
diff.coef.matrix <- matrix(c("sqrt(x1)", "3/5*sqrt(x2)",
 "1/3*sqrt(x3)", "", "4/5*sqrt(x2)","2/3*sqrt(x3)",
 "","","2/3*sqrt(x3)"), 3, 3) 
drift <- c("1-x1","2*(10-x2)","3*(4-x3)")
cor.mod <- setModel(drift = drift, 
 diffusion = diff.coef.matrix,
 solve.variable = c("x1", "x2","x3")) 

set.seed(111) 

## We use a function poisson.random.sampling 
## to get observation by Poisson sampling.
yuima.samp <- setSampling(Terminal = 1, n = 1200) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima,xinit=c(1,7,5)) 

## intentionally displace the second time series

  data2 <- yuima@data@zoo.data[[2]]
  time2 <- time(data2)
  theta2 <- 0.05   # the lag of x2 behind x1
  stime2 <- time2 + theta2  
  time(yuima@data@zoo.data[[2]]) <- stime2

  data3 <- yuima@data@zoo.data[[3]]
  time3 <- time(data3)
  theta3 <- 0.12   # the lag of x3 behind x1
  stime3 <- time3 + theta3 
  time(yuima@data@zoo.data[[3]]) <- stime3

## sampled data by Poisson rules
psample<- poisson.random.sampling(yuima, 
 rate = c(0.2,0.3,0.4), n = 1000) 
 
 ## We search lead-lag parameters on the interval [-0.1, 0.1] with step size 0.01 
G <- seq(-0.1,0.1,by=0.01)

## perform lead-lag test
llag.test(psample, grid = G, R = 999)

## Since the lead-lag parameter for the pair(x1, x3) is not contained in G,
## the null hypothesis is not rejected for this pair
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}% use one of  RShowDoc("KEYWORDS")
