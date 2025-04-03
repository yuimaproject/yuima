\name{mllag}
\alias{mllag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Multiple Lead-Lag Detector
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Detecting the lead-lag parameters of discretely observed processes by picking time shifts at which the Hayashi-Yoshida cross-correlation functions exceed thresholds, which are constructed based on the asymptotic theory of Hayashi and Yoshida (2011).
}
\usage{
mllag(x, from = -Inf, to = Inf, division = FALSE, grid, psd = TRUE, 
      plot = TRUE, alpha = 0.01, fisher = TRUE, bw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}} or \code{yuima.llag-class} (output of \code{\link{llag}}) or \code{yuima.mllag-class} (output of this function).
}
  \item{from}{
%%     ~~Describe \code{from} here~~
passed to \code{\link{llag}}.
}
  \item{to}{
%%     ~~Describe \code{to} here~~
passed to \code{\link{llag}}.
}
  \item{division}{
%%     ~~Describe \code{division} here~~
passed to \code{\link{llag}}.
}
  \item{grid}{
%%     ~~Describe \code{grid} here~~
passed to \code{\link{llag}}.
}
  \item{psd}{
%%     ~~Describe \code{psd} here~~
passed to \code{\link{llag}}.
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
logical. If \code{TRUE}, the estimated cross-correlation functions and the pointwise confidence intervals (under the null hypothesis that the corresponding correlation is zero) as well as the detected lead-lag parameters are plotted. 
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
a posive number indicating the significance level of the confidence intervals for the cross-correlation functions.
}
  \item{fisher}{
%%     ~~Describe \code{fisher} here~~
logical. If \code{TRUE}, the p-values and the confidence intervals for the cross-correlation functions is evaluated after applying the Fisher z transformation.
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
passed to \code{\link{llag}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The computation method of cross-correlation functions and confidence intervals is the same as the one used in \code{\link{llag}}. The exception between this function and \code{\link{llag}} is how to detect the lead-lag parameters. While \code{\link{llag}} only returns the maximizer of the absolute value of the cross-correlations following the theory of Hoffmann et al. (2013),  this function returns all the time shifts at which the cross-correlations exceed (so there is also the possiblity that \emph{no} lead-lag is returned). Note that this approach is mathematically debetable because there would be a multiple testing problem (see also 'Note' of \code{\link{llag}}), so the interpretation of the result from this function should carefully be addressed. In particular, the significance level \code{alpha} probably does not give the "correct" level.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{"yuima.mllag"}, which is a list with the following elements:
\item{mlagcce}{a list of \code{data.frame}-class objects consisting of \code{lagcce} (lead-lag parameters), \code{p.value} and \code{correlation}.}
\item{LLR}{a matrix consisting of lead-lag ratios. See Huth and Abergel (2014) for details.}
\item{ccor}{a list of computed cross-correlation functions.}
\item{avar}{a list of computed asymptotic variances of the cross-correlations (if \code{ci = TRUE}).}
\item{CI}{a list of computed confidence intervals.}
}
\references{
%% ~put references to the literature/web site here ~
Hayashi, T. and Yoshida, N. (2011)
  Nonsynchronous covariation process and limit theorems,
  \emph{Stochastic processes and their applications}, \bold{121}, 2416--2454.
  
Hoffmann, M., Rosenbaum, M. and Yoshida, N. (2013)
  Estimation of the lead-lag parameter from non-synchronous data, 
  \emph{Bernoulli}, \bold{19}, no. 2, 426--461.
  
Huth, N. and Abergel, F. (2014)
  High frequency lead/lag relationships --- Empirical facts,
  \emph{Journal of Empirical Finance}, \bold{26}, 41--58.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA Project Team
}
%\note{
%%  ~~further notes~~
%If \code{x} is not of class \code{"yuima.llag"}, the function first passes its arguments to \code{\link{llag}}, then apply \code{mllag.yuima.llag} to the returned object.

%The object returned by this function is indeed of class \code{c("yuima.mllag", "yuima.llag")}, so \code{mllag} can again be applied to it.
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{llag}}, \code{\link{hyavar}}, \code{\link{llag.test}}
}
\examples{

# The first example is taken from llag

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

## lead-lag estimation by mllag
par(mfcol=c(3,1))
result <- mllag(psample, grid = G)

## Since the lead-lag parameter for the pair(x1, x3) is not contained in G,
## no lead-lag parameter is detected for this pair

par(mfcol=c(1,1))

# The second example is a situation where multiple lead-lag effects exist
set.seed(222)

n <- 3600
Times <- seq(0, 1, by = 1/n)
R1 <- 0.6
R2 <- -0.3

dW1 <- rnorm(n + 10)/sqrt(n)
dW2 <- rnorm(n + 5)/sqrt(n)
dW3 <- rnorm(n)/sqrt(n)

x <- zoo(diffinv(dW1[-(1:10)] + dW2[1:n]), Times)
y <- zoo(diffinv(R1 * dW1[1:n] + R2 * dW2[-(1:5)] + 
                 sqrt(1- R1^2 - R2^2) * dW3), Times)

## In this setting, both x and y have a component leading to the other, 
## but x's leading component dominates y's one

yuima <- setData(list(x, y))

## Lead-lag estimation by llag
G <- seq(-30/n, 30/n, by = 1/n)
est <- llag(yuima, grid = G, ci = TRUE)

## The shape of the plotted cross-correlation is evidently bimodal,
## so there are likely two lead-lag parameters

## Lead-lag estimation by mllag
mllag(est) # succeeds in detecting two lead-lag parameters

## Next consider a non-synchronous sampling case
psample <- poisson.random.sampling(yuima, n = n, rate = c(0.8, 0.7))

## Lead-lag estimation by mllag
est <- mllag(psample, grid = G) 
est # detects too many lead-lag parameters

## Using a lower significant level
mllag(est, alpha = 0.001) # insufficient

## As the plot reveals, one reason is because the grid is too dense
## In fact, this phenomenon can be avoided by using a coarser grid
mllag(psample, grid = seq(-30/n, 30/n, by=5/n)) # succeeds!

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
