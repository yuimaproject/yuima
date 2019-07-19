\name{bns.test}
\alias{bns.test}
\alias{bns.test,yuima-method}
\alias{bns.test,yuima.data-method}
\alias{bns.test,list-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Barndorff-Nielsen and Shephard's Test for the Presence of Jumps Using Bipower Variation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Tests the presence of jumps using the statistic proposed in Barndorff-Nielsen and Shephard (2004,2006) for each component.
}
\usage{
bns.test(yuima, r = rep(1, 4), type = "standard", adj = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{
%%     ~~Describe \code{yuima} here~~
an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.
}
  \item{r}{
%%     ~~Describe \code{r} here~~
a vector of non-negative numbers or a list of vectors of non-negative numbers. Theoretically, it is necessary that \code{sum(r)=4} and \code{max(r)<2}.
}
  \item{type}{
%%     ~~Describe \code{type} here~~
type of the test statistic to use. \code{standard} is default.
}
  \item{adj}{
%%     ~~Describe \code{adj} here~~
logical; if \code{TRUE}, the maximum adjustment suggested in Barndorff-Nielsen and Shephard (2004) is applied to the test statistic when \code{type} is equal to either \dQuote{\code{log}} or \dQuote{\code{ratio}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For the \code{i}-th component, the test statistic is equal to the \code{i}-th component of \code{sqrt(n)*(mpv(yuima,2)-mpv(yuima,c(1,1)))/sqrt(vartheta*mpv(yuima,r))} when \code{type="standard"}, \code{sqrt(n)*log(mpv(yuima,2)/mpv(yuima,c(1,1)))/sqrt(vartheta*mpv(yuima,r)/mpv(yuima,c(1,1))^2)} when \code{type="log"} and \code{sqrt(n)*(1-mpv(yuima,c(1,1))/mpv(yuima,2))/sqrt(vartheta*mpv(yuima,r)/mpv(yuima,c(1,1))^2)} when \code{type="ratio"}. Here, \code{n} is equal to the length of the \code{i}-th component of the \code{zoo.data} of \code{yuima} minus 1 and \code{vartheta} is \code{pi^2/4+pi-5}. When \code{adj=TRUE}, \code{mpv(yuima,r)[i]/mpv(yuima,c(1,1))^2)[i]} is replaced with 1 if it is less than 1.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list with the same length as the \code{zoo.data} of \code{yuima}. Each components of the list has class \dQuote{\code{htest}} and contains the following components:
\item{statistic}{the value of the test statistic of the corresponding component of the \code{zoo.data} of \code{yuima}.}
\item{p.value}{an approximate p-value for the test of the corresponding component.}
\item{method}{the character string \dQuote{\code{Barndorff-Nielsen and Shephard jump test}}.}
\item{data.name}{the character string \dQuote{\code{xi}}, where \code{i} is the number of the component.}
}
\references{
%% ~put references to the literature/web site here ~
Barndorff-Nielsen, O. E. and Shephard, N. (2004)
  Power and bipower variation with stochastic volatility and jumps,
  \emph{Journal of Financial Econometrics}, \bold{2}, no. 1, 1--37.
  
Barndorff-Nielsen, O. E. and Shephard, N. (2006)
  Econometrics of testing for jumps in financial economics using bipower variation,
  \emph{Journal of Financial Econometrics}, \bold{4}, no. 1, 1--30.
  
Huang, X. and Tauchen, G. (2005)
  The relative contribution of jumps to total price variance,
  \emph{Journal of Financial Econometrics}, \bold{3}, no. 4, 456--499.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA Project Team
}
\note{
%%  ~~further notes~~
Theoretically, this test may be invalid if sampling is irregular.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{mpv}}
}
\examples{

set.seed(123)

# One-dimensional case
## Model: dXt=t*dWt+t*dzt, 
## where zt is a compound Poisson process with intensity 5 and jump sizes distribution N(0,0.1).

model <- setModel(drift=0,diffusion="t",jump.coeff="t",measure.type="CP",
                  measure=list(intensity=5,df=list("dnorm(z,0,sqrt(0.1))")),
                  time.variable="t")

yuima.samp <- setSampling(Terminal = 1, n = 390) 
yuima <- setYuima(model = model, sampling = yuima.samp) 
yuima <- simulate(yuima)
plot(yuima) # The path seems to involve some jumps

bns.test(yuima) # standard type

bns.test(yuima,type="log") # log type

bns.test(yuima,type="ratio") # ratio type

# Multi-dimensional case
## Model: dXkt=t*dWk_t (k=1,2,3) (no jump case).

diff.matrix <- diag(3)
diag(diff.matrix) <- c("t","t","t")
model <- setModel(drift=c(0,0,0),diffusion=diff.matrix,time.variable="t",
                  solve.variable=c("x1","x2","x3"))

yuima.samp <- setSampling(Terminal = 1, n = 390) 
yuima <- setYuima(model = model, sampling = yuima.samp) 
yuima <- simulate(yuima)
plot(yuima)

bns.test(yuima)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
% __ONLY ONE__ keyword per line
