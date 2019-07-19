\name{mpv}
\alias{mpv}
%\alias{mpv,list-method}
\alias{mpv,yuima-method}
\alias{mpv,yuima.data-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Realized Multipower Variation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function returns the realized MultiPower Variation (mpv), defined in Barndorff-Nielsen and Shephard (2004), for each components.
}
\usage{
mpv(yuima, r = 2, normalize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{
%%     ~~Describe \code{yuima} here~~
an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.
}
  \item{r}{
%%     ~~Describe \code{r} here~~
a vector of non-negative numbers or a list of vectors of non-negative numbers. 
}
\item{normalize}{logical. See `Details'.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Let \code{d} be the number of the components of the \code{zoo.data} of \code{yuima}.

Let \eqn{X^i_{t_0},X^i_{t_1},\dots,X^i_{t_n}} be the observation data of the \eqn{i}-th component (i.e. the \eqn{i}-th component of the \code{zoo.data} of \code{yuima}).

When \eqn{r} is a \eqn{k}-dimensional vector of non-negative numbers, \code{mpv(yuima,r,normalize=TRUE)} is defined as the \code{d}-dimensional vector with \code{i}-th element equal to
\deqn{\mu_{r[1]}^{-1}\cdots\mu_{r[k]}^{-1}n^{\frac{r[1]+\cdots+r[k]}{2}-1}\sum_{j=1}^{n-k+1}|\Delta X^i_{t_{j}}|^{r[1]}|\Delta X^i_{t_{j+1}}|^{r[2]}\cdots|\Delta X^i_{t_{j+k-1}}|^{r[k]},}
where \eqn{\mu_p} is the p-th absolute moment of the standard normal distribution and \eqn{\Delta X^i_{t_{j}}=X^i_{t_j}-X^i_{t_{j-1}}}. If \code{normalize} is \code{FALSE} the result is not multiplied by \eqn{\mu_{r[1]}^{-1}\cdots\mu_{r[k]}^{-1}}.

When \eqn{r} is a list of vectors of non-negative numbers, \code{mpv(yuima,r,normalize=TRUE)} is defined as the \code{d}-dimensional vector with \code{i}-th element equal to
\deqn{\mu_{r^i_1}^{-1}\cdots\mu_{r^i_{k_i}}^{-1}n^{\frac{r^i_1+\cdots+r^i_{k_i}}{2}-1}\sum_{j=1}^{n-k_i+1}|\Delta X^i_{t_{j}}|^{r^i_1}|\Delta X^i_{t_{j+1}}|^{r^i_2}\cdots|\Delta X^i_{t_{j+k_i-1}}|^{r^i_{k_i}},}
where \eqn{r^i_1,\dots,r^i_{k_i}} is the \code{i}-th component of \code{r}. If \code{normalize} is \code{FALSE} the result is not multiplied by \eqn{\mu_{r^i_1}^{-1}\cdots\mu_{r^i_{k_i}}^{-1}}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A numeric vector with the same length as the \code{zoo.data} of \code{yuima}
}
\references{
%% ~put references to the literature/web site here ~
Barndorff-Nielsen, O. E. and Shephard, N. (2004)
  Power and bipower variation with stochastic volatility and jumps,
  \emph{Journal of Financial Econometrics}, \bold{2}, no. 1, 1--37.
  
Barndorff-Nielsen, O. E. , Graversen, S. E. , Jacod, J. , Podolskij M. and Shephard, N. (2006)
  A central limit theorem for realised power and bipower variations of continuous semimartingales,
  in: Kabanov, Y. , Lipster, R. , Stoyanov J. (Eds.), From Stochastic Calculus to Mathematical Finance: The Shiryaev Festschrift, Springer-Verlag, Berlin, pp. 33--68.
}
\author{
%%  ~~who you are~~
Yuta Koike with YUIMA Project Team
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{setModel}},\code{\link{cce}}
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
plot(yuima)

mpv(yuima) # true value is 1/3
mpv(yuima,1) # true value is 1/2
mpv(yuima,rep(2/3,3)) # true value is 1/3

# Multi-dimensional case
## Model: dXkt=t*dWk_t (k=1,2,3).

diff.matrix <- diag(3)
diag(diff.matrix) <- c("t","t","t")
model <- setModel(drift=c(0,0,0),diffusion=diff.matrix,time.variable="t",
                  solve.variable=c("x1","x2","x3"))

yuima.samp <- setSampling(Terminal = 1, n = 390) 
yuima <- setYuima(model = model, sampling = yuima.samp) 
yuima <- simulate(yuima)
plot(yuima)

mpv(yuima,list(c(1,1),1,rep(2/3,3))) # true varue is c(1/3,1/2,1/3)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
