\name{qmleDegenerate}
\alias{qmleDegenerate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quasi-likelihood Estimation for Degenerate Diffusion Processes
}
\description{
Calculate a quasi maximum likelihood type estimator for the parameters of a degenerate diffusion process based on discrete observation data.
}
\usage{
qmleDegenerate(yuima, start, method = "L-BFGS-B", 
                fixed = list(), print = FALSE, 
                envir = globalenv(), lower, upper, 
                joint = FALSE, ...)
}
\arguments{
  \item{yuima}{a yuima object.}
  \item{start}{initial values to be passed to the optimizer.}
  \item{method}{passed to \code{\link{optim}}.}
  \item{fixed}{a list of fixed parameters.}
  \item{print}{you can see a progress of the estimation when print is \code{TRUE}.}
  \item{envir}{an environment where the model coefficients are evaluated.}
  \item{lower}{a named list for specifying lower bounds of parameters.}
  \item{upper}{a named list for specifying upper bounds of parameters.}
  \item{joint}{perform joint estimation or adaptive estimation with one-step correction? by default \code{joint=FALSE}.}
  \item{\dots}{passed to \code{\link{optim}}.}
}
\details{
When \code{joint=TRUE}, the function calculates the non-adaptive estimator presented in Gloter and Yoshida (2024). 

When \code{joint=FALSE}, the function calculates the adaptive one-step estimator presented in Section 4 of Gloter and Yoshida (2021) with the initial estimator constructed as in Section 3 of Gloter and Yoshida (2021).
}
\value{
An object of class \code{\link{yuima.qmle-class}}.

When \code{joint=FALSE}, the slots \code{min}, \code{details} and \code{minuslogl} are for estimation of the drift parameter of the degenerate components, but the \code{details} slot has an additional component \code{initial} containing the initial estimator. 
}
\references{
Gloter, A. and Yoshida, N. (2020)
  Adaptive and non-adaptive estimation for degenerate diffusion processes, working paper, available at \doi{10.48550/arXiv.2002.10164}.
  
Gloter, A. and Yoshida, N. (2021)
  Adaptive estimation for degenerate diffusion processes, \emph{Electronic Journal of Statistics}, \bold{15}, 1424--1472.
  
Gloter, A. and Yoshida, N. (2024)
  Non-adaptive estimation for degenerate diffusion processes, \emph{Theory of Probability and Mathematical Statistics}, \bold{110}, 75--99.
}
\author{
Yuta Koike with YUIMA Project Team
}
\note{
The model coefficients must be symbolically differentiable. 

If "x" is a component of \code{state.variable}, the character \code{"deltax"} should not be used as any variable name of the model since they are internally used for other variable names. 

Time-inhomogeneous diffusions are not supported.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qmle}}
}
\examples{
## dXt = theta3*Yt*dt
## dYt = -theta2*(Xt+Yt)dt + theta1*dWt

## Setting of the model
drift <- c("theta3*y", "-theta2*(x+y)")
diff.mat <- matrix(c(0, 0, 0, "theta1"), 2, 2)
ymod <- setModel(drift = drift, diffusion = diff.mat,
                 state.variable=c("x", "y"),solve.variable=c("x", "y"))
param <- list(theta1=1, theta2=2, theta3=0.7)

## Setting of the sampling scheme 
n <- 900
ysamp <- setSampling(Terminal = n/100, n = 10*n)
ysubsamp <- setSampling(Terminal = n/100, n = n)
yuima1 <- setYuima(model = ymod, sampling = ysamp)

## Simulation
set.seed(111)
yuima <- simulate(yuima1, xinit=c(1, 1), true.parameter = param)
yuima <- subsampling(yuima, ysubsamp)
plot(yuima)

## Adaptive one-step estimation 
start <- list(theta1=0.5, theta2=0.5, theta3=0.5)
lower <- list(theta1=0.1, theta2=0.1, theta3=0.1)
upper <- list(theta1=5, theta2=5, theta3=5)
res <- qmleDegenerate(yuima, start = start, lower = lower,upper = upper)
print(res)

\dontrun{
## Joint estimation
set.seed(123)
yuima <- simulate(yuima1, xinit=c(1, 1), true.parameter = param)
yuima <- subsampling(yuima, ysubsamp)
res <- qmleDegenerate(yuima, start = start, lower = lower,upper = upper, joint = TRUE)
print(res)

## Fix theta1 to the true value
res <- qmleDegenerate(yuima, start = start, lower = lower,upper = upper, 
                       fixed = list(theta1 =1))
print(res)

## Multi-dimensional case
drift <- c("a1*z", "a2*(z+w)", "-b1*x-z", "-b2*(x+y)-w")
diff.mat <- matrix(c(0,0,0,0,0,0,0,0,0,0,0, "sigma1", 0, 0,"sigma2", "sigma3"), 4, 4,byrow=TRUE)
ymod <- setModel(drift = drift, diffusion = diff.mat,
                 state.variable=c("x", "y", "z", "w"),solve.variable=c("x", "y", "z", "w"))
param <- list(a1=1,a2=0.5,b1=1,b2=0.5,sigma1=1,sigma2=0.5,sigma3=0.5)

n <- 100
ysamp <- setSampling(Terminal = n/10, n = 10*n)
ysubsamp <- setSampling(Terminal = n/10, n = n)
yuima1 <- setYuima(model = ymod, sampling = ysamp)

set.seed(111)
yuima <- simulate(yuima1, xinit=c(1, 1, 1, 1), true.parameter = param)
yuima <- subsampling(yuima, ysubsamp)
plot(yuima)

start <- list(a1=2,a2=2,b1=2,b2=2,sigma1=2,sigma2=2,sigma3=2)
lower <- list(a1=0.1,a2=0.1,b1=0.1,b2=0.1,sigma1=0.1,sigma2=0.1,sigma3=0.1)
upper <- list(a1=5,a2=5,b1=5,b2=5,sigma1=5,sigma2=5,sigma3=5)
res <- qmleDegenerate(yuima, start = start, lower = lower,upper = upper)
print(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ts}
