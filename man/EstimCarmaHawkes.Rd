\name{EstimCarmaHawkes}
\alias{qmle.CarmaHawkes}
\alias{MethodOfMoments.CarmaHawkes}
\alias{EstimCarmaHawkes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation Methods for a CARMA(p,q)-Hawkes Counting Process
}
\description{The function provides two estimation procedures: Maximum Likelihood Estimation and Matching Empirical Correlation}
\usage{
EstimCarmaHawkes(yuima, start, est.method = "qmle", method = "BFGS", 
lower = NULL, upper = NULL, lags = NULL, display = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{a yuima object.}
  \item{start}{initial values to be passed to the optimizer.}
  \item{est.method}{ The method used to estimate the parameters. The default \code{est.method = "qmle"} indicates the MLE while the alternative approach is based on the minimization of the empirical and theoretical autocorrelation.
%%     ~~Describe \code{est.method} here~~
}
  \item{method}{The optimization method to be used. See \code{\link{optim}}.
}
  \item{lower}{Lower Bounds.
}
  \item{upper}{
Upper Bounds.
}
  \item{lags}{
Number of lags used in the autocorrelation.
}
  \item{display}{
you can see a progress of the estimation when \code{display=TRUE}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The output contains the estimated parameters.}
\references{
Mercuri, L., Perchiazzo, A., & Rroji, E. (2022). A Hawkes model with CARMA (p, q) intensity. \doi{10.48550/arXiv.2208.02659}.
}
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
\examples{
\dontrun{
## MLE For A CARMA(2,1)-Hawkes ##

# Inputs:
a <- c(3,2)
b <- c(1,0.3)
mu<-0.30

true.par<-c(mu,a,b)

# step 1) Model Definition => Constructor 'setCarmaHawkes'
p <- 2
q <- 1
mod1 <- setCarmaHawkes(p = p,q = q)

# step 2) Grid Construction => Constructor 'setSampling'
FinalTime <- 5000
t0 <- 0
samp <- setSampling(t0, FinalTime, n = FinalTime)

# step 3) Simulation => method 'simulate'
# We use method 'simulate' to generate our dataset. 
# For the estimation from real data, 
# we use the constructors 'setData' and 
#'setYuima' (input 'model' is an object of 
#           'yuima.CarmaHawkes-class'). 

names(true.par) <- c(mod1@info@base.Int, mod1@info@ar.par, mod1@info@ma.par) 

set.seed(1)
system.time(
sim1 <- simulate(object = mod1, true.parameter = true.par, 
    sampling = samp)
)

plot(sim1)

# step 4) Estimation using the likelihood function.
system.time(
  res <- EstimCarmaHawkes(yuima = sim1, 
    start = true.par)
)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
