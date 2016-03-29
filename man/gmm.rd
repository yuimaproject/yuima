\name{gmm}
\alias{gmm}
\alias{gmm.COGARCH}
\alias{Method of Moment COGARCH}
%\alias{Estimation COGARCH}
%\alias{est.COGARCH}
%\alias{est.yuima Cog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Method of Moments for COGARCH(P,Q).
}
\description{
The function returns the estimated parameters of a COGARCH(P,Q) model. The parameters are abtained by matching theoretical vs empirical autocorrelation function. The theoretical autocorrelation function is computed according the methodology developed in Chadraa (2009).
}
\usage{
gmm(yuima, data = NULL, start,
 method="BFGS", fixed = list(), lower, upper, lag.max = NULL,
 equally.spaced = FALSE, aggregation=TRUE, Est.Incr = "NoIncr", objFun = "L2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{a yuima object or an object of \code{\link{yuima.cogarch-class}}.}
  \item{data}{an object of class \code{\link{yuima.data-class}} contains the observations available at uniformly spaced time. If \code{data=NULL}, the default, the function uses the data in an object of \code{\link{yuima-class}}.}
  \item{start}{a \code{list} containing the starting values for the optimization routine.}
  \item{method}{a string indicating one of the methods available in \code{\link{optim}}.}
  \item{fixed}{a list of fixed parameters in optimization routine.}
  \item{lower}{a named list for specifying lower bounds of parameters.}
  \item{upper}{a named list for specifying upper bounds of parameters.}
  \item{lag.max}{maximum lag at which to calculate the theoretical and empirical acf. Default is \code{sqrt{N}} where \code{N} is the number of observation.}
  \item{equally.spaced}{Logical variable. If \code{equally.spaced = TRUE.}, the function use the returns of COGARCH(P,Q) evaluated at unitary length for the computation of the empirical autocorrelations. If \code{equally.spaced = FALSE}, the increments are evaluated on the interval with frequency specified in an object of class \code{\link{yuima.data-class}} that contains the observed time series.}
  \item{aggregation}{If \code{aggregation=TRUE}, before the estimation of the levy parameters we aggregate the estimated increments}
  \item{Est.Incr}{ a string variable, If \code{Est.Incr = "NoIncr"}, default value, \code{gmm} returns an object of class  \code{\link{cogarch.est-class}} that contains the COGARCH parameters.
  If \code{Est.Incr = "Incr"} or \code{Est.Incr = "IncrPar"} the output is an object of class \code{\link{cogarch.est.incr-class}}. In the first case the object contains the increments of underlying noise while in the second case also the estimated parameter of levy measure.}
  \item{objFun}{a string variable that indentifies the objective function in the optimization step. \code{objFun = "L2"}, default value, the objective function is  a quadratic form where the weighting Matrix is the identity one. \code{objFun = "L2CUE"} the weighting matrix is estimated using Continuously Updating GMM (L2CUE).
  \code{objFun = "L1"}, the objective function is the mean absolute error. In the last case the standard error for estimators are not available.
  }
}
\details{
The routine is based on three steps: estimation of the COGARCH parameters, recovering the increments of the underlying Levy process and estimation of the levy measure parameters. The last two steps are available on request by the user.
}
\value{ The function returns a list with the same components of the object obtained when the function  \code{\link{optim}} is used.
}
\references{
Chadraa, E. (2009) Statistical Modeling with COGARCH(P,Q) Processes. Phd Thesis
}
\author{
The YUIMA Project Team.
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
# Example COGARCH(1,1): the parameters are the same used in Haugh et al. 2005. In this case
# we assume the underlying noise is a symmetric variance gamma.
# As first step we define the COGARCH(1,1) in yuima:

mod1 <- setCogarch(p = 1, q = 1, work = FALSE,
                   measure=list(df="rbgamma(z,1,sqrt(2),1,sqrt(2))"),
                    measure.type = "code", Cogarch.var = "y",
                    V.var = "v", Latent.var="x",XinExpr=TRUE)

param <- list(a1 = 0.038,  b1 =  0.053,
              a0 = 0.04/0.053, x01 = 20)

# We generate a trajectory
samp <- setSampling(Terminal=10000, n=100000)
set.seed(210)
sim1 <- simulate(mod1, sampling = samp, true.parameter = param)

# We estimate the model

res1 <- gmm(yuima = sim1, start = param)

summary(res1)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Method of Moments }
\keyword{ Estimation COGARCH }% __ONLY ONE__ keyword per line
