\name{setHawkes}
\alias{setHawkes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constructor of Hawkes model}
\description{'\code{setHawkes}' constructs an object of class \code{\link{yuima.Hawkes}} that is a mathematical description of a multivariate Hawkes model
}
\usage{
setHawkes(lower.var = "0", upper.var = "t", var.dt = "s",
  process = "N", dimension = 1, intensity = "lambda",
  ExpKernParm1 = "c", ExpKernParm2 = "a", const = "nu",
  measure = NULL, measure.type = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower.var}{Lower bound in the integral}
  \item{upper.var}{Upper bound in the integral}
  \item{var.dt}{Time variable}
  \item{process}{Counting process}
  \item{dimension}{An integer that indicates the components of the counting process}
  \item{intensity}{Intensity Process}
  \item{ExpKernParm1}{Kernel parameters}
  \item{ExpKernParm2}{Kernel parameters}
  \item{const}{Constant term in the intensity process}
  \item{measure}{Jump size. By default 1}
  \item{measure.type}{Type. By default \code{code}.}
}
\details{By default the object is an univariate Hawkes process}
\value{The function returns an object of class \code{\link{yuima.Hawkes}}.}
\author{YUIMA Team}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
# Definition of an univariate hawkes model

provaHawkes2<-setHawkes()
str(provaHawkes2)

# Simulation

true.par <- list(nu1=0.5, c11=3.5,  a11=4.5)

simprv1 <- simulate(object = provaHawkes2, true.parameter = true.par,
  sampling = setSampling(Terminal =70, n=7000))

plot(simprv1)

# Computation of intensity

lambda1 <- Intensity.PPR(simprv1, param = true.par)

plot(lambda1)

# qmle

res1 <- qmle(simprv1, method="Nelder-Mead", start = true.par)

summary(res1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
