\name{DataPPR}
\alias{DataPPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{From \code{zoo} data to \code{yuima.PPR}.}
\description{The function converts an object of class \code{zoo} to an object of class \code{\link{yuima.PPR}}.}
\usage{
DataPPR(CountVar, yuimaPPR, samp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CountVar}{An object of class \code{zoo} that contains counting variables and covariates. \code{index(CountVar)} returns the arrival times.}
  \item{yuimaPPR}{An object of class \code{\link{yuima.PPR}} that contains a mathematical description of the point process regression model assumed to be the generator of the observed data.}
  \item{samp}{An object of class \code{yuima.sampling}.}
}
\value{The function returns an object of class \code{\link{yuima.PPR}} where the slot \code{model} contains the Point process described in \code{yuimaPPR@model}, the slot \code{data} contains the counting variables and the covariates observed on the grid in \code{samp}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
# In this example we generate a dataset contains the Counting Variable N
# and the Covariate X.
# The covariate X is an OU driven by a Gamma process.

# Values of parameters.
mu <- 2
alpha <- 4
beta <-5

# Law definition
my.rKern <- function(n,t){
  res0 <- t(t(rgamma(n, 0.1*t)))
  res1 <- t(t(rep(1,n)))
  res <- cbind(res0,res1)
  return(res)
}

Law.PPRKern <- setLaw(rng = my.rKern)

# Point Process definition
modKern <- setModel(drift = c("0.4*(0.1-X)","0"),
                    diffusion = c("0","0"),
                    jump.coeff = matrix(c("1","0","0","1"),2,2),
                    measure = list(df = Law.PPRKern),
                    measure.type = c("code","code"),
                    solve.variable = c("X","N"),
                    xinit=c("0.25","0"))

gFun <- "exp(mu*log(1+X))"
#
Kernel <- "alpha*exp(-beta*(t-s))"

prvKern <- setPPR(yuima = modKern,
                  counting.var="N", gFun=gFun,
                  Kernel = as.matrix(Kernel),
                  lambda.var = "lambda", var.dx = "N",
                  lower.var="0", upper.var = "t")

# Simulation

Term<-200
seed<-1
n<-20000

true.parKern <- list(mu=mu, alpha=alpha, beta=beta)


set.seed(seed)
# set.seed(1)

time.simKern <-system.time(
  simprvKern <- simulate(object = prvKern, true.parameter = true.parKern,
                         sampling = setSampling(Terminal =Term, n=n))
)


plot(simprvKern,main ="Counting Process with covariates" ,cex.main=0.9)

# Using the function get.counting.data we extract from an object of class
# yuima.PPR the counting process N and the covariate X at the arrival times.

CountVar <- get.counting.data(simprvKern)

plot(CountVar)

# We convert the zoo object in the yuima.PPR object.

sim2 <- DataPPR(CountVar, yuimaPPR=simprvKern, samp=simprvKern@sampling)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
