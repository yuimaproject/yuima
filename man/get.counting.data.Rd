\name{get.counting.data}
\alias{get.counting.data}
\alias{NoisePPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract arrival times from an object of class \code{yuima.PPR}
}
\description{This function extracts arrival times from an object of class \code{\link{yuima.PPR}}.}
\usage{
get.counting.data(yuimaPPR,type="zoo")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuimaPPR}{An object of class \code{\link{yuima.PPR}}.}
  \item{type}{By default \code{type="zoo"} the function returns an object of class \code{zoo}. Other values are \code{yuima.PPR} and \code{matrix}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{By default the function returns an object of class zoo. The arrival times can be extracted by applying the method \code{index} to the output}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
##################
# Hawkes Process #
##################

# Values of parameters.
mu <- 2
alpha <- 4
beta <-5

# Law definition

my.rHawkes <- function(n){
  res <- t(t(rep(1,n)))
  return(res)
}

Law.Hawkes <- setLaw(rng = my.rHawkes)

# Point Process Definition

gFun <- "mu"
Kernel <- "alpha*exp(-beta*(t-s))"

modHawkes <- setModel(drift = c("0"), diffusion = matrix("0",1,1),
  jump.coeff = matrix(c("1"),1,1), measure = list(df = Law.Hawkes),
  measure.type = "code", solve.variable = c("N"),
  xinit=c("0"))

prvHawkes <- setPPR(yuima = modHawkes, counting.var="N", gFun=gFun,
  Kernel = as.matrix(Kernel), lambda.var = "lambda", 
  var.dx = "N", lower.var="0", upper.var = "t")

true.par <- list(mu=mu, alpha=alpha,  beta=beta)

set.seed(1)

Term<-70
n<-7000

# Simulation trajectory

time.Hawkes <-system.time(
  simHawkes <- simulate(object = prvHawkes, true.parameter = true.par,
     sampling = setSampling(Terminal =Term, n=n))
)

# Arrival times of the Counting process.

DataHawkes <- get.counting.data(simHawkes)
TimeArr <- index(DataHawkes)

##################################
# Point Process Regression Model #
##################################

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

Term<-100
seed<-1
n<-10000

true.parKern <- list(mu=mu, alpha=alpha, beta=beta)


set.seed(seed)
# set.seed(1)

time.simKern <-system.time(
  simprvKern <- simulate(object = prvKern, true.parameter = true.parKern,
                         sampling = setSampling(Terminal =Term, n=n))
)


plot(simprvKern,main ="Counting Process with covariates" ,cex.main=0.9)

# Arrival Times
CountVar <- get.counting.data(simprvKern)
TimeArr <- index(CountVar)


}
}
