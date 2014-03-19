\name{asymptotic_term}
\alias{asymptotic_term}
\alias{asymptotic_term,yuima-method}

\title{asymptotic expansion of the expected value of the functional}
\description{calculate the fisrt and second term of asymptotic expansion of the functional mean.}
\usage{
asymptotic_term(yuima, block=100, rho, g, expand.var="e")
}
\arguments{
  \item{yuima}{an yuima object containing model and functional.}
  \item{block}{the number of trapezoids for integrals.}
  \item{rho}{specify discounting factor in mean integral.}
  \item{g}{arbitrary measurable function for mean integral.}
  \item{expand.var}{default expand.var="e".}
}
\details{
Calculate the first and second term of asymptotic expansion of the expected value of the functional associated with a sde. The returned value d0 + epsilon * d1 is approximation of the expected value.
}
\value{
  \item{terms}{list of 1st and 2nd asymptotic terms, terms$d0 and terms$d1.}
}
\author{YUIMA Project Team}
\note{
we need to fix this routine.
}
\examples{
# to the Black-Scholes economy:
# dXt^e = Xt^e * dt + e * Xt^e * dWt
diff.matrix <- "x*e"
model <- setModel(drift = "x", diffusion = diff.matrix)
# call option is evaluated by averating
# max{ (1/T)*int_0^T Xt^e dt, 0}, the first argument is the functional of interest:
Terminal <- 1
xinit <- c(1)
f <- list( c(expression(x/Terminal)), c(expression(0)))
F <- 0
division <- 1000
e <- .3
yuima <- setYuima(model = model, sampling = setSampling(Terminal=Terminal, n=division))
yuima <- setFunctional( yuima, f=f,F=F, xinit=xinit,e=e)

# asymptotic expansion
rho <- expression(0)
F0 <- F0(yuima)
get_ge <- function(x,epsilon,K,F0){
  tmp <- (F0 - K) + (epsilon * x) 
  tmp[(epsilon * x) < (K-F0)] <- 0
  return( tmp )
}
g <- function(x) get_ge(x,epsilon=e,K=1,F0=F0)
set.seed(123)
asymp <- asymptotic_term(yuima, block=10, rho,g)
asymp
sum(asymp$d0 + e * asymp$d1)


### An example of multivariate case: Heston model
## a <- 1;C <- 1;d <- 10;R<-.1
## diff.matrix <- matrix( c("x1*sqrt(x2)*e", "e*R*sqrt(x2)",0,"sqrt(x2*(1-R^2))*e"), 2,2)
## model <- setModel(drift = c("a*x1","C*(10-x2)"), 
## diffusion = diff.matrix,solve.variable=c("x1","x2"),state.variable=c("x1","x2"))
## call option is evaluated by averating
## max{ (1/T)*int_0^T Xt^e dt, 0}, the first argument is the functional of interest:
##
## Terminal <- 1
## xinit <- c(1,1)
## 
## f <- list( c(expression(0), expression(0)),   
## c(expression(0), expression(0)) , c(expression(0), expression(0))  )
## F <- expression(x1,x2)
## 
## division <- 1000
## e <- .3
## 
## yuima <- setYuima(model = model, sampling = setSampling(Terminal=Terminal, n=division))
## yuima <- setFunctional( yuima, f=f,F=F, xinit=xinit,e=e)
## 
## rho <- expression(x1)
## F0 <- F0(yuima)
## get_ge <- function(x){
##  return( max(x[1],0))
## }
## g <- function(x) get_ge(x)
## set.seed(123)
## asymp <- asymptotic_term(yuima, block=10, rho,g)
## sum(asymp$d0 + e * asymp$d1)


}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
