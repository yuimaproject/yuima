\name{adaBayes}
\alias{adaBayes}
\alias{adaBayes,yuima-method}
\title{Adaptive Bayes estimator for the parameters in sde model}
\description{Adaptive Bayes estimator for the parameters in a specific type of sde.}
\usage{
adaBayes( yuima , h , prior2 , prior1 , domain2 , domain1, theta2.offset=0,theta1.offset=0)
}
\arguments{
  \item{yuima}{a 'yuima' object.}
  \item{h}{ sampling interval [positive value].}
  \item{prior2}{ prior distribution function of theta2 [function].}
  \item{prior1}{ prior distribution function of theta1 [function].}
  \item{domain2}{ domain of theta2 [vector or matrix].}
  \item{domain1}{ domain of theta1 [vector or matrix].}
  \item{theta2.offset}{ offset value of theta2 [scalar or vector].}
  \item{theta1.offset}{ offset value of theta1 [scalar or vector].}
}
\details{
  Calculate the values of the parameters theta1 and theta2. When the quasi-likelihood is too large for the given data, the integral in the Bayes estimator may diverges. For such a case, offset values for the parameters can be specified not to diverge.}
\value{
  \item{vector}{a vector of the paramter estimate}
}
\author{YUIMA Project Team}
\note{
This routine is not stable and accurate. Dr.Kamatani is now working for improvements.
}
\examples{
  set.seed(1)
  diff.matrix <- matrix(c("theta1"),1,1)
  mod <- setModel(drift=c("(-1)*theta2*x"),diffusion=diff.matrix,time.variable="t",state.variable="x")
  division <- 100
  Terminal <- (division)^(1/3)
  mod.sampling <- setSampling(Terminal=Terminal, division=division)
  yuima.mod <- setYuima(model=mod, sampling=mod.sampling)
  yuima.mod <- simulate(yuima=yuima.mod,
               xinit=1,
               true.parameter=c(0.6,0.3)
               )
  domain <- c(0,1)
  ## prior distribution of theta1
  prior1 <- function(theta1,domain){
    tmp <- 1/(domain[1,2] - domain[1,1])
    return(tmp)
  }
  ## prior distribution of theta2
  prior2 <- function(theta2,domain){
    tmp <- 1/(domain[1,2]-domain[1,1])
    return( tmp )
  }
  
  rand <- matrix(runif(10,0,1),10,2)
  QL <- 0
  param <- numeric(2)
  for( i in 1:5){
    PAR <- ml.ql(yuima.mod,theta2= rand[i,1],theta1=rand[i,2],h=(division^(-2/3)))
    if(PAR$value > QL){
      QL <- PAR$value
      param <- PAR$par
    }
  }
  BE <- adaBayes(yuima.mod,h=(division-1)^(-2/3),prior2=prior2,prior1=prior1,domain2=domain,domain1=domain,theta2.offset=param[1],theta1.offset=param[2])
  print("True param")
  print("theta2 = .6, theta1 = .3")
  print("ML estimate")
  PAR$par
  print("Bayes estimate")
  BE
}
\keyword{ts}
