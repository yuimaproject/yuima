\name{llag}
\alias{llag}
\alias{llag,list-method}
\title{Lead Lag Estimator}
\description{TBA
}
\usage{
llag(x,verbose=FALSE)
}
\arguments{
  \item{x}{an object of  \code{\link{yuima-class}} or
	\code{\link{yuima.data-class}}.}
  \item{verbose}{whether llag returns matrix or not.}
}
\details{
  TBA
}
\value{
  \item{value}{estimate of the lag, covariance and correlation matrices}
}
\note{TBA
}
\references{
  TBA
}
\author{The YUIMA Project Team}
\examples{


# generate two sample paths

  rho <- 0.4   # correlation coefficiant

  diff.coef.matrix <- matrix(c(1,rho,rho,1),2,2) 
  drift <- c("0*x1","0*x2")
  cor.mod <- setModel(drift = drift, diffusion = diff.coef.matrix,solve.variable = c("x1", "x2")) 

  yuima.samp <- setSampling(Terminal = 1, n = 300) 
  yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
  yuima <- simulate(yuima) 

# intentionally displace the second time series

  data2 <- yuima@data@zoo.data[[2]]
  time <- time(data2)

  theta <- 0.05   # the degree of lag

  stime <- time + theta  
  time(yuima@data@zoo.data[[2]]) <- stime

# estimate lead-lag

 lagcce <- function(theta){
    stime <- time2 + theta  #shifted time
    time(dat@zoo.data[[2]]) <- stime
    return(cce(dat))
  }

# true value must be -theta=0.05
ret <- llag(yuima,verbose=TRUE)
print(ret$opt)
print(ret$covmat)
print(ret$cormat)

plot(ret$mat,xlab="theta",ylab="U(theta)",type="l")


}
\seealso{\link{cce}}
\keyword{ts}
