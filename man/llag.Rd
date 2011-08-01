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

## Set a model
diff.coef.matrix <- matrix(c("sqrt(x1)", "3/5*sqrt(x2)","1/3*sqrt(x3)", 
"", "4/5*sqrt(x2)","2/3*sqrt(x3)","","","2/3*sqrt(x3)"), 3, 3) 
drift <- c("1-x1","2*(10-x2)","3*(4-x3)")
cor.mod <- setModel(drift = drift, diffusion = diff.coef.matrix,solve.variable = c("x1", "x2","x3")) 

set.seed(111) 

## We use a function poisson.random.sampling to get observation by Poisson sampling.
yuima.samp <- setSampling(Terminal = 1, n = 1200) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima,xinit=c(1,7,5)) 

# intentionally displace the second time series

  data2 <- yuima@data@zoo.data[[2]]
  time2 <- time(data2)
  theta2 <- 0.05   # the lag of x2 behind x1
  stime2 <- time2 + theta2  
  time(yuima@data@zoo.data[[2]]) <- stime2

  data3 <- yuima@data@zoo.data[[3]]
  time3 <- time(data3)
  theta3 <- 0.12   # the lag of x3 behind x1
  stime3 <- time3 + theta3 
  time(yuima@data@zoo.data[[3]]) <- stime3




#sampled data by Poisson rules
psample<- poisson.random.sampling(yuima, rate = c(0.2,0.3,0.4), n = 1000) 


# plot
plot(psample)


#cce
cce(psample)

#lead-lag estimation
llag(psample)

}
\seealso{\link{cce}}
\keyword{ts}
