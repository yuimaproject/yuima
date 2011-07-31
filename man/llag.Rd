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
diff.coef.matrix <- matrix(c("sqrt(x1)", "sqrt(3/5)*sqrt(x2)", 
"", "sqrt(4/5)*sqrt(x2)"), 2, 2) 
drift <- c("-x1","2*(10-x2)")
cor.mod <- setModel(drift = drift, diffusion = diff.coef.matrix,solve.variable = c("x1", "x2")) 

set.seed(111) 

## We use a function poisson.random.sampling to get observation by Poisson sampling.
yuima.samp <- setSampling(Terminal = 1, n = 1200) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima,xinit=c(1,7)) 

# intentionally displace the second time series

  data2 <- yuima@data@zoo.data[[2]]
  time <- time(data2)

  theta <- 0.05   # the degree of lag

  stime <- time + theta  
  time(yuima@data@zoo.data[[2]]) <- stime


#sampled data by Poisson rules
psample<- poisson.random.sampling(yuima, rate = c(0.2,0.3), n = 1000) 


# plot
plot(psample)


#cce
cce(psample)

#lead-lag estimation
llag(psample)

}
\seealso{\link{cce}}
\keyword{ts}
