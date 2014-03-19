\name{llag}
\alias{llag}
\alias{llag,list-method}
\title{Lead Lag Estimator}
\description{Estimate the lead-lag parameters of discretely observed processes by maximizing the shifted Hayashi-Yoshida covariation contrast functions, following Hoffmann et al. (2013).
}
\usage{
llag(x,from=-Inf,to=Inf,division=FALSE,verbose=FALSE,grid)
}
\arguments{
  \item{x}{an object of  \code{\link{yuima-class}} or
	\code{\link{yuima.data-class}}.}
  \item{verbose}{whether \code{llag} returns matrices or not.}
\item{from}{a numeric vector each of whose component(s) indicates the lower end of a finite grid on which the contrast function is evaluated, if \code{grid} is missing.}
\item{to}{a numeric vector each of whose component(s) indicates the upper end of a finite grid on which the contrast function is evaluated, if \code{grid} is missing.}
\item{division}{a numeric vector each of whose component(s) indicates the number of the points of a finite grid on which the contrast function is evaluated, if \code{grid} is missing.}
\item{grid}{a numeric vector or a list of numeric vectors. See `Details'.}
}
\details{
Let \eqn{d} be the number of the components of the \code{zoo.data} of the object \code{x}.

Let \eqn{Xi_{ti_0},Xi_{ti_1},\dots,Xi_{ti_n(i)}} be the observation data of the \eqn{i}-th component (i.e. the \eqn{i}-th component of the \code{zoo.data} of the object \code{x}).

The shifted Hayashi-Yoshida covariation contrast function \eqn{Uij(\theta)} of the observations \eqn{Xi} and \eqn{Xj} \eqn{(i<j)} is defined by the same way as in Hoffmann et al. (2013). The lead-lag parameter \eqn{\theta_ij} is defined as a maximizer of \eqn{|Uij(\theta)|}. \eqn{Uij(\theta)} is evaluated on a finite grid \eqn{Gij} defined below. Thus \eqn{\theta_ij} belongs to this grid. If there exist more than two maximizers, the lowest one is selected.

The grid \eqn{Gij} is defined as follows. First, if \code{grid} is missing, \eqn{Gij} is given by  
\deqn{a, a+(b-a)/(N-1), \dots, a+(N-2)(b-a)/(N-1), b,}
where \eqn{a,b} and \eqn{N} are the \eqn{(d(i-1)-(i-1)i/2+(j-i))}-th components of \code{from}, \code{to} and \code{division} respectively.  If the corresponding component of \code{from} (resp. \code{to}) is \code{-Inf} (resp. \code{Inf}), \eqn{a=-(tj_n(j)-ti_0)} (resp. \eqn{b=ti_n(i)-tj_0}) is used, while if the corresponding component of \code{division} is \code{FALSE}, \eqn{N=round(2max(n(i),n(j)))+1} is used. Missing components are filled with \code{-Inf} (resp. \code{Inf}, \code{FALSE}). The default value \code{-Inf} (resp. \code{Inf}, \code{FALSE}) means that all components are \code{-Inf} (resp. \code{Inf}, \code{FALSE}). Next, if \code{grid} is a numeric vector, \eqn{Gij} is given by \code{grid}. If \code{grid} is a list of numeric vectors, \eqn{Gij} is given by the \eqn{(d(i-1)-(i-1)i/2+(j-i))}-th component of \code{grid}.    

The estimated lead-lag parameters are returned as the skew-symmetric matrix \eqn{(\theta_ij)_{i,j=1,\dots,d}}. If \code{verbose} is \code{TRUE}, the covariance matrix \eqn{(Uij(\theta_ij))_{i,j=1,\dots,d}} corresponding to the estimated lead-lag parameters, the corresponding correlation matrix and the computed contrast functions are also returned. The computed contrast functions are returned as a list with the length \eqn{d(d-1)/2}. For \eqn{i<j}, the \eqn{(d(i-1)-(i-1)i/2+(j-i))}-th component of the list consists of an object \eqn{Uij(\theta)} of class \link{zoo} indexed by \eqn{Gij}.}
\value{
If \code{verbose} is \code{FALSE}, a skew-symmetric matrix corresponding to the estimated lead-lag parameters is returned. If \code{verbose} is \code{TRUE}, a list with the following components is returned:
\item{lagcce}{a skew-symmetric matrix corresponding to the estimated lead-lag parameters.}
\item{covmat}{a covariance matrix corresponding to the estimated lead-lag parameters.}
\item{cormat}{a correlation matrix corresponding to the estimated lead-lag parameters.}
\item{crosscov}{a list of computed contrast functions.}
}
%\note{TBA
%}
\references{
Hoffmann, M., Rosenbaum, M. and Yoshida, N. (2013)
  Estimation of the lead-lag parameter from non-synchronous data, 
  \emph{Bernoulli}, \bold{19}, no. 2, 426--461.
}
\author{The YUIMA Project Team}
\examples{

## Set a model
diff.coef.matrix <- matrix(c("sqrt(x1)", "3/5*sqrt(x2)","1/3*sqrt(x3)", 
"", "4/5*sqrt(x2)","2/3*sqrt(x3)","","","2/3*sqrt(x3)"), 3, 3) 
drift <- c("1-x1","2*(10-x2)","3*(4-x3)")
cor.mod <- setModel(drift = drift, diffusion = diff.coef.matrix,
  solve.variable = c("x1", "x2","x3")) 

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
result <- llag(psample,verbose=TRUE)

#estimated lead-lag parameter
result$lagcce

#plotting the computed contrastfunctions
par(mfcol=c(3,1))
plot(result$crosscov[[1]]) #(1,2)-th component
plot(result$crosscov[[2]]) #(1,3)-th component
plot(result$crosscov[[3]]) #(2,3)-th component

}
\seealso{\link{cce}}
\keyword{ts}
