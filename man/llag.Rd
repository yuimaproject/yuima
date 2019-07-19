\name{llag}
\alias{llag}
\alias{llag,list-method}
\title{Lead Lag Estimator}
\description{Estimate the lead-lag parameters of discretely observed processes by maximizing the shifted Hayashi-Yoshida covariation contrast functions, following Hoffmann et al. (2013).
}
\usage{
llag(x, from = -Inf, to = Inf, division = FALSE, verbose = (ci || ccor), 
     grid, psd = TRUE, plot = ci, ccor = ci, ci = FALSE, alpha = 0.01, 
     fisher = TRUE, bw, tol = 1e-7)
}
\arguments{

\item{x}{an object of  \code{\link{yuima-class}} or \code{\link{yuima.data-class}}.}

\item{verbose}{whether \code{llag} returns matrices or not. The default is \code{FALSE}.}

\item{from}{a numeric vector each of whose component(s) indicates the lower end of a finite grid on which the contrast function is evaluated, if \code{grid} is missing.}

\item{to}{a numeric vector each of whose component(s) indicates the upper end of a finite grid on which the contrast function is evaluated, if \code{grid} is missing.}

\item{division}{a numeric vector each of whose component(s) indicates the number of the points of a finite grid on which the contrast function is evaluated, if \code{grid} is missing.}

\item{grid}{a numeric vector or a list of numeric vectors. See 'Details'.}

\item{psd}{logical. If \code{TRUE}, the estimated cross-correlation functions are converted to the interval [-1,1]. See 'Details'.}

\item{plot}{logical. If \code{TRUE}, the estimated cross-correlation functions are plotted. If \code{ci} is also \code{TRUE}, the pointwise confidence intervals (under the null hypothesis that the corresponding correlation is zero) are also plotted. The default is \code{FALSE}.}

\item{ccor}{logical. If \code{TRUE}, the estimated cross-correlation functions are returned. This argument is ignored if \code{verbose} is \code{FALSE}. The default is \code{FALSE}.}

\item{ci}{logical. If \code{TRUE}, (pointwise) confidence intervals of the estimated cross-correlation functions and p-values for the significance of the correlations at the estimated lead-lag parameters are calculated. Note that the confidence intervals are only plotted when \code{plot=TRUE}.}

\item{alpha}{a posive number indicating the significance level of the confidence intervals for the cross-correlation functions.}

\item{fisher}{logical. If \code{TRUE}, the p-values and the confidence intervals for the cross-correlation functions is evaluated after applying the Fisher z transformation. This argument is only meaningful if \code{pval = "corr"}.}

\item{bw}{bandwidth parameter to compute the asymptotic variances. See 'Details' and \code{\link{hyavar}} for details.}

\item{tol}{tolelance parameter to avoid numerical errors in comparison of time stamps. All time stamps are divided by \code{tol} and rounded to integers. Note that the values of \code{grid} are also divided by \code{tol} and rounded to integers. A reasonable choice of \code{tol} is the minimum unit of time stamps. The default value \code{1e-6} supposes that the minimum unit of time stamps is greater or equal to 1 micro-second.}
}

\details{
Let \eqn{d} be the number of the components of the \code{zoo.data} of the object \code{x}.

Let \eqn{Xi_{ti_0},Xi_{ti_1},\dots,Xi_{ti_n(i)}} be the observation data of the \eqn{i}-th component (i.e. the \eqn{i}-th component of the \code{zoo.data} of the object \code{x}).

The shifted Hayashi-Yoshida covariation contrast function \eqn{Uij(\theta)} of the observations \eqn{Xi} and \eqn{Xj} \eqn{(i<j)} is defined by the same way as in Hoffmann et al. (2013), which corresponds to their cross-covariance function. The lead-lag parameter \eqn{\theta_ij} is defined as a maximizer of \eqn{|Uij(\theta)|}. \eqn{Uij(\theta)} is evaluated on a finite grid \eqn{Gij} defined below. Thus \eqn{\theta_ij} belongs to this grid. If there exist more than two maximizers, the lowest one is selected.

If \code{psd} is \code{TRUE}, for any \eqn{i,j} the matrix \eqn{C:=(Ukl(\theta))_{k,l\in{i,j}}} is converted to \code{(C\%*\%C)^(1/2)} for ensuring the positive semi-definiteness, and \eqn{Uij(\theta)} is redefined as the \eqn{(1,2)}-component of the converted \eqn{C}. Here, \eqn{Ukk(\theta)} is set to the realized volatility of \eqn{Xk}. In this case \eqn{\theta_ij} is given as a maximizer of the cross-correlation functions.

The grid \eqn{Gij} is defined as follows. First, if \code{grid} is missing, \eqn{Gij} is given by  
\deqn{a, a+(b-a)/(N-1), \dots, a+(N-2)(b-a)/(N-1), b,}
where \eqn{a,b} and \eqn{N} are the \eqn{(d(i-1)-(i-1)i/2+(j-i))}-th components of \code{from}, \code{to} and \code{division} respectively.  If the corresponding component of \code{from} (resp. \code{to}) is \code{-Inf} (resp. \code{Inf}), \eqn{a=-(tj_n(j)-ti_0)} (resp. \eqn{b=ti_n(i)-tj_0}) is used, while if the corresponding component of \code{division} is \code{FALSE}, \eqn{N=round(2max(n(i),n(j)))+1} is used. Missing components are filled with \code{-Inf} (resp. \code{Inf}, \code{FALSE}). The default value \code{-Inf} (resp. \code{Inf}, \code{FALSE}) means that all components are \code{-Inf} (resp. \code{Inf}, \code{FALSE}). Next, if \code{grid} is a numeric vector, \eqn{Gij} is given by \code{grid}. If \code{grid} is a list of numeric vectors, \eqn{Gij} is given by the \eqn{(d(i-1)-(i-1)i/2+(j-i))}-th component of \code{grid}.    

The estimated lead-lag parameters are returned as the skew-symmetric matrix \eqn{(\theta_ij)_{i,j=1,\dots,d}}. If \code{verbose} is \code{TRUE}, the covariance matrix \eqn{(Uij(\theta_ij))_{i,j=1,\dots,d}} corresponding to the estimated lead-lag parameters, the corresponding correlation matrix and the computed contrast functions are also returned. If further \code{ccor} is \code{TRUE},the computed cross-correlation functions are returned as a list with the length \eqn{d(d-1)/2}. For \eqn{i<j}, the \eqn{(d(i-1)-(i-1)i/2+(j-i))}-th component of the list consists of an object \eqn{Uij(\theta)/sqrt(Uii(\theta)*Ujj(\theta))} of class \link{zoo} indexed by \eqn{Gij}.

If \code{plot} is \code{TRUE}, the computed cross-correlation functions are plotted sequentially.

If \code{ci} is \code{TRUE}, the asymptotic variances of the cross-correlations are calculated at each point of the grid by using the naive kernel approach descrived in Section 8.2 of Hayashi and Yoshida (2011). The implementation is the same as that of \code{\link{hyavar}} and more detailed description is found there.
}

\value{
If \code{verbose} is \code{FALSE}, a skew-symmetric matrix corresponding to the estimated lead-lag parameters is returned. 
Otherwise, an object of class \code{"yuima.llag"}, which is a list with the following components, is returned:
\item{lagcce}{a skew-symmetric matrix corresponding to the estimated lead-lag parameters.}
\item{covmat}{a covariance matrix corresponding to the estimated lead-lag parameters.}
\item{cormat}{a correlation matrix corresponding to the estimated lead-lag parameters.}
\item{LLR}{a matrix consisting of lead-lag ratios. See Huth and Abergel (2014) for details.}

If \code{ci} is \code{TRUE}, the following component is added to the returned list:
\item{p.values}{a matrix of p-values for the significance of the correlations corresponding to the estimated lead-lag parameters.}

If further \code{ccor} is \code{TRUE}, the following components are added to the returned list:
\item{ccor}{a list of computed cross-correlation functions.}
\item{avar}{a list of computed asymptotic variances of the cross-correlations (if \code{ci = TRUE}).}
}

\note{
The default \code{grid} usually contains too many points, so it is better for users to specify this argument in order to reduce the computational time. See 'Examples' below for an example of the specification.

The evaluated p-values should carefully be interpreted because they are calculated based on \emph{pointwise confidence intervals} rather than \emph{simultaneous confidence intervals} (so there would be a multiple testing problem). Evaluation of p-values based on the latter will be implemented in the future extension of this function: Indeed, so far no theory has been developed for this. However, it is conjectured that the error distributions of the estimated cross-correlation functions are asymptotically independent if the grid is not dense too much, so p-values evaluated by this function will still be meaningful as long as sufficiently low significance levels are used.
}

\references{
Hayashi, T. and Yoshida, N. (2011)
  Nonsynchronous covariation process and limit theorems,
  \emph{Stochastic processes and their applications}, \bold{121}, 2416--2454.
  
Hoffmann, M., Rosenbaum, M. and Yoshida, N. (2013)
  Estimation of the lead-lag parameter from non-synchronous data, 
  \emph{Bernoulli}, \bold{19}, no. 2, 426--461.
  
Huth, N. and Abergel, F. (2014)
  High frequency lead/lag relationships --- Empirical facts,
  \emph{Journal of Empirical Finance}, \bold{26}, 41--58.
}
\author{Yuta Koike with YUIMA Project Team}
\examples{

## Set a model
diff.coef.matrix <- matrix(c("sqrt(x1)", "3/5*sqrt(x2)",
 "1/3*sqrt(x3)", "", "4/5*sqrt(x2)","2/3*sqrt(x3)",
 "","","2/3*sqrt(x3)"), 3, 3) 
drift <- c("1-x1","2*(10-x2)","3*(4-x3)")
cor.mod <- setModel(drift = drift, 
 diffusion = diff.coef.matrix,
 solve.variable = c("x1", "x2","x3")) 

set.seed(111) 

## We use a function poisson.random.sampling 
## to get observation by Poisson sampling.
yuima.samp <- setSampling(Terminal = 1, n = 1200) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima,xinit=c(1,7,5)) 

## intentionally displace the second time series

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




## sampled data by Poisson rules
psample<- poisson.random.sampling(yuima, 
 rate = c(0.2,0.3,0.4), n = 1000) 


## plot
plot(psample)


## cce
cce(psample)

## lead-lag estimation (with cross-correlation plots)
par(mfcol=c(3,1))
result <- llag(psample, plot=TRUE)

## estimated lead-lag parameter
result

## computing pointwise confidence intervals
llag(psample, ci = TRUE)

## In practice, it is better to specify the grid because the default grid contains too many points.
## Here we give an example for how to specify it.

## We search lead-lag parameters on the interval [-0.1, 0.1] with step size 0.01 
G <- seq(-0.1,0.1,by=0.01)

## lead-lag estimation (with computing confidence intervals)
result <- llag(psample, grid = G, ci = TRUE)

## Since the true lead-lag parameter 0.12 between x1 and x3 is not contained
## in the searching grid G, we see that the corresponding cross-correlation 
## does not exceed the cofidence interval

## detailed output
## the p-value for the (1,3)-th component is high
result

## Finally, we can examine confidence intervals of other significant levels
## and/or without the Fisher z-transformation via the plot-method defined 
## for yuima.llag-class objects as follows
plot(result, alpha = 0.001)
plot(result, fisher = FALSE)

par(mfcol=c(1,1))

}
\seealso{\code{\link{cce}}, \code{\link{hyavar}}, \code{\link{mllag}}, \code{\link{llag.test}}}
\keyword{ts}
