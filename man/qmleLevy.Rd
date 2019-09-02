\encoding{UTF-8}
\name{qmleLevy}
\alias{qmleLevy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gaussian quasi-likelihood estimation for Levy driven SDE
}
\description{
Calculate the Gaussian quasi-likelihood and Gaussian quasi-likelihood estimators of Levy driven SDE.
}
\usage{
qmleLevy(yuima, start, lower, upper, joint = FALSE, third = FALSE)
}
\arguments{
  \item{yuima}{a yuima object.}
  \item{lower}{a named list for specifying lower bounds of parameters.}
  \item{upper}{a named list for specifying upper bounds of parameters.}
  \item{start}{initial values to be passed to the optimizer.}
  \item{joint}{perform joint estimation or two stage estimation? by default joint=FALSE. If there exists an overlapping parameter, joint=TRUE does not work for the theoretical reason}
  \item{third}{perform third estimation? by default third=FALSE. If there exists an overlapping parameter, third=TRUE does not work for the         
               theoretical reason.}
}
\details{
This function performs Gaussian quasi-likelihood estimation for Levy driven SDE.
}
\value{
\item{first}{estimated values of first estimation (scale parameters)}

\item{second}{estimated values of second estimation (drift parameters)}

\item{third}{estimated values of third estimation (scale parameters)}
}

\note{
The function \code{qmleLevy} uses the function \code{qmle} internally.
It can be applied only for the standardized Levy noise whose moments of any order exist.
In present \code{yuima} package, birateral gamma (bgamma) process, normal inverse Gaussian (NIG) process, variance gamma (VG) process, and normal tempered stable process are such candidates.
In the current version, the standardization condition on the driving noise is internally checked only for the one-dimensional noise.
The standardization condition for the multivariate noise is given in

\href{https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnx5dW1hdWVoYXJhMTkyOHxneDo1OGIxNGQ2YjBlYWIxNzA3}{https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnx5dW1hdWVoYXJhMTkyOHxneDo3ZTdlMTA1OTMyZTBkYjQ2}

or

\href{https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnx5dW1hdWVoYXJhMTkyOHxneDo3ZTdlMTA1OTMyZTBkYjQ2}{https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnx5dW1hdWVoYXJhMTkyOHxneDo3ZTdlMTA1OTMyZTBkYjQ2}.

They also contain more presice explanation of this function.
}
\references{

Masuda, H. (2013). Convergence of Gaussian quasi-likelihood random fields for ergodic Levy driven SDE observed at high frequency. The Annals of Statistics, 41(3), 1593-1641.

Masuda, H. and Uehara, Y. (2017). On stepwise estimation of Levy driven 
stochastic differential equation (Japanese) ., Proc. Inst. Statist. Math., accepted.
}
\author{
The YUIMA Project Team

Contacts: Yuma Uehara \email{y-uehara@math.kyushu-u.ac.jp}
}



%\seealso{
%}
\examples{
## One-dimensional case 
dri<-"-theta0*x" ## set drift
jum<-"theta1/(1+x^2)^(-1/2)" ## set jump
yuima<-setModel(drift = dri
                ,jump.coeff = jum
                ,solve.variable = "x",state.variable = "x"
                ,measure.type = "code"
                ,measure = list(df="rbgamma(z,1,sqrt(2),1,sqrt(2))")) ## set true model
n<-3000
T<-30 ## terminal
hn<-T/n ## stepsize

sam<-setSampling(Terminal = T, n=n) ## set sampling scheme
yuima<-setYuima(model = yuima, sampling = sam) ## model

true<-list(theta0 = 1,theta1 = 2) ## true values
upper<-list(theta0 = 4, theta1 = 4) ## set upper bound
lower<-list(theta0 = 0.5, theta1 = 1) ## set lower bound
set.seed(123)
yuima<-simulate(yuima, xinit = 0, true.parameter = true,sampling = sam) ## generate a path
start<-list(theta0 = runif(1,0.5,4), theta1 = runif(1,1,4)) ## set initial values
qmleLevy(yuima,start=start,lower=lower,upper=upper, joint = TRUE) 

## Multi-dimensional case

lambda<-1/2
alpha<-1
beta<-c(0,0)
mu<-c(0,0)
Lambda<-matrix(c(1,0,0,1),2,2) ## set parameters in noise

dri<-c("1-theta0*x1-x2","-theta1*x2")
jum<-matrix(c("x1*theta2+1","0","0","1"),2,2) ## set coefficients

yuima <- setModel(drift=dri, 
                 solve.variable=c("x1","x2"),state.variable = c("x1","x2"), 
                 jump.coeff=jum, measure.type="code",
                 measure=list(df="rvgamma(z, lambda, alpha, beta, mu, Lambda
                 )"))

n<-3000 ## the number of total samples
T<-30 ## terminal
hn<-T/n ## stepsize

sam<-setSampling(Terminal = T, n=n) ## set sampling scheme
yuima<-setYuima(model = yuima, sampling = sam) ## model

true<-list(theta0 = 1,theta1 = 2,theta2 = 3,lambda=lambda, alpha=alpha, 
beta=beta,mu=mu, Lambda=Lambda) ## true values
upper<-list(theta0 = 4, theta1 = 4, theta2 = 5, lambda=lambda, alpha=alpha, 
beta=beta,mu=mu, Lambda=Lambda) ## set upper bound
lower<-list(theta0 = 0.5, theta1 = 1, theta2 = 1, lambda=lambda, alpha=alpha, 
beta=beta,mu=mu, Lambda=Lambda) ## set lower bound
set.seed(123)
yuima<-simulate(yuima, xinit = c(0,0), true.parameter = true,sampling = sam) ## generate a path
plot(yuima)
start<-list(theta0 = runif(1,0.5,4), theta1 = runif(1,1,4), 
theta2 = runif(1,1,5),lambda=lambda, alpha=alpha, 
beta=beta,mu=mu, Lambda=Lambda) ## set initial values
qmleLevy(yuima,start=start,lower=lower,upper=upper,joint = FALSE,third=TRUE) 
}

\keyword{qmle}% use one of  RShowDoc("KEYWORDS")
\keyword{Estimation}% __ONLY ONE__ keyword per line
