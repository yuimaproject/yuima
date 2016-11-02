\name{Diagnostic.Cogarch}
\alias{Diagnostic.Cogarch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for checking the statistical properties of the COGARCH(p,q) model
}
\description{
The function check the statistical properties of the COGARCH(p,q) model. We verify if the process has a strict positive stationary variance model.}
\usage{
Diagnostic.Cogarch(yuima.cogarch, param = list(), matrixS = NULL, mu = 1, display = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima.cogarch}{ an object of class \code{yuima.cogarch}, \code{yuima} or a class \code{cogarch.gmm-class}
}
  \item{param}{ a list containing the values of the parameters}
  \item{matrixS}{ a Square matrix.}
  \item{mu}{ first moment of the Levy measure.}
  \item{display}{ a logical variable, if \code{TRUE} the function displays the result in the \code{console}. }
}

\value{The functon returns a List with entries:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{meanVarianceProc }{ Unconditional Stationary mean of the variance process. }
\item{meanStateVariable}{ Unconditional Stationary mean of the state process.}
\item{stationary}{ If \code{TRUE}, the COGARCH(p,q) has stationary variance.}
\item{positivity}{ If \code{TRUE}, the variance process is strictly positive.}
}
\author{YUIMA Project Team}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
# Definition of the COGARCH(1,1) process driven by a Variance Gamma nois:
param.VG <- list(a1 = 0.038,  b1 =  0.053,
                  a0 = 0.04/0.053,lambda = 1, alpha = sqrt(2), beta = 0, mu = 0, 
                  x01 = 50.33)

cog.VG <- setCogarch(p = 1, q = 1, work = FALSE,
                      measure=list(df="rvgamma(z, lambda, alpha, beta, mu)"),
                      measure.type = "code", 
                      Cogarch.var = "y",
                      V.var = "v", Latent.var="x",
                      XinExpr=TRUE)

# Verify the stationarity and the positivity of th variance process

test <- Diagnostic.Cogarch(cog.VG,param=param.VG)
show(test)

# Simulate a sample path

set.seed(210)

Term=800
num=24000

samp.VG <- setSampling(Terminal=Term, n=num)

sim.VG <- simulate(cog.VG,
                    true.parameter=param.VG,
                    sampling=samp.VG,
                    method="euler")
plot(sim.VG)

# Estimate the model

res.VG <- gmm(sim.VG, start = param.VG, Est.Incr = "IncrPar")

summary(res.VG)

#  Check if the estimated COGARCH(1,1) has a positive and stationary variance

test1<-Diagnostic.Cogarch(res.VG)
show(test1)

# Simulate a COGARCH sample path using the estimated COGARCH(1,1) 
# and the recovered increments of underlying Variance Gamma Noise

esttraj<-simulate(res.VG)
plot(esttraj)


}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
