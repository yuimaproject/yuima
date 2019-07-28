\name{setCogarch}
\alias{setCogarch}
\alias{Cogarch}
\alias{cogarch}
\alias{COGARCH}
\alias{CoGarch}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Continuous-time GARCH (p,q) process
}
\description{
\code{setCogarch} describes the Cogarch(p,q) model introduced in Brockwell et al. (2006):

\code{dGt = sqrt(Vt)dZt}

\code{Vt = a0 + (a1 Yt(1) + ... + a(p) Yt(p))}

\code{dYt(1) = Yt(2) dt}  
  
  \code{...}
  
  \code{dYt(q-1) = Yt(q) dt}
  
  \code{dYt(q) = (-b(q) Yt(1) - ... - b(1) Yt(q))dt + (a0 + (a1 Yt(1) + ... + a(p) Yt(p))d[ZtZt]^{q}}


}
\usage{
setCogarch(p, q, ar.par = "b", ma.par = "a", loc.par = "a0", Cogarch.var = "g",
   V.var = "v", Latent.var = "y", jump.variable = "z",  time.variable = "t",
   measure = NULL, measure.type = NULL, XinExpr = FALSE, startCogarch = 0,
   work = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{a non negative integer that is the number of the moving average coefficients of the Variance process.}
  \item{q}{a non-negative integer that indicates the number of the autoregressive  coefficients of the Variance process.}
  \item{ar.par}{a character-string that is the label of the autoregressive coefficients.}
  \item{ma.par}{a character-string that is the label of the autoregressive coefficients.}
  \item{loc.par}{the location coefficient.}
  \item{Cogarch.var}{a character-string that is the label of the observed cogarch process.}
  \item{V.var}{a character-string that is the label of the latent variance process.}
  \item{Latent.var}{a character-string that is the label of the latent process in the state space representation for the variance process.}
  \item{jump.variable}{the jump variable.}
  \item{time.variable}{the time variable.}
  \item{measure}{Levy measure of jump variables.}
  \item{measure.type}{type specification for Levy measure.}
  \item{XinExpr}{a vector of \code{expressions} identyfying the starting conditions for Cogarch model.}
  \item{startCogarch}{Start condition for the Cogarch process }
  \item{work}{ Internal Variable. In the final release this input will be removed.}
  \item{\dots}{Arguments to be passed to \code{setCogarch} such as the slots of the \code{\link{yuima.model-class}}}
}
\details{
We remark that \code{yuima} describes a Cogarch(p,q) model using the formulation proposed in Brockwell et al. (2006). This representation has the Cogarch(1,1) model introduced in Kluppelberg et al. (2004) as a special case. Indeed, by choosing \code{beta = a0 b1, eta = b1} and \code{phi = a1}, we obtain the Cogarch(1,1) model proposed in Kluppelberg et al. (2004) defined as the solution of the SDEs:

\code{dGt = sqrt(Vt)dZt}

\code{dVt = (beta - eta Vt) dt + phi Vt d[ZtZt]^{q}}

Please refer to the vignettes and the examples.


An object of \code{\link{yuima.cogarch-class}} contains:

\describe{    
\item{\code{info}:}{It  is an object 
 of  \code{\link{cogarch.info-class}} which is a list of arguments that identifies the Cogarch(p,q) model}
}
and the same slots in an object of \code{\link{yuima.model-class}} .
}
\value{
  \item{model}{an object of  \code{\link{yuima.cogarch-class}}.}
}
\author{The YUIMA Project Team}
\note{There may be missing information in the model description. 
Please contribute with suggestions and fixings.
}

\references{
               Brockwell, P., Chadraa, E. and Lindner, A. (2006)
Continuous-time GARCH processes, \emph{The Annals of Applied Probability}, \bold{16}, 790-826.

Kluppelberg, C., Lindner, A., and Maller, R. (2004) A continuous-time GARCH process driven by a Levy process: Stationarity and second-order behaviour, \emph{Journal of Applied Probability}, \bold{41}, 601-622.

Stefano M. Iacus, Lorenzo Mercuri, Edit Rroji (2017) COGARCH(p,q): Simulation and Inference with the yuima Package, \emph{Journal of Statistical Software}, \bold{80}(4), 1-49.
}

\examples{
# Ex 1. (Continuous time GARCH process driven by a compound poisson process)
prova<-setCogarch(p=1,q=3,work=FALSE,
                  measure=list(intensity="1", df=list("dnorm(z, 0, 1)")),
                  measure.type="CP",
                  Cogarch.var="y",
                  V.var="v",
                  Latent.var="x")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{COGARCH}
\keyword{yuima.cogarch}% __ONLY ONE__ keyword per line
