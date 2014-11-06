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
\code{setCogarch} describes the following model:

\code{dGt = sqrt(Vt)dZt}

\code{Vt = mu0 + (b0 Yt(0) + ... + b(q) Yt(q))}

\code{dYt(0) = Yt(1) dt}  
  
  \code{...}
  
  \code{dYt(p-2) = Yt(p-1) dt}
  
  \code{dYt(p-1) = (-a(p) Yt(0) - ... - a(1) Yt(p-1))dt + (mu0 + (b0 Yt(0) + ... + b(q) Yt(q))d[ZtZt]_{q}}


}
\usage{
setCogarch(p, q, ar.par = "a", ma.par = "b", loc.par = "mu", Cogarch.var = "g", V.var = "v", Latent.var = "y", jump.variable = "z",  time.variable = "t", measure = NULL, measure.type = NULL, XinExpr = FALSE, startCogarch = 0, work = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{a non negative integer that is the number of the autoregressive coefficients.}
  \item{q}{a non-negative integer that indicates the number of the moving average coefficients.}
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
  \item{XinExpr}{a vector of \code{expressions} identyfying the starting conditions for CaGarch model.}
  \item{startCogarch}{Start condition for the CoGarch process }
  \item{work}{ Internal Variable. In the final release this input will be removed.}
  \item{\dots}{Arguments to be passed to \code{setCogarch} such as the slots of the \code{\link{yuima.model-class}}}
}
\details{
Please refer to the vignettes and the examples or to the \pkg{yuimadocs} package.

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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
