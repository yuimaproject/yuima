\name{setCarma}
\alias{Carma}
\alias{CARMA}
\alias{setCarma}
%\alias{Continuous}

\title{
  Continuous Autoregressive Moving Average (p, q) model 
}

\description{
  'setCarma' describes the following model: 
  
  \code{Vt = c0 + sigma (b0 Xt(0) + ... + b(q) Xt(q))} 
  
  \code{dXt(0) = Xt(1) dt}  
  
  \code{...}
  
  \code{dXt(p-2) = Xt(p-1) dt}
  
  \code{dXt(p-1) = (-a(p) Xt(0) - ... - a(1) Xt(p-1))dt + (gamma(0) + gamma(1) Xt(0) + ... + gamma(p) Xt(p-1))dZt}
%  \code{dX_{q,t} = (-beta_{q} X_{1,t} - ... - beta_{1} X_{1,t})dt + dL_{t}}

The continuous ARMA process using the state-space representation as in Brockwell (2000) is obtained by choosing: 

  \code{ gamma(0) = 1, gamma(1) = gamma(2) = ... = gamma(p) = 0.} 
  

Please refer to the vignettes and the examples or the \bold{yuima} documentation for details.  

}

\usage{
setCarma(p,q,loc.par=NULL,scale.par=NULL,ar.par="a",ma.par="b",
lin.par=NULL,Carma.var="v",Latent.var="x",XinExpr=FALSE, Cogarch=FALSE, ...)
}

\arguments{
  \item{p}{a non-negative integer that indicates the number of the autoregressive coefficients.}
  \item{q}{a non-negative integer that indicates the number of the moving average coefficients.}
  \item{loc.par}{location coefficient. The default value \code{loc.par=NULL} implies that \code{c0=0}.}
  \item{scale.par}{scale coefficient. The default value \code{scale.par=NULL} implies that \code{sigma=1}.}
  \item{ar.par}{a character-string that is the label of the autoregressive coefficients. The default Value is \code{ar.par="a"}.}
  \item{ma.par}{a character-string that is the label of the moving average coefficients. The default Value is \code{ma.par="b"}.}
  \item{Carma.var}{a character-string that is the label of the observed process. Defaults to \code{"v"}.}
  \item{Latent.var}{a character-string that is the label of the unobserved process. Defaults to \code{"x"}.}
  \item{lin.par}{a character-string that is the label of the linear coefficients. If \code{lin.par=NULL}, the default, the 'setCarma' builds the CARMA(p, q) model defined as in Brockwell (2000). }
  \item{XinExpr}{a logical variable. The default value \code{XinExpr=FALSE} implies that the starting condition for \code{Latent.var} is zero. If \code{XinExpr=TRUE}, each component of \code{Latent.var} has a parameter as a initial value.}
  \item{Cogarch}{a logical variable. The default value \code{Cogarch=FALSE} implies that the parameters are specified according to Brockwell (2000).  }
  \item{...}{Arguments to be passed to 'setCarma', such as the slots of \code{\link{yuima.model-class}}
  \describe{
  \item{\code{measure}}{Levy measure of jump variables.}
  \item{\code{measure.type}}{type specification for Levy measure.}
  \item{\code{xinit}}{a vector of \code{expression}s identyfying the starting conditions for CARMA model.}
  }
  }
}

\details{
Please refer to the vignettes and the examples or to the \pkg{yuimadocs} package.

An object of \code{\link{yuima.carma-class}} contains:

\describe{  	
\item{\code{info}:}{It  is an object 
 of  \code{\link{carma.info-class}} which is a list of arguments that identifies the carma(p,q) model}
}
and the same slots in an object of \code{\link{yuima.model-class}} .
}
%\describe{		
%\item{\code{p}:}{an \R ... which specifies ... .}
%\item{\code{q}:}{an \R ... which specifies ... .}
%\item{\code{ma.par}:}{ ... .}
%\item{\code{ar.par}:}{ ... .}
%\item{\code{lin.par}:}{... .}
%\item{\code{...}:}{... \code{\link{setModel}}  
%... }
%  }

\value{
  \item{model}{an object of  \code{\link{yuima.carma-class}}.}
}

\author{The YUIMA Project Team}
\note{There may be missing information in the model description. 
Please contribute with suggestions and fixings.
}

 \references{
               Brockwell, P. (2000)
Continuous-time ARMA processes, \emph{Stochastic Processes: Theory and Methods. Handbook of Statistics}, \bold{19}, (C. R. Rao and D. N. Shandhag, eds.) 249-276. North-Holland, Amsterdam.
}

\examples{
# Ex 1. (Continuous ARMA process driven by a Brownian Motion)
# To describe the state-space representation of a CARMA(p=3,q=1) model:
# Vt=c0+alpha0*X0t+alpha1*X1t
# dX0t = X1t*dt 
# dX1t = X2t*dt
# dX2t = (-beta3*X0t-beta2*X1t-beta1*X2t)dt+dWt
# we set
mod1<-setCarma(p=3, 
               q=1, 
               loc.par="c0")
# Look at the model structure by
str(mod1)

# Ex 2. (General setCarma model driven by a Brownian Motion)
# To describe the model defined as:
# Vt=c0+alpha0*X0t+alpha1*X1t
# dX0t = X1t*dt 
# dX1t = X2t*dt
# dX2t = (-beta3*X0t-beta2*X1t-beta1*X2t)dt+(c0+alpha0*X0t)dWt
# we set 
mod2 <- setCarma(p=3,
                 q=1,
                 loc.par="c0",
                 ma.par="alpha",
                 ar.par="beta",
                 lin.par="alpha")
# Look at the model structure by
str(mod2)

# Ex 3. (Continuous Arma model driven by a Levy process)
# To specify the CARMA(p=3,q=1) model driven by a Compound Poisson process defined as:
# Vt=c0+alpha0*X0t+alpha1*X1t
# dX0t = X1t*dt 
# dX1t = X2t*dt
# dX2t = (-beta3*X0t-beta2*X1t-beta1*X2t)dt+dzt
# we set the Levy measure as in setModel
mod3 <- setCarma(p=3,
                 q=1,
                 loc.par="c0",
                 measure=list(intensity="1",df=list("dnorm(z, 0, 1)")),
                 measure.type="CP")
# Look at the model structure by
str(mod3)

# Ex 4. (General setCarma  model driven by a Levy process)
# Vt=c0+alpha0*X0t+alpha1*X1t
# dX0t = X1t*dt 
# dX1t = X2t*dt
# dX2t = (-beta3*X1t-beta2*X2t-beta1*X3t)dt+(c0+alpha0*X0t)dzt
mod4 <- setCarma(p=3,
                 q=1,
                 loc.par="c0",
                 ma.par="alpha",
                 ar.par="beta",
                 lin.par="alpha",
                 measure=list(intensity="1",df=list("dnorm(z, 0, 1)")),
                 measure.type="CP")
# Look at the model structure by
str(mod4)
}
\keyword{ts}


