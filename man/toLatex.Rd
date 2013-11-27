\name{toLatex}
\title{Additional Methods for LaTeX Representations for Yuima objects}
\alias{toLatex}
\alias{toLatex.yuima}
\alias{toLatex.yuima.model}
\alias{toLatex.yuima.carma}
\description{
  These methods convert \code{\link{yuima-class}},
   \code{\link{yuima.model-class}} or \code{\link{yuima.carma-class}} objects to character vectors with
   LaTeX markup.
}
\usage{
\method{toLatex}{yuima}(object,\dots)
\method{toLatex}{yuima.model}(object,\dots)
\method{toLatex}{yuima.carma}(object,\dots)
}
\arguments{
  \item{object}{object of a class yuima, yuima.model or yuima.carma.}
  \item{\dots}{currently not used.}
}
\details{
  This methos tries to convert a formal description of the model slot
  of the yuima object into a LaTeX formula.
  This is just a simple proof of concept and probably further LaTex 
  manipulations for use in papers.
  Copy and paste of the output of \code{toLatex} into a real LaTeX file
  should do the job.
} 
\examples{
# dXt = theta*Xt*dt + dWt
mod1 <- setModel(drift="theta*y", diffusion=1, solve.variable=c("y"))
str(mod1)
toLatex(mod1)

# A multi-dimensional (correlated) diffusion process. 
# To describe the following model: 
# X=(X1,X2,X3); dXt = U(t,Xt)dt + V(t)dWt
# For drift coeffcient
U <- c("-x1","-2*x2","-t*x3")
# For diffusion coefficient of X1 
v1 <- function(t) 0.5*sqrt(t)
# For diffusion coefficient of X2
v2 <- function(t) sqrt(t)
# For diffusion coefficient of X3
v3 <- function(t) 2*sqrt(t)
# correlation
rho <- function(t) sqrt(1/2)
# coefficient matrix for diffusion term
V <- matrix( c( "v1(t)",
                "v2(t) * rho(t)",
                "v3(t) * rho(t)",
                "",
                "v2(t) * sqrt(1-rho(t)^2)",
                "",
                "",
                "",
                "v3(t) * sqrt(1-rho(t)^2)" 
               ), 3, 3)
# Model sde using "setModel" function
cor.mod <- setModel(drift = U, diffusion = V,
state.variable=c("x1","x2","x3"), 
solve.variable=c("x1","x2","x3") )
str(cor.mod)
toLatex(cor.mod)

# A CARMA(p=3,q=1) process.
carma1<-setCarma(p=3,q=1,loc.par="c",scale.par="s")
str(carma1)
toLatex(carma1)
}
\keyword{misc}
