\name{setCarmaHawkes}
\alias{CarmaHawkes}
\alias{Hawkes.Carma.Intensity}
\alias{setCarmaHawkes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hawkes Process with a Continuous Autoregressive Moving Average(p, q) Intensity
}
\description{
'setCarmaHawkes' describes a self-exciting Hawkes process where the intensity is a CARMA(p,q) process. The model admits the Hawkes process with exponential kernel as a special case but it is able to reproduce a more complex time-dependence structure
}
\usage{
setCarmaHawkes(p, q, law = NULL, base.Int = "mu0", 
ar.par = "a", ma.par = "b", Counting.Process = "N", 
Intensity.var = "lambda", Latent.var = "x", time.var = "t", 
Type.Jump = FALSE, XinExpr = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{a non-negative integer that indicates the number of the autoregressive coefficients.
}
  \item{q}{a non-negative integer that indicates the number of the moving average coefficients.
}
  \item{law}{ An object of \code{\link{yuima.law-class}} that describes the jump size. If it is \code{NULL}, the jump size is unitary.
}
  \item{base.Int}{a character-string that is the label of the baseline Intensity parameter. Defaults to \code{"mu0"}.
%%     ~~Describe \code{base.Int} here~~
}
  \item{ar.par}{a character-string that is the label of the autoregressive coefficients. The default Value is \code{ar.par="a"}.
}
  \item{ma.par}{a character-string that is the label of the moving average coefficients. The default Value is \code{ma.par="b"}.
}
  \item{Counting.Process}{a character-string that is the label of the Counting process. Defaults to \code{"N"}.
}
  \item{Intensity.var}{a character-string that is the label of the Intensity process. Defaults to \code{"lambda"}
}
  \item{Latent.var}{
a character-string that is the label of the unobserved process. Defaults to \code{"x"}.
}
  \item{time.var}{
the name of the time variable.
}
  \item{Type.Jump}{
a logical value. If it is \code{TRUE}, the jump size is deterministic. 
}
  \item{XinExpr}{ 
a logical variable. The default value \code{XinExpr=FALSE} implies that the starting condition for \code{Latent.var} is zero. If \code{XinExpr=TRUE}, each component of \code{Latent.var} has a parameter as a initial value.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Model an object of \code{\link{yuima.carmaHawkes-class}}.
}
\references{
Mercuri, L., Perchiazzo, A., & Rroji, E. (2022). A Hawkes model with CARMA (p, q) intensity. \emph{Insurance: Mathematics and Economics}, \bold{116}, 1-26. \doi{10.1016/j.insmatheco.2024.01.007}.
}
\author{
The YUIMA Project Team

Contacts: Lorenzo Mercuri \email{lorenzo.mercuri@unimi.it}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
# Definition of an Hawkes with exponential Kernel
mod1 <- setCarmaHawkes(p = 1, q = 0)

# Definition of an Hawkes with a CARMA(2,1) Intensity process
mod2 <- setCarmaHawkes(p = 2, q = 1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
