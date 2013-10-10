\name{mmfrac}
\alias{mmfrac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mmfrac
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mmfrac(yuima, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yuima}{
%%     ~~Describe \code{yuima} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (yuima, ...) 
{
    call <- match.call()
    if (missing(yuima)) {
        yuima.stop("yuima object is missing.")
    }
    if (class(yuima) != "yuima") {
        yuima.stop("an object of class yuima is needed.")
    }
    Ddiffx0 <- D(yuima@model@diffusion[[1]], yuima@model@state.variable) == 
        0
    Ddifft0 <- D(yuima@model@diffusion[[1]], yuima@model@time.variable) == 
        0
    bconst <- FALSE
    Hconst <- FALSE
    diffnoparam <- length(yuima@model@parameter@diffusion) == 
        0
    if (Ddiffx0 & Ddifft0 & diffnoparam) {
        if (eval(yuima@model@diffusion[[1]]) == 0) {
            yuima.stop("the model has no fractional part.")
        }
    }
    if (length(yuima@model@parameter@drift) != 1) {
        yuima.stop("the drift is malformed.")
    }
    dx <- D(yuima@model@drift, yuima@model@state.variable)
    dbx <- D(yuima@model@diffusion[[1]], yuima@model@state.variable) == 
        0
    dbt <- D(yuima@model@diffusion[[1]], yuima@model@time.variable) == 
        0
    bxx <- D(dx, yuima@model@state.variable) == 0
    bmix <- D(dx, yuima@model@time.variable) == 0
    isfOU <- (bxx || bmix) && (dbx && dbt)
    if (!isfOU) {
        yuima.stop("estimation not available for this model")
    }
    process <- yuima@data@zoo.data[[1]]
    T <- yuima@sampling@Terminal
    est <- qgv(yuima, ...)
    H <- est$coeff[1]
    sigma <- est$coeff[2]
    if ((!is.na(H)) & (!is.na(sigma))) {
        theta <- (2 * mean(process^2)/(sigma^2 * gamma(2 * H + 
            1)))^(-1/2/H)
        sH <- sqrt((4 * H - 1) * (1 + (gamma(1 - 4 * H) * gamma(4 * 
            H - 1))/(gamma(2 - 2 * H) * gamma(2 * H))))
        sdtheta <- sqrt(theta) * (sH/(2 * H))/sqrt(T)
    }
    else {
        yuima.warn("Diffusion estimation not available, can not estimate the drift parameter")
    }
    x <- c(est$coeff, theta)
    names(x) <- c(names(est$coeff), yuima@model@parameter@drift)
    sdx <- matrix(, 3, 3)
    diag(sdx) <- c(diag(est$vcov), sdtheta)
    colnames(sdx) <- names(x)
    rownames(sdx) <- names(x)
    obj <- list(coefficients = x, vcov = sdx, call = call)
    class(obj) <- "mmfrac"
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
