\name{yuima.adabayes-class}
\docType{class}
\alias{yuima.adabayes-class}

\title{Class for Adaptive Bayes Estimation of Stochastic Differential Equations}

\description{
    The \code{yuima.adabayes} class is used to describe the output of the functions \code{\link{adaBayes}}.
}
\section{Slots}{
    \describe{
        \item{\code{model}}{is an object of class \code{yuima.model-class}.}
        \item{\code{mcmc}:}{is a list of MCMC objects for all estimated parameters.}
        \item{\code{accept_rate}:}{is a list acceptance rates for diffusion and drift parts.}
        \item{\code{call}}{is an object of class \code{language}.}
        \item{\code{vcov}}{is an object of class \code{matrix}.}
        \item{\code{fullcoef}}{is an object of class \code{numeric} that contains estimated and fixed parameters.} 
        \item{\code{coef}}{is an object of class \code{numeric} that contains estimated parameters.}
        \item{\code{fixed}}{is an object of class \code{numeric} that contains fixed parameters.}
    }
}
