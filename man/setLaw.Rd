\name{setLaw}
\alias{setLaw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random variable constructor }
\description{Constructor of a random variable}
\usage{
setLaw(rng = function(n, ...) {
    NULL
}, density = function(x, ...) {
    NULL
}, cdf = function(q, ...) {
    NULL
}, quant = function(p, ...) {
    NULL
}, characteristic = function(u, ...) {
    NULL
}, time.var = "t", dim = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rng}{
function}
  \item{density}{
function}
  \item{cdf}{
function}
  \item{characteristic}{
function}
  \item{quant}{
function}
  \item{time.var}{
label}
  \item{dim}{
label}
}
\details{
Insert additional info}
\value{
object of class yuima.law}
\author{
YUIMA TEAM}
\note{
Insert additional info}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\examples{}
