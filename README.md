# yuima: Simulation and Inference for SDEs and Other Stochastic Processes

[![](https://cranlogs.r-pkg.org/badges/grand-total/yuima?color=blue)](https://cran.r-project.org/package=yuima) [![](https://img.shields.io/badge/doi-10.18637%2Fjss.v057.i04-orange)](https://doi.org/10.18637/jss.v057.i04)

This R package performs various central statistical analyses such as quasi maximum likelihood estimation, adaptive Bayes estimation, structural change point analysis, hypotheses testing, asynchronous covariance estimation, lead-lag estimation, LASSO model selection, and so on. It also supports stochastic numerical analysis by fast computation of the expected value of functionals of stochastic processes through automatic asymptotic expansion by means of the Malliavin calculus. All models can be multidimensional, multiparametric or non parametric.

## Installation

Install the stable release from CRAN:

```R
install.packages("yuima")
```

## For developers

When you add Rcpp or other C code, please always do, within R-devel Console:

```R
tools::package_native_routine_registration_skeleton('yuima',,,FALSE)
```

take the output and update the file [`src/yuima_init.c`](https://github.com/yuimaproject/yuima/blob/main/src/yuima_init.c) with the above output.

### Help

See the help of "package_native_routine_registration_skeleton"

```R
library(tools)
?package_native_routine_registration_skeleton
```

[This page](https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols) is also of interest.

## Acknowledgments

The Project has been funded up to 2010 by the Japan Science Technology (JST) Basic Research Programs PRESTO, Grants-in-Aid for Scientific Research No. 19340021. Presently, the YUIMA Project is supported by the Japan Science and Technology Agency CREST.

