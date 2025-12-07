# Class Definitions
# This source MUST be loaded first

# Class 'yuima.pars'

# parameter object included in 'yuima.model'
setClass("model.parameter", representation(
  all = "character",
  common = "character",
  diffusion = "character",
  drift = "character",
  jump = "character",
  measure = "character",
  # Insert parameters for starting conditions
  xinit = "character"
))

# Class 'yuima.model'
setClass("yuima.model", representation(
  drift = "expression",
  diffusion = "list",
  hurst = "ANY",
  jump.coeff = "list",
  # jump.coeff="expression",
  measure = "list",
  measure.type = "character",
  parameter = "model.parameter",
  state.variable = "character",
  jump.variable = "character",
  time.variable = "character",
  noise.number = "numeric",
  equation.number = "numeric",
  dimension = "numeric",
  solve.variable = "character",
  #                                       xinit="numeric",
  xinit = "expression",
  J.flag = "logical"
))

# Class 'carma.info'
setClass(
  "carma.info",
  representation(
    p = "numeric",
    q = "numeric",
    loc.par = "character",
    scale.par = "character",
    ar.par = "character",
    ma.par = "character",
    lin.par = "character",
    Carma.var = "character",
    Latent.var = "character",
    XinExpr = "logical"
  )
)

# Class 'yuima.carma'

setClass("yuima.carma",
  representation(info = "carma.info"),
  contains = "yuima.model"
)

# Class Compound Poisson
setClass("yuima.poisson", contains = "yuima.model")


# Class 'yuima.data'

# we want yuimaS4 to use any class of data as input
# the original data will be stored in OrigData
# we convert these objects internally to "zoo" object
# in the future, we may want to use more flexible
# classes

setClass("yuima.data", representation(
  original.data = "ANY",
  zoo.data = "ANY"
))


# Class 'yuima.sampling'

# sampling is now empty, but should give informations on the sampling
# type, rate, deltas, etc.

setClass("yuima.sampling", representation(
  Initial = "numeric",
  Terminal = "numeric",
  n = "numeric",
  delta = "numeric",
  grid = "ANY",
  random = "ANY",
  regular = "logical",
  sdelta = "numeric",
  sgrid = "ANY",
  oindex = "ANY",
  interpolation = "character"
))

# Class 'yuima.functional'

# functional model used in 'asymptotic term' procedure

setClass("yuima.functional", representation(
  F = "ANY",
  f = "list",
  xinit = "numeric",
  e = "numeric"
))


# Class 'yuima'

# this is the principal class of yuima project. It may contain up to
# three slots for now: the data, the model and the sampling

setClass("yuima.characteristic", representation(
  equation.number = "numeric",
  time.scale = "numeric"
))


setClass("yuima", representation(
  data = "yuima.data",
  model = "yuima.model",
  sampling = "yuima.sampling",
  characteristic = "yuima.characteristic",
  functional = "yuima.functional"
))

# Class yuima.carma.qmle
setClass("yuima.carma.qmle", representation(
  Incr.Lev = "ANY",
  model = "yuima.carma",
  logL.Incr = "ANY"
),
contains = "mle"
)



setClass("yuima.qmle", representation(
  model = "yuima.model"
),
contains = "mle"
)

setClass("yuima.CP.qmle", representation(
  Jump.times = "ANY",
  Jump.values = "ANY",
  X.values = "ANY",
  model = "yuima.model",
  threshold = "ANY"
),
contains = "mle"
)

setClass("summary.yuima.carma.qmle", representation(
  MeanI = "ANY",
  SdI = "ANY",
  logLI = "ANY",
  TypeI = "ANY",
  NumbI = "ANY",
  StatI = "ANY",
  model = "yuima.carma",
  Additional.Info = "ANY"
),
contains = "summary.mle"
)

setClass("summary.yuima.CP.qmle",
  representation(
    NJ = "ANY",
    MeanJ = "ANY",
    SdJ = "ANY",
    MeanT = "ANY",
    Jump.times = "ANY",
    Jump.values = "ANY",
    X.values = "ANY",
    model = "yuima.model",
    threshold = "ANY"
  ),
  contains = "summary.mle"
)


setClass("summary.yuima.qmle",
  representation(
    model = "yuima.model",
    threshold = "ANY",
    Additional.Info = "ANY"
  ),
  contains = "summary.mle"
)



# The yuima.carma.qmle extends the S4 class "mle". It contains three slots: Estimated Levy,
# The description of the carma model and the mle.

setClass("yuima.state_space_model",
  representation(is.observed = "logical"),
  contains = "yuima.model"
)

setClass("yuima.linear_state_space_model",
  representation(
    drift_slope = "list", # list of expressions
    drift_intercept = "list" # list of expressions
  ),
  contains = "yuima.state_space_model"
)

setClass("yuima.linear_state_space_qmle",
  slots = c(
    model = "yuima.linear_state_space_model",
    drop_terms = "numeric",
    explicit = "logical",
    mean_init = "numeric"
  ),
  contains = "yuima.qmle"
)


# kalmanBucyFilter related
setClass(
  "yuima.kalmanBucyFilter",
  representation(
    model = "yuima.linear_state_space_model",
    mean = "ts",
    vcov = "array",
    mean.init = "numeric",
    vcov.init = "matrix",
    delta = "numeric",
    data = "yuima.data"
  )
)

setMethod(
  "show", "yuima.kalmanBucyFilter",
  function(object) {
    start <- start(object@mean)[1]
    end <- end(object@mean)[1]
    freq <- frequency(object@mean)
    time_points <- seq(start, end, by = 1/freq)
    cat("Kalman-Bucy Filter\n")
    if (dim(object@mean)[2] == 1) {
      cat("Mean and variance values:\n")
      var_name <- colnames(object@mean)
      mean_variance_mat <- cbind(as.matrix(object@mean), as.matrix(object@vcov))
      colnames(mean_variance_mat) <- paste(c("Mean of", "Variance of"), var_name)
      rownames(mean_variance_mat) <- time_points
      n <- dim(mean_variance_mat)[1]
      if (n <= 12) {
        print(mean_variance_mat)
      } else {
        print(rbind(head(mean_variance_mat), "...", tail(mean_variance_mat)), quote = FALSE)
      }
    } else {
      cat("Mean values:\n")
      n <- dim(object@mean)[1]
      if (n <= 12) {
        print(object@mean)
      } else {
        mean_mat <- as.matrix(object@mean)
        rownames(mean_mat) <- time_points
        print(rbind(head(mean_mat), "...", tail(mean_mat)), quote = FALSE)
      }
      cat("Variance-covariance matrices\n")
      if (n < 12) {
        print(object@vcov)
      } else {
        for (i in 1:6) {
          cat("\n, , ", i, "\n", sep = "") 
          print(object@vcov[, , i], quote = TRUE) 
        }
        cat("...")
        for (i in (dim(object@vcov)[3]-5):dim(object@vcov)[3]) {
          cat("\n, , ", i, "\n", sep = "") 
          print(object@vcov[, , i], quote = TRUE) 
        }
      }
    }
  }
)

setMethod("summary", "yuima.kalmanBucyFilter", function(object) {
    cat("Summary of estimation by Kalman-Bucy Filter\n")
    cat("Model:\n")
    cat("  A linear state space model.\n")
    cat("  State Variables:", object@model@state.variable[!object@model@is.observed], "\n")
    cat("Mean:\n")
    cat("  A ts object of", dim(object@mean)[2], "variables.\n")
    cat("  Start:      ", start(object@mean)[1], "\n")
    cat("  End:        ", end(object@mean)[1], "\n")
    cat("  Frequency:  ", frequency(object@mean), "\n")
    cat("  Time points:", dim(object@mean)[1], "\n")
    cat("Variance-Covariance Matrix:\n")
    cat(paste0("  A 3D array of (", dim(object@vcov)[1], ", ", dim(object@vcov)[2], ", ", dim(object@vcov)[3], ").\n"))
  }
)

# adaBayes related
setClass(
  "adabayes",
  #contains = "mle",
  slots = c(
    mcmc = "list",
    accept_rate = "list",
    coef = "numeric",
    call = "call",
    vcov = "matrix",
    fullcoef = "numeric",
    fixed = "numeric"
  )
)