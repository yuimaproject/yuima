setClass("info.PPR",
         representation(allparam = "character",
                        allparamPPR = "character",
                        common ="character",
                        counting.var = "character",
                        var.dx = "character",
                        upper.var = "character",
                        lower.var = "character",
                        covariates = "character",
                        var.dt = "character",
                        additional.info = "character",
                        Info.measure = "list",
                        RegressWithCount = "logical",
                        IntensWithCount = "logical")
         )

setClass("yuima.PPR",
         representation(PPR = "info.PPR",
                        gFun = "info.Map",
                        Kernel = "Integral.sde"),
         contains="yuima"
)

setMethod("initialize",
          "info.PPR",
          function(.Object,
                   allparam = character(),
                   allparamPPR = character(),
                   common = character(),
                   counting.var = character(),
                   var.dx = character(),
                   upper.var = character(),
                   lower.var = character(),
                   covariates = character(),
                   var.dt = character(),
                   additional.info = character(),
                   Info.measure = list(),
                   RegressWithCount = FALSE,
                   IntensWithCount = TRUE){
            .Object@allparam <- allparam
            .Object@allparamPPR <- allparamPPR
            .Object@common <- common
            .Object@counting.var <- counting.var
            .Object@var.dx <- var.dx
            .Object@upper.var <- upper.var
            .Object@lower.var <- lower.var
            .Object@covariates <- covariates
            .Object@var.dt <- var.dt
            .Object@additional.info <- additional.info
            .Object@Info.measure <- Info.measure
            .Object@RegressWithCount <- RegressWithCount
            .Object@IntensWithCount <- IntensWithCount
            return(.Object)
          }
)

setMethod("initialize",
          "yuima.PPR",
          function(.Object,
                   PPR = new("info.PPR"),
                   gFun = new("info.Map"),
                   Kernel = new("Integral.sde"),
                   yuima = new("yuima")){
            #.Object@param <- param
            .Object@PPR <- PPR
            .Object@gFun <- gFun
            .Object@Kernel <- Kernel
            .Object@data <- yuima@data
            .Object@model <- yuima@model
            .Object@sampling <- yuima@sampling
            .Object@characteristic <- yuima@characteristic
            .Object@functional <- yuima@functional
            return(.Object)
          }
)

setClass("yuima.Hawkes",
         contains="yuima.PPR"
)

# Class yuima.PPR.qmle

setClass("yuima.PPR.qmle",representation(
  model = "yuima.PPR"),
  contains="mle"
)


setClass("summary.yuima.PPR.qmle",
         representation(
           model = "yuima.PPR"),
         contains="summary.mle"
)

setMethod("show", "summary.yuima.PPR.qmle",
          function (object)
          {
            
            cat("Quasi-Maximum likelihood estimation\n\nCall:\n")
            print(object@call)
            cat("\nCoefficients:\n")
            print(coef(object))
            cat("\n-2 log L:", object@m2logL, "\n")
            
          }
)
