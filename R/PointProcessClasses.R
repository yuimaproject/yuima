setClass("info.Ppr",
         representation(allparam = "character",
                        allparamPpr = "character",
                        common ="character",
                        counting.var = "character",
                        var.dx = "character",
                        upper.var = "character",
                        lower.var = "character",
                        covariates = "character",
                        var.dt = "character",
                        additional.info = "character",
                        Info.measure = "list")
         )

setClass("yuima.Ppr",
         representation(Ppr = "info.Ppr",
                        gFun = "info.Output",
                        Kernel = "Integral.sde"),
         contains="yuima"
)

setMethod("initialize",
          "info.Ppr",
          function(.Object,
                   allparam = character(),
                   allparamPpr = character(),
                   common = character(),
                   counting.var = character(),
                   var.dx = character(),
                   upper.var = character(),
                   lower.var = character(),
                   covariates = character(),
                   var.dt = character(),
                   additional.info = character(),
                   Info.measure = list()){
            .Object@allparam <- allparam
            .Object@allparamPpr <- allparamPpr
            .Object@common <- common
            .Object@counting.var <- counting.var
            .Object@var.dx <- var.dx
            .Object@upper.var <- upper.var
            .Object@lower.var <- lower.var
            .Object@covariates <- covariates
            .Object@var.dt <- var.dt
            .Object@additional.info <- additional.info
            .Object@Info.measure <- Info.measure
            return(.Object)
          }
)

setMethod("initialize",
          "yuima.Ppr",
          function(.Object,
                   Ppr = new("info.Ppr"),
                   gFun = new("info.Output"),
                   Kernel = new("Integral.sde"),
                   yuima = new("yuima")){
            #.Object@param <- param
            .Object@Ppr <- Ppr
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
         contains="yuima.Ppr"
)
