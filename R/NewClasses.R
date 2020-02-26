# Here we insert new classes for extending the object of classes yuima
setClass("param.Map",
         representation(out.var = "character",
                        allparam = "character",
                        allparamMap = "character",
                        common = "character",
                        Input.var = "character",
                        time.var = "character"))

setClass("info.Map",
         representation(formula="vector",
                        dimension="numeric",
                        type="character",
                        param = "param.Map"))


setClass("yuima.Map",
         representation(Output = "info.Map"),
         contains="yuima"
           )

# Initialization

setMethod("initialize",
           "param.Map",
           function(.Object, out.var = character(),
                    allparam = character(),
                    allparamMap = character(),
                    common = character(),
                    Input.var = character(),
                    time.var = character()){
             .Object@out.var <- out.var
             .Object@allparam <- allparam
             .Object@allparamMap <- allparamMap
             .Object@common <- common
             .Object@Input.var <-Input.var
             .Object@time.var <- time.var
             return(.Object)
           }
)
#
setMethod("initialize",
          "info.Map", function(.Object,
                                  formula = vector(mode = expression),
                                  dimension = numeric(),
                                  type = character(),
                                  param = new("param.Map")){
                            .Object@formula <- formula
                            .Object@dimension <- dimension
                            .Object@type <- type
                            .Object@param <- param
                            return(.Object)
                          }
          )

setMethod("initialize",
          "yuima.Map",
          function(.Object,
                   #param = new("param.Map"),
                   Output = new("info.Map"),
                   yuima = new("yuima")){
            #.Object@param <- param
            .Object@Output <- Output
            .Object@data <- yuima@data
            .Object@model <- yuima@model
            .Object@sampling <- yuima@sampling
            .Object@characteristic <- yuima@characteristic
            .Object@functional <- yuima@functional
            return(.Object)

          }
)
#
# Class for yuima.integral  is structured as follows:

#   param.Integral
#     Integral$param$allparam
#     Integral$param$common
#     Integral$param$IntegrandParam

setClass("param.Integral",representation(allparam = "character",
  common = "character", Integrandparam = "character")
)
#
setMethod("initialize","param.Integral",
          function(.Object, allparam = character(),
                   common = character(),
                   Integrandparam = character()){
            .Object@allparam <- allparam
            .Object@common <- common
            .Object@Integrandparam <- Integrandparam
            return(.Object)
          }
)
#
# #   variable.Integral
# #     Integral$var.dx
# #     Integral$lower.var
# #     Integral$upper.var
# #     Integral$out.var
# #     Integral$var.time <-"s"
#
setClass("variable.Integral",
         representation(var.dx = "character",
                        lower.var = "character",
                        upper.var = "character",
                        out.var = "character",
                        var.time = "character")
)
#
setMethod("initialize","variable.Integral",
          function(.Object,
                   var.dx = character(),
                   lower.var = character(),
                   upper.var = character(),
                   out.var = character(),
                   var.time = character()){
            .Object@var.dx <- var.dx
            .Object@lower.var <- lower.var
            .Object@upper.var <- upper.var
            .Object@out.var <- out.var
            .Object@var.time <- var.time
            return(.Object)
          }
          )

#   Integrand
#     Integral$IntegrandList
#     Integral$dimIntegrand

setClass("Integrand",
         representation(IntegrandList = "list",
                        dimIntegrand = "numeric")
         )
setMethod("initialize","Integrand",
          function(.Object,
                   IntegrandList = list(),
                   dimIntegrand = numeric()){
            .Object@IntegrandList <- IntegrandList
            .Object@dimIntegrand <- dimIntegrand
            return(.Object)
          }
          )
#
# #   Integral.sde
#
setClass("Integral.sde", representation(param.Integral = "param.Integral",
                                        variable.Integral = "variable.Integral", Integrand = "Integrand")
)
#
setMethod("initialize", "Integral.sde",
          function(.Object,
                   param.Integral = new("param.Integral"),
                   variable.Integral = new("variable.Integral"),
                   Integrand = new("Integrand")){
            .Object@param.Integral <- param.Integral
            .Object@variable.Integral <- variable.Integral
            .Object@Integrand <- Integrand
            return(.Object)
          }
)
#
# # yuima.Integral
#
setClass("yuima.Integral", representation(
  Integral = "Integral.sde"),
  contains = "yuima"
)
#
setMethod("initialize", "yuima.Integral",
          function(.Object,
                   Integral = new("Integral.sde"),
                   yuima = new("yuima")){
            .Object@Integral <- Integral
            #.Object@param <- param
            #.Object@Output <- Output
            .Object@data <- yuima@data
            .Object@model <- yuima@model
            .Object@sampling <- yuima@sampling
            .Object@characteristic <- yuima@characteristic
            .Object@functional <- yuima@functional
            return(.Object)
          }
)
#


# yuima.multimodel. We replacate the yuima.model class in order to
# describe from mathematical point of view the multi dimensional jump
# diffusion model
setClass("yuima.multimodel",
         contains="yuima.model")

setClass("yuima.snr", representation(call = "call", coef = "numeric", snr = "numeric", model = "yuima.model"), prototype = list(call = NULL, coef = NULL, snr = NULL, model = NULL))

## yuima.qmle.incr-class
#setClassUnion("yuima.qmle.incr", members=c("yuima.carma.qmle","cogarch.est.incr"))
setClass("yuima.qmleLevy.incr",representation(Incr.Lev = "ANY",
                                              logL.Incr = "ANY",
                                              minusloglLevy="function",
                                              Levydetails= "list",
                                              Data = "ANY"),
         contains="yuima.qmle")

setMethod("initialize", "yuima.qmleLevy.incr",
          function(.Object,
                   Incr.Lev = NULL,
                   logL.Incr = NULL,
                   minusloglLevy=function(){NULL},
                   Levydetails= list(),
                   Data=NULL,
                   yuima = new("yuima.qmle")){
            .Object@Incr.Lev <- Incr.Lev
            #.Object@param <- param
            #.Object@Output <- Output
            .Object@logL.Incr <- logL.Incr
            .Object@Levydetails<- Levydetails
            .Object@minusloglLevy <- minusloglLevy 
            .Object@Data <- Data
            .Object@model <- yuima@model
            .Object@call <- yuima@call
            .Object@coef <- yuima@coef
            .Object@fullcoef <- yuima@fullcoef
            .Object@vcov <- yuima@vcov
            .Object@min <-yuima@min
            .Object@details<- yuima@details
            .Object@minuslogl<-yuima@minuslogl
            .Object@nobs<-yuima@nobs
            .Object@method<-yuima@method
            return(.Object)
          }
)

