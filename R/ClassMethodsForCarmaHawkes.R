# Class 'yuima.carmaHawkes'

setClass("carmaHawkes.info",
         representation(p = "numeric",
                        q = "numeric",
                        Counting.Process = "character",
                        base.Int = "character",
                        ar.par = "character",
                        ma.par = "character",
                        Intensity.var = "character",
                        Latent.var = "character",
                        XinExpr = "logical",
                        Type.Jump = "logical")
)

# Initialization Carma.Hawkes Info
setMethod("initialize", "carmaHawkes.info",
          function(.Object,
                   p=numeric(),
                   q=numeric(),
                   Counting.Process = character(),
                   base.Int=character(),
                   ar.par=character(),
                   ma.par=character(),
                   Intensity.var=character(),
                   Latent.var=character(),
                   XinExpr=logical(),
                   Type.Jump=logical()){
            .Object@p <- p
            .Object@q <- q
            .Object@Counting.Process <- Counting.Process
            .Object@base.Int <- base.Int
            .Object@ar.par <- ar.par
            .Object@ma.par <- ma.par
            .Object@Intensity.var <- Intensity.var
            .Object@Latent.var <- Latent.var
            .Object@XinExpr <- XinExpr
            .Object@Type.Jump <- Type.Jump
            return(.Object)
          })

setClass("yuima.carmaHawkes",
         representation(info="carmaHawkes.info"),
         contains="yuima.model")

setMethod("initialize", "yuima.carmaHawkes",
          function(.Object,
                   info = new("carmaHawkes.info"),
                   drift = expression() ,
                   diffusion = list() ,
                   hurst = 0.5,
                   jump.coeff = expression(),
                   measure=list(),
                   measure.type=character(),
                   parameter = new("model.parameter"),
                   state.variable = "lambda",
                   jump.variable = "N",
                   time.variable = "t",
                   noise.number = numeric(),
                   equation.number = numeric(),
                   dimension = numeric(),
                   solve.variable = character(),
                   xinit = expression(),
                   J.flag = logical()){
            .Object@info <- info
            .Object@drift <- drift
            .Object@diffusion <- diffusion
            .Object@hurst <- hurst
            .Object@jump.coeff <- jump.coeff
            .Object@measure <- measure
            .Object@measure.type <- measure.type
            .Object@parameter <- parameter
            .Object@state.variable <- state.variable
            .Object@jump.variable <- jump.variable
            .Object@time.variable <- time.variable
            .Object@noise.number <- noise.number
            .Object@equation.number <- equation.number
            .Object@dimension <- dimension
            .Object@solve.variable <- solve.variable
            .Object@xinit <- xinit
            .Object@J.flag <- J.flag
            return(.Object)
          })

# setCarmaHawkes function

setCarmaHawkes <- function(p, q, law = NULL, base.Int = "mu0", 
                           ar.par="a", ma.par = "b", 
                           Counting.Process = "N",
                           Intensity.var = "lambda",
                           Latent.var = "x",
                           time.var = "t",
                           Type.Jump = FALSE,
                           XinExpr = FALSE){
  if(is.null(law)){
    mylaw<-setLaw(rng = function(n){rconst(n,1)}, density=function(x){dconst(x,1)},
                  time.var = time.var)
    Type.Jump <- TRUE
  }else{
    mylaw <- law
  }
  ar.par1 <- paste0(ar.par,1:p)
  ma.par1 <- paste0(ma.par,0:q)
  carmaHawkesInfo <- new("carmaHawkes.info",
                         p=p,
                         q=q,
                         Counting.Process = Counting.Process,
                         base.Int=base.Int,
                         ar.par=ar.par1,
                         ma.par=ma.par1,
                         Intensity.var=Intensity.var,
                         Latent.var=Latent.var,
                         XinExpr=XinExpr,
                         Type.Jump=Type.Jump)
  InternalCarma <- setCarma(p,q, loc.par=base.Int, 
                     Carma.var=Intensity.var, 
                     ar.par = ar.par,
                     ma.par = ma.par,
                     Latent.var=Latent.var, diffusion=NULL,
                     # time.variable = time.var, 
                     # jump.variable = Counting.Process,
                     measure= list(df=setLaw()),
                     measure.type="code",
                     XinExpr = XinExpr)
  # prova1 <- setCarma(p,q, loc.par=base.Int,
  #                    Carma.var=Intensity.var,
  #                    ar.par = ar.par,
  #                    ma.par = ma.par,
  #                    Latent.var=Intensity.var, diffusion=NULL,
  #                    time.variable = "s", jump.variable = Counting.Process,
  #                    measure= list(df= mylaw),
  #                    measure.type="code",
  #                    XinExpr = XinExpr)
  InternalCarma@measure<-list(df=mylaw)
  InternalCarma@time.variable<-time.var
  InternalCarma@jump.variable <- Counting.Process
  CARMA_HAWKES <- new("yuima.carmaHawkes",
                    info = carmaHawkesInfo,
                    drift = InternalCarma@drift,
                    diffusion = InternalCarma@diffusion,
                    hurst = InternalCarma@hurst,
                    jump.coeff = InternalCarma@jump.coeff,
                    measure = InternalCarma@measure,
                    measure.type = InternalCarma@measure.type,
                    parameter = InternalCarma@parameter,
                    state.variable = InternalCarma@state.variable,
                    jump.variable = InternalCarma@jump.variable,
                    time.variable = InternalCarma@time.variable,
                    noise.number =  InternalCarma@noise.number,
                    equation.number = InternalCarma@equation.number,
                    dimension = InternalCarma@dimension,
                    solve.variable = InternalCarma@solve.variable,
                    xinit = InternalCarma@xinit,
                    J.flag = InternalCarma@J.flag)
  return(CARMA_HAWKES)
}