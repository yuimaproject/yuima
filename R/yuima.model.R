setMethod("initialize", "model.parameter",
          function(.Object,
                   all,
                   common,
                   diffusion,
                   drift,
                   jump,
                   measure){
            .Object@all <- all
            .Object@common <- common
            .Object@diffusion <- diffusion
            .Object@drift <- drift
            .Object@jump <- jump
            .Object@measure <- measure
            return(.Object)
          })

setMethod("initialize", "yuima.model",
          function(.Object,
                   drift,
                   diffusion,
				   hurst,
                   jump.coeff,
                   measure,
                   measure.type,
                   parameter,
                   state.variable,
                   jump.variable,
                   time.variable,
                   noise.number,
                   equation.number,
                   dimension,
                   solve.variable,
                   xinit,
                   J.flag){
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

## setModel
## setter of class 'yuima.model'
## set yuima model from SDE
setModel <- function(drift=NULL,
                     diffusion=NULL,
					 hurst=0.5,
                     jump.coeff=character(),
                     measure=list(),
                     measure.type=character(),
                     state.variable="x",
                     jump.variable="z",
                     time.variable="t",
                     solve.variable,
                     xinit=NULL){
  ##::measure and jump term #####################################
  
  ##::initialize objects ########
  MEASURE <- list()
  ##::end initialize objects ########
  
  ##::make type function list ####
  CPlist <- c("dnorm", "dgamma", "dexp")
  codelist <- c("rIG", "rNIG", "rgamma", "rbgamma", "rngamma", "rstable")
  ##::end make type function list ####
  
  if(!length(measure.type)){
    if( length(jump.coeff) || length(measure) ){
      cat("\nmeasure type isn't matched with jump term.\n")
      return(NULL)
    }
    jump.variable <- character()
    measure.par <- character()
  }else{
    if( !length(jump.coeff) || !length(measure) ){
      cat("\nmeasure type isn't matched with jump term.\n")
      return(NULL)
    }else if(length(jump.coeff)!=1){
      cat("\nmulti dimentional jump term is not supported yet.\n")
      return(NULL)
    }
    
    if(measure.type=="CP"){ ##::CP
      if(length(measure)!=2){
        cat(paste("\nlength of measure must be two on type", measure.type, ".\n"))
        return(NULL)
      }
      
      if(!is.list(measure)){          
        measure <- list(intensity=measure[1], df=measure[2])
      }else{
        if(length(measure[[1]])!=1 || length(measure[[2]])!=1){
          cat("\nmulti dimentional jump term is not supported yet.\n")
          return(NULL)
        }
        ##::naming measure list ########
        tmpc <- names(measure)
        if(is.null(tmpc)){
          names(measure) <- c("intensity", "df")
        }else{
          whichint <- match("intensity", tmpc)            
          whichdf <- match("df", tmpc)
          if(!is.na(whichint)){
            if(names(measure)[-whichint]=="" || names(measure)[-whichint]=="df"){
              names(measure)[-whichint] <- "df"
            }else{
              cat("\nnames of measure are incorrect.\n")
              return(NULL)
            }
          }else if(!is.na(whichdf)){
            if(names(measure)[-whichdf]=="" || names(measure)[-whichdf]=="intensity"){
              names(measure)[-whichdf] <- "intensity"
            }else{
              cat("\nnames of measure are incorrect.\n")
              return(NULL)
            }
          }else{
            cat("\nnames of measure are incorrect.\n")
            return(NULL)
          }
        }
        ##::end naming measure list ########
      }
      
      ##::check df name ####################
      tmp <- regexpr("\\(", measure$df)[1]
      measurefunc <- substring(measure$df, 1, tmp-1)
      if(!is.na(match(measurefunc, codelist))){
        cat(paste("\ndistribution function", measurefunc, "should be defined as type code.\n"))
        return(NULL)
      }else if(is.na(match(measurefunc, CPlist))){
        warning(paste("\ndistribution function", measurefunc, "is not officialy supported as type CP.\n"))
      }
      MEASURE$df$func <- eval(parse(text=measurefunc))
      MEASURE$df$expr <- parse(text=measure$df)
      MEASURE$intensity <- parse(text=measure$intensity)
      
      measure.par <- unique( c( all.vars(MEASURE$intensity), all.vars(MEASURE$df$expr) ) ) 
      ##measure.par$intensity <- unique(all.vars(MEASURE$intensity))
      ##::end check df name ####################
      ##::end CP
    }else if(measure.type=="code"){ ##::code
      if(length(measure)!=1){
        cat(paste("\nlength of measure must be one on type", measure.type, ".\n"))
        return(NULL)
      }
      
      if(!is.list(measure)){
        measure <- list(df=measure)
      }else{
        if(length(measure[[1]])!=1){
          cat("\nmulti dimentional jump term is not supported yet.\n")
          return(NULL)
        }
        ##::naming measure list #############
        if(is.null(names(measure)) || names(measure)=="df"){
          names(measure) <- "df"
        }else{
          cat("\nname of measure is incorrect.\n")
          return(NULL)
        }
        ##::end naming measure list #############
      }
      
      ##::check df name ####################
      tmp <- regexpr("\\(", measure$df)[1]
      measurefunc <- substring(measure$df, 1, tmp-1)
      if(!is.na(match(measurefunc, CPlist))){
        cat(paste("\ndistribution function", measurefunc, "should be defined as type CP.\n"))
        return(NULL)
      }else if(is.na(match(measurefunc, codelist))){
        warning(paste("\ndistribution function", measurefunc, "is not officialy supported as type code.\n"))
      }
      ##MEASURE$df$func <- eval(parse(text=measurefunc))
      MEASURE$df$expr <- parse(text=measure$df)
      
      measure.par <- unique(all.vars(MEASURE$df$expr))
      ##::end check df name ####################
      ##::end code
    }else if(measure.type=="density"){ ##::density
      if(length(measure)!=1){
        cat(paste("\nlength of measure must be one on type", measure.type, ".\n"))
        return(NULL)
      }
      
      if(!is.list(measure)){
        measure <- list(df=measure)
      }else{
        if(length(measure[[1]])!=1){
          cat("\nmulti dimentional jump term is not supported yet.\n")
          return(NULL)
        }
        
        ##::naming measure list #############
        if(is.null(names(measure))){
          names(measure) <- "df"
        }else if(names(measure)!="density" && names(measure)!="df"){
          cat("\nname of measure is incorrect.\n")
          return(NULL)
        }
        ##::end naming measure list #############
      }
      
      ##::check df name ####################
      tmp <- regexpr("\\(", measure[[names(measure)]])[1]
      measurefunc <- substring(measure[[names(measure)]], 1, tmp-1)
      if(!is.na(match(measurefunc, CPlist))){
        cat(paste("\ndistribution function", measurefunc, "should be defined as type CP.\n"))
        return(NULL)
      }else if(!is.na(match(measurefunc, codelist))){
        cat(paste("\ndistribution function", measurefunc, "should be defined as type code.\n"))
        return(NULL)
      }
      MEASURE[[names(measure)]]$func <- eval(parse(text=measurefunc))
      MEASURE[[names(measure)]]$expr <- parse(text=measure[[names(measure)]])
      
      measure.par <- unique(all.vars(MEASURE[[names(measure)]]$expr))
      ##::end check df name ####################
      ##::end density
    }else{ ##::else
      cat(paste("\nmeasure type", measure.type, "isn't supported.\n"))
      return(NULL)
    }
    n.eqn3 <- 1
    n.jump <- 1
  }
  
  ##::end measure and jump term #####################################
  
  
  ##:: check for errors and reform values
  if(any(time.variable %in% state.variable)){
    cat("\ntime and state(s) variable must be different\n")
    return(NULL)
  }
  if(is.null(dim(drift))){ # this is a vector
    n.eqn1 <- length(drift)
    n.drf <- 1
  }else{ # it is a matrix
    n.eqn1 <- dim(drift)[1]
    n.drf <- dim(drift)[2]
  }
  
  if(is.null(dim(diffusion))){ # this is a vector
    n.eqn2 <- length(diffusion)
    n.noise <- 1
  }else{ # it is a matrix
    n.eqn2 <- dim(diffusion)[1]
    n.noise <- dim(diffusion)[2]
  }
  
  if(is.null(diffusion)){
    diffusion <- rep("0", n.eqn1)
    n.eqn2 <- n.eqn1
    n.noise <- 1
  }

  ## TBC
  n.eqn3 <- n.eqn1
  
  if(!length(measure)){
    n.eqn3 <- n.eqn1
  }

  if(n.eqn1 != n.eqn2 || n.eqn1 != n.eqn3){
    cat("\nMalformed model, number of equations in the drift and diffusion do not match\n")
    return(NULL)
  }
  n.eqn <- n.eqn1
  
  if(is.null(xinit)){
    xinit <- numeric(n.eqn)
  }else if(length(xinit) != n.eqn){
    if(length(xinit)==1){
      xinit <- rep(xinit, n.eqn)
    }else{
      cat("\nDimension of xinit variables missmuch.\n")
      return(NULL)
    }
  }
  
  if(missing(solve.variable)){
    cat("\nSolution variable (lhs) not specified. Trying to use state variables\n")
    solve.variable <- state.variable
  }
  if(n.eqn != length(solve.variable)){
    cat("\nMalformed model, number of solution variables (lhs) do no match number of equations (rhs)\n")
    return(NULL)
  }
  
  loc.drift <- matrix(drift, n.eqn, n.drf)
  loc.diffusion <- matrix(diffusion, n.eqn, n.noise)
  
  ##:: allocate vectors
  DRIFT <- vector(n.eqn, mode="expression")
  DIFFUSION <- vector(n.eqn, mode="list")
  JUMP <- vector(n.eqn, mode="expression")
  
  ##:: function to make expression from drift characters
  pre.proc <- function(x){
    for(i in 1:length(x)){
      if(length(parse(text=x[i]))==0){
        x[i] <- "0"
      }
    }
    parse(text=paste(sprintf("(%s)", x), collapse="+"))
  }
  
  ##:: make expressions of drifts and diffusions
  for(i in 1:n.eqn){
    DRIFT[i] <- pre.proc(loc.drift[i,])
    for(j in 1:n.noise){
      expr <- parse(text=loc.diffusion[i,j])
      if(length(expr)==0){
        expr <- expression(0)  # expr must have something
      }
      DIFFUSION[[i]][j] <- expr
    }
  }
  JUMP <- parse(text=jump.coeff)
  
  ##:: get parameters in drift expression
  drift.par <- unique(all.vars(DRIFT))
  drift.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), drift.par)))
  if(length(drift.idx)>0){
    drift.par <- drift.par[-drift.idx]
  }
  
  ##:: get parameters in diffusion expression
  diff.par <- unique(unlist(lapply(DIFFUSION, all.vars)))
  diff.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), diff.par)))
  if(length(diff.idx)>0){
    diff.par <- diff.par[-diff.idx]
  }

  ##:: get parameters in jump expression
  J.flag <- FALSE
  jump.par <- unique(all.vars(JUMP))
  if(length(na.omit(match(jump.par, jump.variable)))){
    J.flag <- TRUE
  }
  jump.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), jump.par)))
  if(length(jump.idx)>0){
    jump.par <- jump.par[-jump.idx]
  }

  ##:: get parameters in measure expression
  measure.idx <- as.numeric(na.omit(match(c(state.variable, time.variable, jump.variable, solve.variable), measure.par)))
  if(length(measure.idx)>0){
    measure.par <- measure.par[-measure.idx]
  }
  
  ##:: order parameters for 'yuima.pars'
  ##id1 <- which(diff.par %in% drift.par)
  ##id2 <- which(drift.par %in% diff.par)
  ##common <- unique(c(diff.par[id1], drift.par[id2]))
  common <- c(drift.par, diff.par)
  common <- common[duplicated(common)]
  if(length(measure)){
    common <- c(common, jump.par)
    common <- common[duplicated(common)]
    common <- c(common, measure.par)
    common <- common[duplicated(common)]
  }
  all.par <- unique(c(drift.par, diff.par, jump.par, measure.par))
    
  ##:: instanciate class
  tmppar <- new("model.parameter",
                all= all.par,
                common= common,
                diffusion= diff.par,
                drift= drift.par,
                jump= jump.par,
                measure= measure.par)
  tmp <- new("yuima.model",
             drift= DRIFT,
             diffusion= DIFFUSION,
			 hurst=as.numeric(hurst),
             jump.coeff=JUMP,
             measure= MEASURE,
             measure.type= measure.type,
             parameter= tmppar, 
             state.variable= state.variable,
             jump.variable= jump.variable,
             time.variable= time.variable,
             noise.number= n.noise,
             equation.number= n.eqn,
             dimension= c(
               length(tmppar@all),
               length(tmppar@common),
               length(tmppar@diffusion),
               length(tmppar@drift),
               length(tmppar@jump),
               length(tmppar@measure)
               ),
             solve.variable= solve.variable,
             xinit= xinit,
             J.flag <- J.flag)
  return(tmp)
}
