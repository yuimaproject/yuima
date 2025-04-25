## information criteria

IC <- function(drif = NULL, diff = NULL, jump.coeff = NULL, data = NULL, Terminal = 1, add.settings = list(), start, lower, upper, ergodic = TRUE, stepwise = FALSE, weight = FALSE, rcpp = FALSE, ...){ 

  Levy <- FALSE
  if(length(jump.coeff) > 0){
    torf.jump <- (jump.coeff == "0")
    for(i in 1:length(jump.coeff)){
      if(torf.jump[i] == FALSE){
        Levy <- TRUE
        stepwise <- TRUE
        break
      }
    }
  }
  
  if(Levy == FALSE && ergodic == FALSE){
    stepwise <- FALSE
  }
  
  settings <- list(hurst = 0.5, measure = list(), measure.type = character(), state.variable = "x", jump.variable = "z", time.variable = "t", solve.variable = "x")
  if(length(add.settings) > 0){
    match.settings <- match(names(add.settings), names(settings))
    for(i in 1:length(match.settings)){
      settings[[match.settings[i]]] <- add.settings[[i]]
    }
  }
  
  if(stepwise == FALSE){
    # Joint
    ## Candidate models
    yuimas <- NULL
    if(ergodic == TRUE){
      joint <- TRUE
      for(i in 1:length(diff)){
        for(j in 1:length(drif)){
          mod <- setModel(drift = drif[[j]], diffusion = diff[[i]], hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
          if(is.matrix(data) == FALSE){
            n <- length(data)-1
            modsamp <- setSampling(Terminal = Terminal, n = n)
            modyuima <- setYuima(model = mod, sampling = modsamp)
            sub.zoo.data <- list(zoo(x = data, order.by = modyuima@sampling@grid[[1]]))
            names(sub.zoo.data)[1] <- "Series 1"
          }else{
            n <- nrow(data)-1
            modsamp <- setSampling(Terminal = Terminal, n = n)
            modyuima <- setYuima(model = mod, sampling = modsamp)
            sub.zoo.data <- list()
            for(j in 1:ncol(data)){
              sub.zoo.data <- c(sub.zoo.data, list(zoo(x = data[,j], order.by = modyuima@sampling@grid[[1]])))
              names(sub.zoo.data)[j] <- paste("Series", j)
            }
          }
          modyuima@data@zoo.data <- sub.zoo.data
          
          yuimas <- c(yuimas, list(modyuima))
        }
      }
    }else{
      joint <- FALSE
      for(i in 1:length(diff)){
        if(is.matrix(data) == FALSE){
          mod <- setModel(drift = "0", diffusion = diff[[i]], hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
          n <- length(data)-1
          modsamp <- setSampling(Terminal = Terminal, n = n)
          modyuima <- setYuima(model = mod, sampling = modsamp)
          sub.zoo.data <- list(zoo(x = data, order.by = modyuima@sampling@grid[[1]]))
          names(sub.zoo.data)[1] <- "Series 1"
        }else{
          zerovec <- rep("0", length=ncol(data))
          mod <- setModel(drift = zerovec, diffusion = diff[[i]], hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
          n <- nrow(data)-1
          modsamp <- setSampling(Terminal = Terminal, n = n)
          modyuima <- setYuima(model = mod, sampling = modsamp)
          sub.zoo.data <- list()
          for(j in 1:ncol(data)){
            sub.zoo.data <- c(sub.zoo.data, list(zoo(x = data[,j], order.by = modyuima@sampling@grid[[1]])))
            names(sub.zoo.data)[j] <- paste("Series", j)
          }
        }
        modyuima@data@zoo.data <- sub.zoo.data
        
        yuimas <- c(yuimas, list(modyuima))
      }
    }
    mod.num <- length(yuimas)
    
    ## Model comparison
    Esti <- BIC <- QBIC <- AIC <- NULL
    for(i in 1:mod.num){
      yuima <- yuimas[[i]]
      #alpha <- yuima@model@parameter@drift
      #beta <- yuima@model@parameter@diffusion
      
      para.num.init  <- match(yuima@model@parameter@all, names(start))
      para.num.low  <- match(yuima@model@parameter@all, names(lower))
      para.num.upp  <- match(yuima@model@parameter@all, names(upper))
      para.start <- NULL
      para.lower <- NULL
      para.upper <- NULL
      for(j in 1:length(yuima@model@parameter@all)){
        para.start <- c(para.start, list(start[[para.num.init[j]]]))
        para.lower <- c(para.lower, list(lower[[para.num.low[j]]]))
        para.upper <- c(para.upper, list(upper[[para.num.upp[j]]]))
      }
      names(para.start) <- yuima@model@parameter@all
      names(para.lower) <- yuima@model@parameter@all
      names(para.upper) <- yuima@model@parameter@all
      
      mle <- qmle(yuima, start = para.start, lower = para.lower, upper = para.upper, method = "L-BFGS-B", joint = joint, rcpp = rcpp)
      hess <- list(mle@details$hessian)
      hess.diff <- subset(hess[[1]], rownames(hess[[1]])%in%yuima@model@parameter@diffusion, select=yuima@model@parameter@diffusion)
      hess.drif <- subset(hess[[1]], rownames(hess[[1]])%in%yuima@model@parameter@drift, select=yuima@model@parameter@drift)
      
      esti <- list(coef(mle))
      names(esti[[1]]) <- c(yuima@model@parameter@diffusion, yuima@model@parameter@drift)
      bic <- summary(mle)@m2logL+length(yuima@model@parameter@drift)*log(Terminal)+length(yuima@model@parameter@diffusion)*log(n)
      if(det(hess.diff) > 0 && det(hess.drif) > 0){
        qbic <- summary(mle)@m2logL+log(det(hess.diff))+log(det(hess.drif))
      }else{
        qbic <- summary(mle)@m2logL+length(yuima@model@parameter@drift)*log(Terminal)+length(yuima@model@parameter@diffusion)*log(n)
      }
      aic <- summary(mle)@m2logL+2*(length(yuima@model@parameter@drift)+length(yuima@model@parameter@diffusion))
      
      Esti <- c(Esti, esti)
      BIC <- c(BIC, bic)
      QBIC <- c(QBIC, qbic)
      AIC <- c(AIC, aic)
    }
    BIC.opt <- which.min(BIC)
    QBIC.opt <- which.min(QBIC)
    AIC.opt <- which.min(AIC)
  
    ## Names
    if(ergodic == TRUE){
      for(i in 1:length(diff)){
        for(j in 1:length(drif)){
          names(Esti)[(length(drif)*(i-1)+j)] <- paste("scale_", i, " & drift_", j, sep = "") 
        }
      }
      
      BIC <- matrix(BIC, length(drif), length(diff))
      QBIC <- matrix(QBIC, length(drif), length(diff))
      AIC <- matrix(AIC, length(drif), length(diff))
      
      diff.name <- numeric(length(diff))
      drif.name <- numeric(length(drif))
      for(i in 1:length(diff)){
        diff.name[i] <- paste("scale", i, sep = "_") 
      }
      colnames(BIC) <- colnames(QBIC) <- colnames(AIC) <- diff.name
      for(i in 1:length(drif)){
        drif.name[i] <- paste("drift", i, sep = "_")
      }
      rownames(BIC) <- rownames(QBIC) <- rownames(AIC) <- drif.name
    }else{
      for(i in 1:length(diff)){
        names(Esti)[i] <- paste("scale", i, sep = "_") 
      }
      
      diff.name <- numeric(length(diff))
      for(i in 1:length(diff)){
        diff.name[i] <- paste("scale", i, sep = "_") 
      }
      names(BIC) <- names(QBIC) <- diff.name
    }
    
    ## Model weights
    if(weight == TRUE){
      BIC.weight <- exp(-(1/2)*(BIC-BIC[BIC.opt]))/sum(exp(-(1/2)*(BIC-BIC[BIC.opt])))
      QBIC.weight <- exp(-(1/2)*(QBIC-QBIC[QBIC.opt]))/sum(exp(-(1/2)*(QBIC-QBIC[QBIC.opt])))
      AIC.weight <- exp(-(1/2)*(AIC-AIC[AIC.opt]))/sum(exp(-(1/2)*(AIC-AIC[AIC.opt])))
      
      if(ergodic == TRUE){
        BIC.weight <- matrix(BIC.weight, length(drif), length(diff))
        QBIC.weight <- matrix(QBIC.weight, length(drif), length(diff))
        AIC.weight <- matrix(AIC.weight, length(drif), length(diff))
        
        colnames(BIC.weight) <- colnames(QBIC.weight) <- colnames(AIC.weight) <- diff.name
        rownames(BIC.weight) <- rownames(QBIC.weight) <- rownames(AIC.weight) <- drif.name
      }else{
        names(BIC.weight) <- names(QBIC.weight) <- diff.name
      }
    }
    
    ## Results
    diff.copy <- diff
    drif.copy <- drif
    for(i in 1:length(diff)){
      names(diff.copy)[i] <- paste("scale", i, sep = "_") 
    }
    if(ergodic == TRUE){
      for(i in 1:length(drif)){
        names(drif.copy)[i] <- paste("drift", i, sep = "_") 
      }
      diff.BIC.opt <- (BIC.opt-1)%/%length(drif)+1
      diff.QBIC.opt <- (QBIC.opt-1)%/%length(drif)+1
      diff.AIC.opt <- (AIC.opt-1)%/%length(drif)+1
      drif.BIC.opt <- (BIC.opt+(length(drif)-1))%%length(drif)+1
      drif.QBIC.opt <- (QBIC.opt+(length(drif)-1))%%length(drif)+1
      drif.AIC.opt <- (AIC.opt+(length(drif)-1))%%length(drif)+1
    }else{
      drif <- NULL
    }
    
    call <- match.call()
    model.coef <- list(drift = drif.copy, scale = diff.copy)
    if(length(drif) >0){
      bic.selected.coeff <- list(drift = drif[[drif.BIC.opt]], scale = diff[[diff.BIC.opt]])
      qbic.selected.coeff <- list(drift = drif[[drif.QBIC.opt]], scale = diff[[diff.QBIC.opt]])
      aic.selected.coeff <- list(drift = drif[[drif.AIC.opt]], scale = diff[[diff.AIC.opt]])
    }else{
      bic.selected.coeff <- list(drift = NULL, scale = diff[[BIC.opt]])
      qbic.selected.coeff <- list(drift = NULL, scale = diff[[QBIC.opt]])
      aic.selected.coeff <- list(drift = NULL, scale = NULL)
      AIC <- NULL
      AIC.weight <- NULL
    }
    ic.selected <- list(BIC = bic.selected.coeff, QBIC = qbic.selected.coeff, AIC = aic.selected.coeff)
    if(weight == TRUE){
      ak.weight <- list(BIC = BIC.weight, QBIC = QBIC.weight, AIC = AIC.weight)
    }else{
      ak.weight <- NULL
    }
    final_res <- list(call = call, model = model.coef, par = Esti, BIC = BIC, QBIC = QBIC, AIC = AIC, weight = ak.weight, selected = ic.selected)
    
  }else{
    # Stepwise
    pena.aic <- function(yuimaaic, data, pdiff, moment){
      tmp.env <- new.env()
      aic1para <- yuimaaic@model@parameter@diffusion
      for(i in 1:length(aic1para)){
        aic1match <- match(aic1para[i], names(pdiff)[i])
        assign(aic1para[i], pdiff[aic1match], envir=tmp.env)
      }
      aic1state <- yuimaaic@model@state.variable
      aic1ldata <- length(data)-1
      aic1dx <- diff(data)
      assign(aic1state, data, envir=tmp.env)
      
      aic1ter <- yuimaaic@sampling@Terminal
      aic1diff <- eval(yuimaaic@model@diffusion[[1]], envir=tmp.env)
      if(length(aic1diff) == 1){
        aic1diff <- rep(aic1diff, aic1ldata)
      }
      aic1sum <- 0
      for(i in 1:aic1ldata){
        subaic1sum <- (aic1dx[i]/aic1diff[i])^moment
        aic1sum <- aic1sum + subaic1sum
      }
      return(aic1sum/aic1ter)
    }
    
    Esti1 <- BIC1 <- QBIC1 <- AIC1 <- NULL
    Esti2.bic <- Esti2.qbic <- Esti2.aic <- BIC2 <- QBIC2 <- AIC2 <- NULL
    
    # First step
    yuimas1 <- swbeta <- NULL
    if(Levy == TRUE){
      diff <- jump.coeff
    }
    
    for(i in 1:length(diff)){
      ## Candidate models
      if(is.matrix(data) == FALSE){
        mod <- setModel(drift = "0", diffusion = diff[[i]], hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
        n <- length(data)-1
        modsamp <- setSampling(Terminal = Terminal, n = n)
        modyuima <- setYuima(model = mod, sampling = modsamp)
        sub.zoo.data <- list(zoo(x = data, order.by = modyuima@sampling@grid[[1]]))
        names(sub.zoo.data)[1] <- "Series 1"
      }else{
        zerovec <- rep("0", length=ncol(data))
        mod <- setModel(drift = zerovec, diffusion = diff[[i]], hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
        n <- nrow(data)-1
        modsamp <- setSampling(Terminal = Terminal, n = n)
        modyuima <- setYuima(model = mod, sampling = modsamp)
        sub.zoo.data <- list()
        for(j in 1:ncol(data)){
          sub.zoo.data <- c(sub.zoo.data, list(zoo(x = data[,j], order.by = modyuima@sampling@grid[[1]])))
          names(sub.zoo.data)[j] <- paste("Series", j)
        }
      }
      modyuima@data@zoo.data <- sub.zoo.data
      yuimas1 <- c(yuimas1, list(modyuima))
      
      ## Model comparison
      yuima <- modyuima
      swbeta <- c(swbeta, list(yuima@model@parameter@diffusion))
      
      para.num.init  <- match(swbeta[[i]], names(start))
      para.num.low  <- match(swbeta[[i]], names(lower))
      para.num.upp  <- match(swbeta[[i]], names(upper))
      para.start <- NULL
      para.lower <- NULL
      para.upper <- NULL
      for(j in 1:length(swbeta[[i]])){
        para.start <- c(para.start, list(start[[para.num.init[j]]]))
        para.lower <- c(para.lower, list(lower[[para.num.low[j]]]))
        para.upper <- c(para.upper, list(upper[[para.num.upp[j]]]))
      }
      names(para.start) <- swbeta[[i]]
      names(para.lower) <- swbeta[[i]]
      names(para.upper) <- swbeta[[i]]
      
      mle <- qmle(yuima, start = para.start, lower = para.lower, upper = para.upper, method = "L-BFGS-B", joint = FALSE, rcpp = rcpp)
      hess <- mle@details$hessian
      
      esti <- list(coef(mle))
      names(esti[[1]]) <- swbeta[[i]]
      if(is.matrix(data) == FALSE && Levy == TRUE){
        bic <- summary(mle)@m2logL+(length(swbeta[[i]])/(yuima@sampling@delta))*log(yuima@sampling@Terminal)
        qbic <- summary(mle)@m2logL+(length(swbeta[[i]])/(yuima@sampling@delta))*log(yuima@sampling@Terminal)
      }else{
        bic <- summary(mle)@m2logL+length(swbeta[[i]])*log(n)
        if(det(hess) > 0){
          qbic <- summary(mle)@m2logL+log(det(hess))
        }else{
          qbic <- summary(mle)@m2logL+length(swbeta[[i]])*log(n)
        }
      }
      if(is.matrix(data) == FALSE && Levy == TRUE){
        aic <- summary(mle)@m2logL+length(swbeta[[i]])*((1/yuima@sampling@delta)*pena.aic(yuima,data,coef(mle),4)-(pena.aic(yuima,data,coef(mle),2))^2)
      }else{
        aic <- summary(mle)@m2logL+2*length(swbeta[[i]])
      }
      
      Esti1 <- c(Esti1, esti)
      BIC1 <- c(BIC1, bic)
      QBIC1 <- c(QBIC1, qbic)
      AIC1 <- c(AIC1, aic)
      
    }
    BIC.opt1 <- which.min(BIC1)
    QBIC.opt1 <- which.min(QBIC1)
    AIC.opt1 <- which.min(AIC1)
    
    ## Names
    for(i in 1:length(diff)){
      names(Esti1)[i] <- paste("scale", i, sep = "_") 
      names(BIC1)[i] <- paste("scale", i, sep = "_") 
      names(QBIC1)[i] <- paste("scale", i, sep = "_")
      names(AIC1)[i] <- paste("scale", i, sep = "_")
    }
    
    ## Model weights
    if(weight == TRUE){
      BIC.weight1 <- exp(-(1/2)*(BIC1-BIC1[BIC.opt1]))/sum(exp(-(1/2)*(BIC1-BIC1[BIC.opt1])))
      QBIC.weight1 <- exp(-(1/2)*(QBIC1-QBIC1[QBIC.opt1]))/sum(exp(-(1/2)*(QBIC1-QBIC1[QBIC.opt1])))
      AIC.weight1 <- exp(-(1/2)*(AIC1-AIC1[AIC.opt1]))/sum(exp(-(1/2)*(AIC1-AIC1[AIC.opt1])))
      for(i in 1:length(diff)){
        names(BIC.weight1)[i] <- paste("scale", i, sep = "_") 
        names(QBIC.weight1)[i] <- paste("scale", i, sep = "_")
        names(AIC.weight1)[i] <- paste("scale", i, sep = "_") 
      }
    }
    
    # Second step
    ## Use the selection results of first step
    diff.row.bic <- length(yuimas1[[BIC.opt1]]@model@diffusion)
    Diff.esti.bic <- NULL
    Esti1.chr.bic <- as.character(Esti1[[BIC.opt1]])
    Diff.esti.bic <- diff[[BIC.opt1]]
    for(i in 1:diff.row.bic){
      if(length(Esti1.chr.bic) == 1){
        Diff.esti.bic.sub <- gsub(swbeta[[BIC.opt1]][1], Esti1.chr.bic[1], yuimas1[[BIC.opt1]]@model@diffusion[[i]])
      }else{
        Diff.esti.bic.sub <- gsub(swbeta[[BIC.opt1]][1], Esti1.chr.bic[1], yuimas1[[BIC.opt1]]@model@diffusion[[i]])
        for(j in 1:(length(Esti1.chr.bic)-1)){
          Diff.esti.bic.sub <- gsub(swbeta[[BIC.opt1]][(j+1)], Esti1.chr.bic[(j+1)], Diff.esti.bic.sub)
        }
      }
    #  if(class(Diff.esti.bic) == "character"){
      if(inherits(Diff.esti.bic, "character")){
          Diff.esti.bic <- Diff.esti.bic.sub
      } else {
          Diff.esti.bic[i,] <- Diff.esti.bic.sub
      }
    }
    
    diff.row.qbic <- length(yuimas1[[QBIC.opt1]]@model@diffusion)
    Diff.esti.qbic <- NULL
    Esti1.chr.qbic <- as.character(Esti1[[QBIC.opt1]])
    Diff.esti.qbic <- diff[[QBIC.opt1]]
    for(i in 1:diff.row.qbic){
      if(length(Esti1.chr.qbic) == 1){
        Diff.esti.qbic.sub <- gsub(swbeta[[QBIC.opt1]][1], Esti1.chr.qbic[1], yuimas1[[QBIC.opt1]]@model@diffusion[[i]])
      }else{
        Diff.esti.qbic.sub <- gsub(swbeta[[QBIC.opt1]][1], Esti1.chr.qbic[1], yuimas1[[QBIC.opt1]]@model@diffusion[[i]])
        for(j in 1:(length(Esti1.chr.qbic)-1)){
          Diff.esti.qbic.sub <- gsub(swbeta[[QBIC.opt1]][(j+1)], Esti1.chr.qbic[(j+1)], Diff.esti.qbic.sub)
        }
      }
      #if(class(Diff.esti.qbic) == "character"){
      if(inherits(Diff.esti.qbic,"character")){
        Diff.esti.qbic <- Diff.esti.qbic.sub
      } else {
        Diff.esti.qbic[i,] <- Diff.esti.qbic.sub
      }
    }
    
    diff.row.aic <- length(yuimas1[[AIC.opt1]]@model@diffusion)
    Diff.esti.aic <- NULL
    Esti1.chr.aic <- as.character(Esti1[[AIC.opt1]])
    Diff.esti.aic <- diff[[AIC.opt1]]
    for(i in 1:diff.row.aic){
      if(length(Esti1.chr.aic) == 1){
        Diff.esti.aic.sub <- gsub(swbeta[[AIC.opt1]][1], Esti1.chr.aic[1], yuimas1[[AIC.opt1]]@model@diffusion[[i]])
      }else{
        Diff.esti.aic.sub <- gsub(swbeta[[AIC.opt1]][1], Esti1.chr.aic[1], yuimas1[[AIC.opt1]]@model@diffusion[[i]])
        for(j in 1:(length(Esti1.chr.aic)-1)){
          Diff.esti.aic.sub <- gsub(swbeta[[AIC.opt1]][(j+1)], Esti1.chr.aic[(j+1)], Diff.esti.aic.sub)
        }
      }
      #if(class(Diff.esti.aic) == "character"){
      if(inherits(Diff.esti.aic, "character")){
        Diff.esti.aic <- Diff.esti.aic.sub
      } else {
        Diff.esti.aic[i,] <- Diff.esti.aic.sub
      }
    }
    
    yuimas2.bic <- yuimas2.qbic <- yuimas2.aic <- swalpha <- NULL
    for(i in 1:length(drif)){
      ## Candidate models
      if(is.matrix(data) == FALSE){
        mod.bic <- setModel(drift = drif[[i]], diffusion = Diff.esti.bic, hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
        mod.qbic <- setModel(drift = drif[[i]], diffusion = Diff.esti.qbic, hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
        mod.aic <- setModel(drift = drif[[i]], diffusion = Diff.esti.aic, hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
        n <- length(data)-1
        modsamp <- setSampling(Terminal = Terminal, n = n)
        modyuima.bic <- setYuima(model = mod.bic, sampling = modsamp)
        modyuima.qbic <- setYuima(model = mod.qbic, sampling = modsamp)
        modyuima.aic <- setYuima(model = mod.aic, sampling = modsamp)
        sub.zoo.data.bic <- list(zoo(x = data, order.by = modyuima.bic@sampling@grid[[1]]))
        sub.zoo.data.qbic <- list(zoo(x = data, order.by = modyuima.qbic@sampling@grid[[1]]))
        sub.zoo.data.aic <- list(zoo(x = data, order.by = modyuima.aic@sampling@grid[[1]]))
        names(sub.zoo.data.bic)[1] <- names(sub.zoo.data.qbic)[1] <- names(sub.zoo.data.aic)[1] <- "Series 1"
      }else{
        mod.bic <- setModel(drift = drif[[i]], diffusion = Diff.esti.bic, hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
        mod.qbic <- setModel(drift = drif[[i]], diffusion = Diff.esti.qbic, hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
        mod.aic <- setModel(drift = drif[[i]], diffusion = Diff.esti.aic, hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
        n <- nrow(data)-1
        modsamp <- setSampling(Terminal = Terminal, n = n)
        modyuima.bic <- setYuima(model = mod.bic, sampling = modsamp)
        modyuima.qbic <- setYuima(model = mod.qbic, sampling = modsamp)
        modyuima.aic <- setYuima(model = mod.bic, sampling = modsamp)
        sub.zoo.data.bic <- sub.zoo.data.qbic <- sub.zoo.data.aic <-  list()
        for(j in 1:ncol(data)){
          sub.zoo.data.bic <- c(sub.zoo.data.bic, list(zoo(x = data[,j], order.by = modyuima.bic@sampling@grid[[1]])))
          sub.zoo.data.qbic <- c(sub.zoo.data.qbic, list(zoo(x = data[,j], order.by = modyuima.qbic@sampling@grid[[1]])))
          sub.zoo.data.aic <- c(sub.zoo.data.aic, list(zoo(x = data[,j], order.by = modyuima.aic@sampling@grid[[1]])))
          names(sub.zoo.data.bic)[j] <- names(sub.zoo.data.qbic)[j] <- names(sub.zoo.data.aic)[j] <- paste("Series", j)
        }
      }
      modyuima.bic@data@zoo.data <- sub.zoo.data.bic
      modyuima.qbic@data@zoo.data <- sub.zoo.data.qbic
      modyuima.aic@data@zoo.data <- sub.zoo.data.aic
      yuimas2.bic <- c(yuimas2.bic, list(modyuima.bic))
      yuimas2.qbic <- c(yuimas2.qbic, list(modyuima.qbic))
      yuimas2.aic <- c(yuimas2.aic, list(modyuima.aic))
      
      ## Model comparison
      swalpha <- c(swalpha, list(modyuima.bic@model@parameter@drift))
      
      para.number.init  <- match(swalpha[[i]], names(start))
      para.number.low  <- match(swalpha[[i]], names(lower))
      para.number.upp  <- match(swalpha[[i]], names(upper))
      para.start <- NULL
      para.lower <- NULL
      para.upper <- NULL
      for(j in 1:length(swalpha[[i]])){
        para.start <- c(para.start, list(start[[para.number.init[j]]]))
        para.lower <- c(para.lower, list(lower[[para.number.low[j]]]))
        para.upper <- c(para.upper, list(upper[[para.number.upp[j]]]))
      }
      names(para.start) <- swalpha[[i]]
      names(para.lower) <- swalpha[[i]]
      names(para.upper) <- swalpha[[i]]
      
      mle.bic <- qmle(modyuima.bic, start = para.start, lower = para.lower, upper = para.upper, method = "L-BFGS-B", rcpp = rcpp)
      mle.qbic <- qmle(modyuima.qbic, start = para.start, lower = para.lower, upper = para.upper, method = "L-BFGS-B", rcpp = rcpp)
      mle.aic <- qmle(modyuima.aic, start = para.start, lower = para.lower, upper = para.upper, method = "L-BFGS-B", rcpp = rcpp)
      hess2 <- mle.qbic@details$hessian
      
      esti.bic <- list(coef(mle.bic))
      esti.qbic <- list(coef(mle.qbic))
      esti.aic <- list(coef(mle.aic))
      names(esti.bic[[1]]) <- names(esti.qbic[[1]]) <- names(esti.aic[[1]]) <- swalpha[[i]]
      bic <- summary(mle.bic)@m2logL+length(swalpha[[i]])*log(Terminal)
      if(det(hess2) > 0){
        qbic <- summary(mle.qbic)@m2logL+log(det(hess2))
      }else{
        qbic <- summary(mle.qbic)@m2logL+length(swalpha[[i]])*log(Terminal)
      }
      aic <- summary(mle.aic)@m2logL+2*length(swalpha[[i]])
      
      Esti2.bic <- c(Esti2.bic, esti.bic)
      Esti2.qbic <- c(Esti2.qbic, esti.qbic)
      Esti2.aic <- c(Esti2.aic, esti.aic)
      BIC2 <- c(BIC2, bic)
      QBIC2 <- c(QBIC2, qbic)
      AIC2 <- c(AIC2, aic)
    }
    BIC.opt2 <- which.min(BIC2)
    QBIC.opt2 <- which.min(QBIC2)
    AIC.opt2 <- which.min(AIC2)
    
    ## Names
    for(i in 1:length(drif)){
      names(Esti2.bic)[i] <- paste("drift", i, sep = "_")
      names(Esti2.qbic)[i] <- paste("drift", i, sep = "_")
      names(Esti2.aic)[i] <- paste("drift", i, sep = "_")
      names(BIC2)[i] <- paste("drift", i, sep = "_") 
      names(QBIC2)[i] <- paste("drift", i, sep = "_")
      names(AIC2)[i] <- paste("drift", i, sep = "_")
    }
    
    ## Model weights
    if(weight == TRUE){
      BIC.weight.full <- QBIC.weight.full <- AIC.weight.full <- matrix(0, length(drif), length(diff))
      for(i in 1:length(diff)){
        diff.row <- length(yuimas1[[i]]@model@diffusion)
        Esti1.chr <- as.character(Esti1[[i]])
        Diff.esti <- diff[[i]]
        for(j in 1:diff.row){
          if(length(Esti1.chr) == 1){
            Diff.esti.sub <- gsub(swbeta[[i]][1], Esti1.chr[1], yuimas1[[i]]@model@diffusion[[j]])
          }else{
            Diff.esti.sub <- gsub(swbeta[[i]][1], Esti1.chr[1], yuimas1[[i]]@model@diffusion[[j]])
            for(k in 1:(length(Esti1.chr)-1)){
              Diff.esti.sub <- gsub(swbeta[[i]][(k+1)], Esti1.chr[(k+1)], Diff.esti.sub)
            }
          }
          #if(class(Diff.esti) == "character"){
          if(inherits(Diff.esti, "character")){
            Diff.esti <- Diff.esti.sub
          } else {
            Diff.esti[j,] <- Diff.esti.sub
          }
        }
        
        BIC2.sub <- QBIC2.sub <- AIC2.sub <- NULL
        for(j in 1:length(drif)){
          if(is.matrix(data) == FALSE){
            mod <- setModel(drift = drif[[j]], diffusion = Diff.esti, hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
            n <- length(data)-1
            modsamp <- setSampling(Terminal = Terminal, n = n)
            modyuima <- setYuima(model = mod, sampling = modsamp)
            sub.zoo.data <- list(zoo(x = data, order.by = modyuima@sampling@grid[[1]]))
            names(sub.zoo.data)[1] <- "Series 1"
          }else{
            mod <- setModel(drift = drif[[j]], diffusion = Diff.esti, hurst = settings[[1]], measure = settings[[2]], measure.type = settings[[3]], state.variable = settings[[4]], jump.variable = settings[[5]], time.variable = settings[[6]], solve.variable = settings[[7]])
            n <- nrow(data)-1
            modsamp <- setSampling(Terminal = Terminal, n = n)
            modyuima <- setYuima(model = mod, sampling = modsamp)
            sub.zoo.data <-  list()
            for(k in 1:ncol(data)){
              sub.zoo.data <- c(sub.zoo.data, list(zoo(x = data[,k], order.by = modyuima@sampling@grid[[1]])))
              names(sub.zoo.data)[k] <- paste("Series", k)
            }
          }
          modyuima@data@zoo.data <- sub.zoo.data
          
          para.number.init  <- match(swalpha[[j]], names(start))
          para.number.low  <- match(swalpha[[j]], names(lower))
          para.number.upp  <- match(swalpha[[j]], names(upper))
          para.start <- NULL
          para.lower <- NULL
          para.upper <- NULL
          for(k in 1:length(swalpha[[j]])){
            para.start <- c(para.start, list(start[[para.number.init[k]]]))
            para.lower <- c(para.lower, list(lower[[para.number.low[k]]]))
            para.upper <- c(para.upper, list(upper[[para.number.upp[k]]]))
          }
          names(para.start) <- swalpha[[j]]
          names(para.lower) <- swalpha[[j]]
          names(para.upper) <- swalpha[[j]]
          
          mle.weight <- qmle(modyuima, start = para.start, lower = para.lower, upper = para.upper, method = "L-BFGS-B", rcpp = rcpp)
          hess.weight <- mle.weight@details$hessian
          
          esti.weight <- list(coef(mle.weight))
          names(esti.weight[[1]]) <- swalpha[[j]]
          bic <- summary(mle.weight)@m2logL+length(swalpha[[j]])*log(Terminal)
          if(det(hess.weight) > 0){
            qbic <- summary(mle.weight)@m2logL+log(det(hess.weight))
          }else{
            qbic <- summary(mle.weight)@m2logL+length(swalpha[[j]])*log(Terminal)
          }
          aic <- summary(mle.weight)@m2logL+2*length(swalpha[[j]])
          
          BIC2.sub <- c(BIC2.sub, bic)
          QBIC2.sub <- c(QBIC2.sub, qbic)
          AIC2.sub <- c(AIC2.sub, aic)
        }
        BIC2.sub.opt <- which.min(BIC2.sub)
        QBIC2.sub.opt <- which.min(QBIC2.sub)
        AIC2.sub.opt <- which.min(AIC2.sub)
        
        BIC.weight2 <- exp(-(1/2)*(BIC2.sub-BIC2.sub[BIC2.sub.opt]))/sum(exp(-(1/2)*(BIC2.sub-BIC2.sub[BIC2.sub.opt])))
        QBIC.weight2 <- exp(-(1/2)*(QBIC2.sub-QBIC2.sub[BIC2.sub.opt]))/sum(exp(-(1/2)*(QBIC2.sub-QBIC2.sub[QBIC2.sub.opt])))
        AIC.weight2 <- exp(-(1/2)*(AIC2.sub-AIC2.sub[AIC2.sub.opt]))/sum(exp(-(1/2)*(AIC2.sub-AIC2.sub[AIC2.sub.opt])))

        BIC.weight.full[,i] <- BIC.weight1[i]*BIC.weight2
        QBIC.weight.full[,i] <- QBIC.weight1[i]*QBIC.weight2
        AIC.weight.full[,i] <- AIC.weight1[i]*AIC.weight2
      }
      
      colname.weight <- numeric(length(diff))
      rowname.weight <- numeric(length(drif))
      for(i in 1:length(diff)){
        colname.weight[i] <- paste("scale", i, sep = "_") 
      }
      colnames(BIC.weight.full) <- colname.weight
      colnames(QBIC.weight.full) <- colname.weight
      colnames(AIC.weight.full) <- colname.weight
      for(i in 1:length(drif)){
        rowname.weight[i] <- paste("drift", i, sep = "_")
      }
      rownames(BIC.weight.full) <- rowname.weight
      rownames(QBIC.weight.full) <- rowname.weight
      rownames(AIC.weight.full) <- rowname.weight
    }
    
    ## Results
    diff.copy <- diff
    drif.copy <- drif
    for(i in 1:length(diff)){
      names(diff.copy)[i] <- paste("scale", i, sep = "_") 
    }
    for(i in 1:length(drif)){
      names(drif.copy)[i] <- paste("drift", i, sep = "_") 
    }
    BIC <- list(first = BIC1, second = BIC2)
    QBIC <- list(first = QBIC1, second = QBIC2)
    AIC <- list(first = AIC1, second = AIC2)
    if(is.matrix(data) == FALSE && Levy == TRUE){
      Esti <- list(first = Esti1, second.bic = NULL, second.qbic = Esti2.qbic, second.aic = Esti2.aic)
    }else{
      Esti <- list(first = Esti1, second.bic = Esti2.bic, second.qbic = Esti2.qbic, second.aic = Esti2.aic)
    }
    
    call <- match.call()
    model.coef <- list(drift = drif.copy, scale = diff.copy)
    bic.selected.coeff <- list(drift = drif[[BIC.opt2]], scale = diff[[BIC.opt1]])
    qbic.selected.coeff <- list(drift = drif[[QBIC.opt2]], scale = diff[[QBIC.opt1]])
    aic.selected.coeff <- list(drift = drif[[AIC.opt2]], scale = diff[[AIC.opt1]])
    if(is.matrix(data) == FALSE && Levy == TRUE){
      ic.selected <- list(BIC = NULL, QBIC = bic.selected.coeff, AIC = aic.selected.coeff)
    }else{
      ic.selected <- list(BIC = bic.selected.coeff, QBIC = qbic.selected.coeff, AIC = aic.selected.coeff)
    }
    if(weight == TRUE){
      if(is.matrix(data) == FALSE && Levy == TRUE){
        ak.weight <- list(BIC = NULL, QBIC = BIC.weight.full, AIC = AIC.weight.full)
      }else{
        ak.weight <- list(BIC = BIC.weight.full, QBIC = QBIC.weight.full, AIC = AIC.weight.full)
      }
    }else{
      ak.weight <- NULL
    }
    if(is.matrix(data) == FALSE && Levy == TRUE){
      final_res <- list(call = call, model = model.coef, par = Esti, BIC = NULL, QBIC = BIC, AIC = AIC, weight = ak.weight, selected = ic.selected)
    }else{
      final_res <- list(call = call, model = model.coef, par = Esti, BIC = BIC, QBIC = QBIC, AIC = AIC, weight = ak.weight, selected = ic.selected)
    }
  }

  class(final_res) <- "yuima.ic"
  return(final_res)
}

 print.yuima.ic <- function(x, ...){
  	cat("\nCall:\n")
  	print(x$call)
  	cat("\nInformation criteria:\n")
  	cat("\nBIC:\n")
  	print(x$BIC)
  	cat("\nQBIC:\n")
  	print(x$QBIC)
  	cat("\nAIC:\n")
  	print(x$AIC)
  	#if(class(x$AIC) == "matrix"){
  	#if(is.matrix(x$AIC)){ # fixed by YK
  	#	if(!is.null(x$AIC)){
  	#		cat("\nAIC:\n")
  	#		print(x$AIC)
  	#    }
  	#}
  	#if(class(x$AIC) == "list"){
  	#if(is.list(class(x$AIC))){ # fixed by YK
  	#	if(!is.null(x$AIC$first)){
  	#		cat("\nAIC:\n")
  	#		print(x$AIC)
  	#	}
    #}
    invisible(x)
 }
