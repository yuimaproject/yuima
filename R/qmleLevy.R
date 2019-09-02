########################################################################
#    Stepwise estimation for ergodic Levy driven SDE
########################################################################


qmleLevy<-function(yuima,start,lower,upper,joint = FALSE,third = FALSE)
{
  if(missing(yuima))
    yuima.stop("yuima object is missing.")
  
  if(!is(yuima,"yuima"))
    yuima.stop("This function is for yuima-class.")
  
  sdeModel<-yuima@model
  if(length(sdeModel@parameter@measure)!=0){
    nPar<-length(sdeModel@parameter@measure)
    for(i in c(1:nPar)){
      assign(x = sdeModel@parameter@measure[i],
             value = start[[sdeModel@parameter@measure[i]]])
    }
    names1 <- names(start)
    index <- which(names1 %in% sdeModel@parameter@measure)
    start <- start[-index]
    names1 <- names(lower)
    index <- which(names1 %in% sdeModel@parameter@measure)
    lower <- lower[-index]
    names1 <- names(upper)
    index <- which(names1 %in% sdeModel@parameter@measure)
    upper <- upper[-index]
  }
  if(class(sdeModel@measure$df)!="yuima.law"){
    code <- suppressWarnings(sub("^(.+?)\\(.+", "\\1", sdeModel@measure$df$expr, perl=TRUE))
    
    candinoise<-c("rNIG","rvgamma","rnts","rbgamma")
    
    
    
    if(is.na(match(code,candinoise))){
      yuima.stop("This function works only for the standardized normal inverse Gaussian process, variance gamma process, bilateral gamma process, and normal tempered stable process now.")
    }
    
    
    if(length(sdeModel@xinit) == 1){
      args <- unlist(strsplit(suppressWarnings(sub("^.+?\\((.+)\\)", "\\1", sdeModel@measure$df$expr, perl=TRUE)), ","))
    if(code == "rNIG"){
        if(!((abs(eval(parse(text = paste("(",args[5] ,")+(", args[3], ")*(", args[4],
                                            ")/sqrt((", args[2], ")^2-(", args[3], ")^2)" )))) < 10^(-10))
             && (abs(eval(parse(text = paste("(",args[2], ")^2*(", args[4], ")/(sqrt((",  args[2],
                                               ")^2-(", args[3], ")^2))^3"))) -1) < 10^(-10))))
        {
          yuima.stop("This function is only for standardized Levy noises.")
        }
      }
      else if(code == "rvgamma"){
        if(!((abs(eval(parse(text = paste("(", args[5], ")+2*(", args[2], ")*(", args[4], ")/((", args[3], ")^2-("
                                            , args[4], ")^2)" )))) < 10^(-10))
             && (abs(eval(parse(text = paste("2*((", args[2], ")*(", args[3], ")^2+(", args[4], ")^2)", "/(("
                                               , args[3], ")^2-(", args[4], ")^2)^2"))) - 1) < 10^(-10))))
        {
          yuima.stop("This function is only for standardized Levy noises")
        }
      }
    else if((code == "rnts")){
        if(!((abs(eval(parse(text = paste("(", args[6], ")-(", args[2], ")*(",
                                            args[3], ")*(", args[4], ")^((", args[2],
                                            ")-1)*gamma(1-(", args[2], "))*(-1/(", args[2],"))*(", args[5],
                                            ")"
        )))) < 10^(-10))
        &&(abs(eval(parse(text = paste("(", args[3], ")*(", args[2],")*((", args[2], ")-1)*gamma(1-(", args[2], "))*(-1/(", args[2],"))*(", args[4], ")^((",args[2], ")-2)*(", args[5], ")^2-(",
                                         args[2], ")*(", args[3], ")*(", args[4], ")^((", args[2], ")-1)*gamma(1-(", args[2], "))*(-1/(", args[2],"))"
        ))) - 1) < 10^(-10))))
        {
          yuima.stop("This function is only for standardized Levy processes.")
      }
    }else if(code == "rbgamma"){
      if(!((abs(eval(parse(text = paste("(", args[2], ")/(", args[3],")-(", args[4],")/(", args[5], ")")))) < 10^(-10))
           && (abs(eval(parse(text = paste("(", args[2], ")/(", args[3], ")^2","+(", args[4],")/(", args[5],")^2"
           ))) - 1) < 10^(-10) )))
      {
        yuima.stop("This function is only for standardized Levy processes.")
      }
    }
    }else{
      warning("In this version, the standardized conditions on multidimensional noises can not be verified.
              The expressions of mean and variance are given in help page.")
      # The noise condition checker below does not work now (YU: 3/23).
      
      # args <- suppressWarnings(sub("^.+?\\((.+)\\)", "\\1", sdeModel@measure$df$expr, perl=TRUE))
      # yuimaEnv <- new.env()
      # yuimaEnv$mean <- switch(code,
      #                         rNIG = function(x=1,alpha,beta,delta0,mu,Lambda){mu+as.vector(delta0/(sqrt(alpha^2-t(beta)%*%Lambda%*%beta)))*Lambda%*%beta},
      #                         rnts = function(x=1,alpha,a,b,beta,mu,Lambda){mu+gamma(1-alpha)*a*b^(alpha-1)*Lambda%*%beta},
      #                         rvgamma = function(x=1,lambda,alpha,beta,mu,Lambda){mu+as.vector(2*lambda/(alpha^2-t(beta)%*%Lambda%*%beta)^2)*beta}
      #                         )
      # 
      # yuimaEnv$covariance <- switch(code,
      #                               rNIG = function(x=1,alpha,beta,delta0,mu,Lambda){as.vector(delta0/(sqrt(alpha^2-t(beta)%*%Lambda%*%beta))^3)*Lambda%*%beta%*%t(beta)%*%Lambda+as.vector(delta0/sqrt(alpha^2-t(beta)%*%Lambda%*%beta))*Lambda},
      #                               rnts = function(x=1,alpha,a,b,beta,mu,Lambda){a*(1-alpha)*gamma(1-alpha)*b^(alpha-2)*Lambda%*%beta%*%t(beta)%*%Lambda+a*gamma(1-alpha)*b^(alpha-1)*Lambda},
      #                               rvgamma = function(x=1,lambda,alpha,beta,mu,Lambda){as.vector(4*lambda/(alpha^2-t(beta)%*%Lambda%*%beta)^2)*Lambda%*%beta%*%t(beta)%*%Lambda+as.vector(2*lambda/alpha^2-t(beta)%*%Lambda%*%beta)*Lambda}
      #                               )
      # judgemean<-sum(eval(parse(text = paste("mean","(",args,")")),yuimaEnv)==numeric(length(sdeModel@xinit)))
      # judgecovariance<-sum(eval(eval(parse(text = paste("covariance","(",args,")")),yuimaEnv)==diag(1,length(sdeModel@xinit))))
      # if(!((judgemean==length(sdeModel@xinit))&&(judgecovariance==length(sdeModel@xinit)*length(sdeModel@xinit))))
      #    {
      #   yuima.stop("This function is only for standardized Levy processes.")
      # }
    }
  }else{fullcoef<-NULL}
  yuima@sampling@delta <- yuima@sampling@delta[1]
  yuima@model@noise.number <- as.integer(yuima@model@equation.number)
  if(!joint){
    DRIFT <- yuima@model@drift
    DRPAR <- yuima@model@parameter@drift
    
    if(length(yuima@model@parameter@jump)>0)
      fullcoef <- yuima@model@parameter@jump
    
    if(length(DRPAR)>0)
      fullcoef <- c(fullcoef, DRPAR)
    
    oo <- match(yuima@model@parameter@all, fullcoef)
    yuima@model@parameter@all <- yuima@model@parameter@all[order(oo)]
    
    oo <- match(names(start), fullcoef)
    start <- start[order(oo)]
    
    oo <- match(names(upper), fullcoef)
    upper <- upper[order(oo)]
    
    oo <- match(names(lower), fullcoef)
    lower <- lower[order(oo)]
    
    yuima@model@diffusion <- yuima@model@jump.coeff
    yuima@model@parameter@diffusion <- yuima@model@parameter@jump[1:length(yuima@model@parameter@jump)]
    yuima@model@parameter@all <- yuima@model@parameter@diffusion
    
    for(i in 1:length(yuima@model@drift)){
      yuima@model@drift[i] <- expression((0))
    }
    
    yuima@model@jump.coeff <- list()
    yuima@model@parameter@drift <- character(0)
    yuima@model@measure <- list()
    yuima@model@jump.variable <- character(0)
    yuima@model@measure.type <- character(0)
    yuima@model@parameter@jump <- character(0)
    yuima@model@parameter@measure <- character(0)
    
    diffstart <- start[1:length(yuima@model@parameter@diffusion)]
    diffupper <- upper[1:length(yuima@model@parameter@diffusion)]
    difflower <- lower[1:length(yuima@model@parameter@diffusion)]
    fres <- qmle(yuima=yuima,start=diffstart,lower=difflower,upper=diffupper,rcpp = TRUE,joint = FALSE,method = "L-BFGS-B")
    
    DIPAR <- yuima@model@parameter@diffusion
    DIFFUSION <- yuima@model@diffusion
    yuima@model@parameter@all <- DRPAR
    yuima@model@parameter@drift <- DRPAR
    yuima@model@drift <- DRIFT
    
    dristart <- start[-(1:length(yuima@model@parameter@diffusion))]
    driupper <- upper[-(1:length(yuima@model@parameter@diffusion))]
    drilower <- lower[-(1:length(yuima@model@parameter@diffusion))]
    
    partcoef <- yuima@model@parameter@diffusion
    ov <- match(yuima@model@parameter@drift,partcoef)
    ovp <- which(!is.na(ov))
    if(length(ovp)>0)
    {yuima@model@parameter@drift <- yuima@model@parameter@drift[-ovp]}
    yuima@model@parameter@all <- yuima@model@parameter@drift
    
    ma <- match(names(fres@coef),partcoef)
    sort <- fres@coef[order(ma)]
    esti <- numeric(length(partcoef))
    newdiff <- yuima@model@diffusion
    newdri <- yuima@model@drift
    
    for(i in 1:length(partcoef))
    {
      esti[i] <- as.character(fres@coef[[i]])
    }
    
    if(length(yuima@model@drift) == 1){
      for(i in 1:length(partcoef))
      {
        newdri <- gsub(partcoef[i],esti[i],newdri)
        yuima@model@drift[1] <- parse(text = newdri[1])
        newdiff[[1]] <- gsub(partcoef[i],esti[i],newdiff[[1]])
        yuima@model@diffusion[[1]] <- parse(text = newdiff[[1]])
      }
    }else{
      for(i in 1:length(partcoef))
      {
        for(j in 1:length(yuima@model@drift))
        {
          newdri[j] <- gsub(partcoef[i],esti[i],newdri[j])
          yuima@model@drift[j] <- parse(text = newdri[j])
        }
        for(k in 1:length(yuima@model@diffusion)){
          for(l in 1:length(yuima@model@diffusion[[1]])){
            newdiff[[k]][l] <- gsub(partcoef[i],esti[i],newdiff[[k]][l])
            yuima@model@diffusion[[k]][l] <- parse(text = newdiff[[k]][l])
          }
        }
      }
    }
    yuima@model@parameter@diffusion <- character(0)
    
    sres<-qmle(yuima=yuima,start=dristart,lower=drilower,upper=driupper,rcpp = TRUE,method = "L-BFGS-B")
    
    if((length(ovp) == 0) && (third)){
      yuima@model@diffusion <- DIFFUSION
      yuima@model@drift <- DRIFT
      yuima@model@parameter@diffusion <- DIPAR
      yuima@model@parameter@all <- DIPAR
      newdri <- yuima@model@drift
      
      for(i in 1:length(sres@coef))
      {
        esti[i] <- as.character(sres@coef[[i]])
      }
      
      if(length(yuima@model@drift)==1){
        for(i in 1:length(sres@coef))
        {
          newdri <- gsub(yuima@model@parameter@drift[i],esti[i],newdri)
          yuima@model@drift[1] <- parse(text=newdri[1])
        }
      }else{
        for(i in 1:length(sres@coef))
        {
          for(j in 1:length(yuima@model@drift))
          {
            newdri[j] <- gsub(yuima@model@parameter@drift[i],esti[i],newdri[j])
            yuima@model@drift[j] <- parse(text = newdri[j])
          }
        }
      }
      
      yuima@model@parameter@drift <- character(0)
      too <- match(names(fres@coef),names(diffstart))
      diffstart <- diffstart[order(too)]
      for(i in 1:length(diffstart)){
        diffstart[[i]] <- fres@coef[[i]]
      }
      tres <- qmle(yuima=yuima,start=diffstart,lower=difflower,upper=diffupper,rcpp = TRUE,joint = FALSE,method = "L-BFGS-B")
      
      res <- list(first = fres@coef, second = sres@coef, third = tres@coef)
    }else if((length(ovp) > 0) || !(third)){
      res <- list(first = fres@coef, second = sres@coef)}
    else{
      yuima.stop("third estimation may be theoretical invalid under the presence of an overlapping parameter.")
    }
  }else{
    if(third){
      yuima.stop("third estimation does not make sense in joint estimation.")
    }
    if(length(yuima@model@parameter@jump)>0)
      fullcoef <- yuima@model@parameter@jump
    if(length(yuima@model@parameter@drift)>0)
      fullcoef <- c(fullcoef, yuima@model@parameter@drift)
    oo <- match(yuima@model@parameter@all, fullcoef)
    yuima@model@parameter@all <- yuima@model@parameter@all[order(oo)]
    yuima@model@parameter@all <- yuima@model@parameter@all[1:length(which(!is.na(oo)))]
    yuima@model@diffusion <- yuima@model@jump.coeff
    yuima@model@jump.coeff <- list()
    yuima@model@parameter@diffusion <- yuima@model@parameter@jump[1:length(yuima@model@parameter@jump)]
    yuima@model@measure <- list()
    yuima@model@jump.variable <- character(0)
    yuima@model@measure.type <- character(0)
    yuima@model@parameter@jump <- character(0)
    yuima@model@parameter@measure <- character(0)
    jres<-qmle(yuima,start = start,lower = lower,upper = upper,rcpp = TRUE, joint = TRUE,method = "L-BFGS-B")
    res<-list(joint = jres@coef)
  }
  res
}


