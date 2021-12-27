########################################################################
#    Stepwise estimation for ergodic Levy driven SDE
########################################################################


#qmleLevy<-function(yuima,start,lower,upper,joint = FALSE,third = FALSE,
#                   Est.Incr = c("NoIncr","Incr","IncrPar"),
#                   aggregation = TRUE)
qmleLevy<-function(yuima,start,lower,upper,joint = FALSE,third = FALSE,
                   Est.Incr = "NoIncr",
                   aggregation = TRUE)
{
  oldyuima<-yuima #line1
  myjumpname <- yuima@model@jump.variable
  mymeasureparam <- yuima@model@parameter@measure
  if(!(Est.Incr %in% c("NoIncr","Incr","IncrPar")))
  	stop("Argument'Est.Incr' must be one of \"NoIncr\",\"Incr\" or \"IncrPar\"")
  call <- match.call()
  truestart<-start
  cat("\nStarting QGMLE for SDE ... \n")
  parameter<-yuima@model@parameter@all
  orig.mylaw<-yuima@model@measure
  mylaw<-yuima@model@measure$df
  numbLev<-length(yuima@model@measure.type)
  if(missing(yuima))
    yuima.stop("yuima object is missing.")
  
  if(!is(yuima,"yuima"))
    yuima.stop("This function is for yuima-class.")
  
  if(length(yuima@model@parameter@jump)>0)
    paracoef <- yuima@model@parameter@jump
  
  if(length(yuima@model@parameter@drift)>0)
    paracoef <- c(paracoef, yuima@model@parameter@drift)
  if(Est.Incr == "IncrPar"){
    start0<-start
    lower0<-lower
    upper0<-upper
    lev.names<-yuima@model@parameter@measure
  }
  
  DRIFT <- yuima@model@drift
  JUMP <- yuima@model@jump.coeff
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
    DiffHessian<- fres@details$hessian #182 
    
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
    DriftHessian <- sres@details$hessian #239
    
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
      coef<-c(sres@coef,fres@coef)
      mycoef<-unlist(truestart)
      #mycoef1<-mycoef[names(coef)]
      mycoef2<-mycoef[!names(mycoef)%in%names(coef)]
      mycoef<-c(coef,mycoef2)
      vcov0<-matrix(NA,nrow = length(coef),ncol=length(coef))
      rownames(vcov0)<-names(coef)
      colnames(vcov0)<-names(coef)
      min0<- c(fres@min,sres@min)
      details0<-list(sres@details,fres@details)
      nobs0<-sres@nobs
      res<-new("yuima.qmle", call = call, coef = coef, fullcoef = mycoef,
                     vcov = vcov0, min = min0, details = details0, minuslogl = minusquasilogl,
                     method = sres@method, nobs=nobs0, model=sdeModel)
      # res <- list(first = fres@coef, second = sres@coef)}
      if(length(oldyuima@data@zoo.data)==1){
        myGamhat <- matrix(0,length(coef),length(coef))
        myGamhat[1:dim(DiffHessian)[1],1:dim(DiffHessian)[2]]<-DiffHessian
        myGamhat[dim(DiffHessian)[1]+1:dim(DriftHessian)[1],dim(DiffHessian)[1]+1:dim(DriftHessian)[2]]<-DriftHessian#293
        myGamhat<-myGamhat/oldyuima@sampling@Terminal
      }
    }else{
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
    # jres<-qmle(yuima,start = start,lower = lower,upper = upper,rcpp = TRUE, joint = TRUE,method = "L-BFGS-B")
    # res<-list(joint = jres@coef)
    res<-qmle(yuima,start = start,lower = lower,upper = upper,rcpp = TRUE, joint = TRUE,
              method = "L-BFGS-B")
  }
  if(Est.Incr == "NoIncr"){
    return(res)
    }
  cat("\nStarting Estimation of Levy Increments ... \n")
  data <- get.zoo.data(yuima)
  s.size<-yuima@sampling@n 
  if(length(data)==1){
    X<-as.numeric(data[[1]])
    pX<-X[1:(s.size-1)]
    inc<-double(s.size-1)
    inc<-X[2:(s.size)]-pX
  }else{
    pX<- simplify2array(lapply(X = data, 
      FUN = as.numeric))
    pX<-pX[-nrow(pX),]
    inc<- simplify2array(lapply(X = data, 
      FUN = function(X){diff(as.numeric(X))}))
    
  }
  
  modeltime<-yuima@model@time.variable
  modelstate<-yuima@model@solve.variable
  tmp.env<-new.env()
  
  #if(joint){
    coeffic<- coef(res) 
  # }else{
  #   coeffic<- res[[1]]
  #   if(length(res)>1){
  #     for(j in c(2:length(res))){
  #       coeffic<-c(coeffic,res[[j]])
  #     }
  #   }
  #   
  #}
  mp<-match(names(coeffic),parameter)
  esort <- coeffic[order(mp)]
  for(i in 1:length(coeffic))
  {
    assign(parameter[i],esort[[i]],envir=tmp.env)
  }
  
  # DRIFT <- yuima@model@drift
  # JUMP <- yuima@model@jump.coeff
  
  if(length(yuima@model@solve.variable)==1){
    #parameter<-yuima@model@parameter@all
    
    
    #resi<-double(s.size-1) 
    assign(modeltime,yuima@sampling@delta,envir=tmp.env)
    h<-yuima@sampling@delta
    assign(modelstate,pX,envir=tmp.env)
    jump.term<-eval(JUMP[[1]],envir=tmp.env)
    drif.term<-eval(DRIFT,envir=tmp.env)
    if(length(jump.term)==1){
      jump.term <- rep(jump.term, s.size-1)
    }
    if(length(drif.term)==1){
      drif.term <- rep(drif.term, s.size-1)
    } # vectorization (note. if an expression type object does not include state.variable, the length of the item after "eval" operation is 1.)
    
    # for(s in 1:(s.size-1)){
    #   nova<-sqrt((jump.term)^2) # normalized variance
    #   resi[s]<-(1/(nova[s]))*(inc[s]-h*drif.term[s])
    # }
    nova<-sqrt((jump.term)^2)
    resi<-(1/(nova[1:(s.size-1)]))*(inc[1:(s.size-1)]-h*drif.term[1:(s.size-1)])
    if(length(oldyuima@data@zoo.data)==1){
      coefSigdiff<- 1/h*sum(resi^4)  #389 resi
      coefDriftSig <- 1/h*sum(resi^3)
    }
    
    if(aggregation){
      Ter <- yuima@sampling@Terminal
      ures <- numeric(floor(Ter))
      for(l in 1:floor(Ter)){
        ures[l] <- sum(resi[(floor((l-1)/h)):(floor(l/h)-1)]) 
      }
      res.incr<-ures
    }else{
      res.incr<-resi
      }
    
  }else{
    h<-yuima@sampling@delta
    
    Tbig<-dim(inc)[1]
    assign(modeltime,h,envir=tmp.env)
    numbofvar<- length(modelstate)
    for(j in c(1:numbofvar)){
      assign(modelstate[j],pX[,j],envir=tmp.env)
    }

    
    drif.term<-array(0,c(Tbig,numbofvar))
    for(i in c(1:numbofvar))
      drif.term [,i]<- eval(DRIFT[i],envir=tmp.env)
    # Check using variable in the drift
    
    # jump.term<-sapply(1:numbofvar,function(i){ 
    #    sapply(1:numbLev, function(j){
    #       eval(JUMP[[i]][j],envir=tmp.env)
    #      },simplify = TRUE)
    #    },simplify = TRUE)
    
    jump.term<-array(0,c(numbofvar,numbLev,Tbig))
    for(i in c(1:numbofvar)){
      for(j in c(1:numbLev))
        jump.term[i,j, ] <-  eval(JUMP[[i]][j],envir=tmp.env)
    }
    

   # if(dim(jump.term)[1]==numbofvar){
   #    if(dim(jump.term)[2]==numbofvar){
   #      if(det(jump.term)==0){
   #        Invjump.term<-solve(t(jump.term)%*%jump.term)%*%t(jump.term)
   #      }else{
   #        Invjump.term<-solve(jump.term)
   #      }
   #    }else{
   #      Invjump.term<-solve(t(jump.term)%*%jump.term)%*%t(jump.term)
   #    }
   #  Invjump.term <- Invjump.term%o%rep(1,Tbig)
   # }else{
   #   
   # }
   # 
  
  DeltaInc<-(inc-drif.term*h)
  if(dim(jump.term)[1]==numbofvar){
    if(dim(jump.term)[2]==numbofvar){
        resi<-t(sapply(i:Tbig,function(i){
            if(det(jump.term[,,i])==0){
              step1<-t(jump.term[,,i])%*%jump.term[,,i]
              if(det(step1)==0){
                Invjump.term<-diag(rep(1,dim(step1)[1]))
              }else{
                Invjump.term<-solve(step1)%*%t(jump.term[,,i])
              }
            }else{
              Invjump.term<-solve(jump.term[,,i])
            }
               Invjump.term%*% DeltaInc[i,]
          }
          )
          )
    }else{
      resi<-t(sapply(i:Tbig,function(i){
          step1<-t(jump.term[,,i])%*%jump.term[,,i]
          if(det(step1)==0){
            Invjump.term<-diag(rep(1,dim(step1)[1]))
          }else{
            Invjump.term<-solve(step1)%*%t(jump.term[,,i])
          }
          Invjump.term%*% DeltaInc[i,]
          }
        )
      )
    }
  }

    if(aggregation){
      Ter <- min(floor(yuima@sampling@Terminal))
      
      res.incr<-t(sapply(1:Ter,function(i) colSums(resi[(floor((i-1)/h)):(floor(i/h)-1),]) ))
      
    }else{
      res.incr<-resi
    }
  }
  
  
  if(aggregation){
    if(!is.matrix(res.incr)){
      res.incr<- as.matrix(res.incr)
    }
    if(dim(res.incr)[2]==1){
      colnames(res.incr)<-sdeModel@jump.variable
    }else{
      colnames(res.incr)<-paste0(sdeModel@jump.variable,c(1:dim(res.incr)[2]))
    }
    Incr.Lev <- zooreg(data=res.incr)
    Incr.Lev<- setData(original.data = Incr.Lev)
  }else{
    Incr.Lev <- zoo(res.incr,order.by=yuima@sampling@grid[[1]][-1])
    Incr.Lev <- setData(original.data=Incr.Lev)
  }
  
  if(length(oldyuima@data@zoo.data)==1){
    mydiff<-oldyuima@model@jump.coeff[[1]]
    mydiffDer <-deriv(mydiff,oldyuima@model@parameter@jump)
    myenvdiff<- new.env()
    if(length(oldyuima@model@parameter@jump)>=1){
      for(i in c(1:length(oldyuima@model@parameter@jump))){
        assign(value=coef[oldyuima@model@parameter@jump[i]],x=oldyuima@model@parameter@jump[i],envir=myenvdiff)
      }
    }
    EvalPartDiff <- Vectorize(FUN= function(myenvdiff,mydiffDer, data){
      assign(x=oldyuima@model@solve.variable, value=data,envir = myenvdiff)
      return(attr(eval(mydiffDer, envir=myenvdiff),"gradient"))},vectorize.args = "data") 
    
    DiffJumpCoeff<-EvalPartDiff(myenvdiff,mydiffDer, data=pX)
    if(!is.matrix(DiffJumpCoeff)){
      sigmadiffhat<- as.matrix(sum(DiffJumpCoeff^2/jump.term[1:(oldyuima@sampling@n-1)]^2)/(oldyuima@sampling@n))*coefSigdiff
      DiffJumpCoeff<- t(DiffJumpCoeff)
    }else{
      sigmadiffhat<- matrix(0,dim(DiffJumpCoeff)[1],dim(DiffJumpCoeff)[1])
      for(t in c(1:dim(DiffJumpCoeff)[2])){
        sigmadiffhat <-sigmadiffhat+as.matrix(DiffJumpCoeff[,t])%*%DiffJumpCoeff[,t]/jump.term[t]^2
      }
      sigmadiffhat<-  sigmadiffhat/(oldyuima@sampling@n)*coefSigdiff
    }
    
    mydrift<-oldyuima@model@drift[[1]]
    mydriftDer <-deriv(mydrift,oldyuima@model@parameter@drift)
    
    myenvdrift<- new.env()
    if(length(oldyuima@model@parameter@drift)>=1){
      for(i in c(1:length(oldyuima@model@parameter@drift))){
        assign(value=coef[oldyuima@model@parameter@drift[i]],x=oldyuima@model@parameter@drift[i],envir=myenvdrift)
      }
    }
    
    DriftDerCoeff<-EvalPartDiff(myenvdrift,mydriftDer, data=pX)
    if(!is.matrix(DriftDerCoeff)){
     # sigmadrifthat<- as.matrix(sum(DriftDerCoeff^2/jump.term[1:(oldyuima@sampling@n-1)]^2)/(oldyuima@sampling@n))*coefDriftSig
      sigmadrifthat<- as.matrix(sum(DriftDerCoeff^2/jump.term[1:(oldyuima@sampling@n-1)]^2)/(oldyuima@sampling@n))
      DriftDerCoeff <- t(DriftDerCoeff)
    }else{
      sigmadrifthat<- matrix(0,dim(DriftDerCoeff)[1],dim(DriftDerCoeff)[1])
      for(t in c(1:dim(DriftDerCoeff)[2])){
        sigmadrifthat <-sigmadrifthat+as.matrix(DriftDerCoeff[,t])%*%DriftDerCoeff[,t]/jump.term[t]^2
      }
      sigmadrifthat<- sigmadrifthat/(oldyuima@sampling@n)
    }
    sigmadriftdiff <- matrix(0, dim(sigmadrifthat)[2],dim(sigmadiffhat)[1])
    for(t in c(1:dim(DriftDerCoeff)[2]))
      sigmadriftdiff<-sigmadriftdiff+DriftDerCoeff[,t]%*%t(DiffJumpCoeff[,t])
    
    #sigmadriftdiff<-sigmadriftdiff/oldyuima@sampling@n
    sigmadriftdiff<-sigmadriftdiff/oldyuima@sampling@n*coefDriftSig
    
    MatSigmaHat <- rbind(cbind(sigmadiffhat,t(sigmadriftdiff)),cbind(sigmadriftdiff,sigmadrifthat))
    
    res@coef<-res@coef[c(oldyuima@model@parameter@jump,oldyuima@model@parameter@drift)]
    res@vcov<-solve(t(myGamhat)%*%solve(MatSigmaHat)%*%myGamhat*oldyuima@sampling@Terminal)
  }
  
  if(Est.Incr == "Incr"){
    
      
    result<- new("yuima.qmleLevy.incr",Incr.Lev=Incr.Lev,
                 Data = yuima@data,  yuima=res)
    return(result)
  }
  cat("\nEstimation Levy parameters ... \n")
  
  if(class(mylaw)=="yuima.law"){
    if(aggregation){
      minusloglik <- function(para){
        para[length(para)+1]<-1
        names(para)[length(para)]<-yuima@model@time.variable
        -sum(dens(object=mylaw, x=res.incr, param = para, log = TRUE), 
           na.rm = T)
        }
      }else{
        minusloglik <- function(para){
          para[length(para)+1] <- yuima@sampling@delta
          names(para)[length(para)]<-yuima@model@time.variable
          -sum(dens(object=mylaw, x=res.incr, param = para, log = TRUE), 
               na.rm = T)
        }
      }
    para <- start0[lev.names]
    lowerjump <- lower0[lev.names]
    upperjump <- upper0[lev.names]
    esti <- optim(fn = minusloglik, lower = lowerjump, upper = upperjump, 
                  par = para, method = "L-BFGS-B")

    HessianEta <- optimHess(par=esti$par, fn=minusloglik)
    res@coef<-c(res@coef,esti$par)
    res@fullcoef[names(para)]<-esti$par
    if(length(oldyuima@data@zoo.data)==1 & is(oldyuima@model@measure$df, "yuima.law")){
      if(!aggregation){
        Ter <- yuima@sampling@Terminal
        ures <- numeric(floor(Ter))
        for(l in 1:floor(Ter)){
          ures[l] <- sum(resi[(floor((l-1)/h)):(floor(l/h)-1)]) 
        } 
      }
      mypar<-res@coef[oldyuima@model@parameter@measure]
      mypar[oldyuima@model@measure$df@time.var]<-1
      fdataeta<-dens(object=oldyuima@model@measure$df,x=ures,param=mypar)# f(eps, eta)
      del <- 10^-3
      fdatadeltaeta<- dens(object=oldyuima@model@measure$df,x=ures+del,param=mypar) # f(eps + delta, eta)
      dummy <- t(rep(1,length(mypar[oldyuima@model@measure$df@param.measure])))
      myeta<- as.matrix(mypar[oldyuima@model@measure$df@param.measure])%*%dummy
      myetapert <- myeta+del*diag(rep(1,dim(myeta)[1])) 
      fdataetadelta<-sapply(X=1:dim(myeta)[1],FUN = function(X){
        par<-myetapert[,X]
        par[oldyuima@model@measure$df@time.var]<-1
        dens(object=oldyuima@model@measure$df,x=ures,param=par)
      }
      )# f(eps, eta+delta)
      fdatadeltaetadelta<-sapply(X=1:dim(myeta)[1],FUN = function(X){
        par<-myetapert[,X]
        par[oldyuima@model@measure$df@time.var]<-1
        dens(object=oldyuima@model@measure$df,x=ures+del,param=par)
      }
      ) # f(eps +deta, eta+delta) 
      term1<-1/(fdataeta)
      term2 <- fdatadeltaeta*term1%*%dummy
      
      term2 <- fdatadeltaetadelta-term2*fdataetadelta
      mixpartial<-t((as.matrix(term1)%*%dummy)/del^2*term2/oldyuima@sampling@Terminal)
      
      
      # construction of b_i 
      # DiffJumpCoeff, DriftDerCoeff, jump.term length(resi)
      # step1 <- t(DiffJumpCoeff)%*%DiffJumpCoeff
      DerMeta <- 1/del*(fdataetadelta - as.matrix(fdataeta)%*%rep(1,dim(fdataetadelta)[2]))*(as.matrix(term1)%*%rep(1,dim(fdataetadelta)[2]))
      SigmaEta <- t(DerMeta)%*%DerMeta/oldyuima@sampling@Terminal
      
      #SigmaEtaAlpha<- 1/oldyuima@sampling@n*DriftDerCoeff%*%(t(DiffJumpCoeff)/(as.matrix(jump.term[-length(jump.term)]^2)%*%rep(1,dim(DiffJumpCoeff)[1])) )
      SigmaEtaAlpha<- 1/oldyuima@sampling@n*DriftDerCoeff%*%(t(DiffJumpCoeff)/(as.matrix(jump.term^2)%*%rep(1,dim(DiffJumpCoeff)[1])) )
      SigmaEtaAlpha <- SigmaEtaAlpha*sum(resi^3)/oldyuima@sampling@delta
      
      b_i <- matrix(0,floor(Ter),length(c(oldyuima@model@parameter@drift, oldyuima@model@parameter@jump)))
      Coef1 <- matrix(0,floor(Ter),length( oldyuima@model@parameter@jump))
      Coef2 <- matrix(0,floor(Ter),length( oldyuima@model@parameter@drift))
      
      for(l in 1:floor(Ter)){
        pos <- (floor((l-1)/h)):(floor(l/h)-1)
        if(length(oldyuima@model@parameter@jump)==1){
          b_i[l,1:length(oldyuima@model@parameter@jump)] <- sum(-DiffJumpCoeff[,pos]*resi[pos]/jump.term[pos])
          Coef1[l,]<-sum(DiffJumpCoeff[,pos]/jump.term[pos]*(resi[pos]^2-h))
        }else{
          interm <- as.matrix(resi[pos]/jump.term[pos])
          b_i[l,1:length(oldyuima@model@parameter@jump)] <-  -t(DiffJumpCoeff[,pos]%*%interm)
          interm2 <- as.matrix((resi[pos]^2-h)/jump.term[pos])
          Coef1[l,] <-DiffJumpCoeff[,pos]%*%interm2
        }
        if(length(oldyuima@model@parameter@drift)==1){
          b_i[l,1:length(oldyuima@model@parameter@drift)+length(oldyuima@model@parameter@jump)] <- -sum(h* DriftDerCoeff[,pos]%*%jump.term[pos])
          Coef2[l,]<-sum(DriftDerCoeff[,pos]/jump.term[pos]*(resi[pos]))
        }else{
          b_i[l,1:length(oldyuima@model@parameter@drift)+length(oldyuima@model@parameter@jump)] <--h* DriftDerCoeff[,pos]%*%jump.term[pos]
          interm3 <- as.matrix((resi[pos])/jump.term[pos])
          Coef2[l,] <-DriftDerCoeff[,pos]%*%interm3
        }
      }
      MatrUnder <- mixpartial%*%b_i
      I_n <- cbind(rbind(myGamhat,MatrUnder),rbind(matrix(0, dim(myGamhat)[1],dim(HessianEta)[2]),HessianEta/oldyuima@sampling@Terminal))
      SigmaGammaEta <- t(DerMeta)%*%Coef1/Ter
      SigmaAlphaEta <- t(DerMeta)%*%Coef2/Ter
      # dim(fdataetadelta), length(fdataeta),  length(term1)
      dum <- cbind(SigmaGammaEta , SigmaAlphaEta)
      MatSigmaHat1<- rbind(cbind(MatSigmaHat,t(dum)),cbind(dum,SigmaEta)) 
      InvIn<- solve(I_n)
      res@vcov <-InvIn%*%MatSigmaHat1%*%t(InvIn)/Ter
    }else{
      res@vcov<-cbind(res@vcov,matrix(NA,ncol=length(esti$par),nrow=dim(res@vcov)[1]))
      res@vcov<-rbind(res@vcov,matrix(NA,nrow=length(esti$par),ncol=dim(res@vcov)[2]))
    }
    colnames(res@vcov)<-names(res@fullcoef)
    #res@vcov<-rbind(res@vcov,matrix(NA,nrow=length(esti$par),ncol=dim(res@vcov)[2]))
    rownames(res@vcov)<-names(res@fullcoef)
    res@min<-c(res@min,esti$value)
    res@nobs<-c(res@nobs,length(Incr.Lev@zoo.data[[1]]))
    
    result<- new("yuima.qmleLevy.incr",Incr.Lev=Incr.Lev,
                 minusloglLevy = minusloglik,logL.Incr=-esti$value,
                 Data = yuima@data,  yuima=res, Levydetails=esti)
    
    return(result)
    }else{
      dist <- substr(as.character(orig.mylaw$df$expr), 2, 10^3)
  
      startjump <- start0[lev.names]
      lowerjump <- lower0[lev.names]
      upperjump <- upper0[lev.names]
  
      if(length(startjump) == -1){
        logdens <- function(para){
          exlogdens <- parse(text = sprintf("log(d%s)", dist))
          assign(myjumpname, ures, envir = tmp.env)
          assign(mymeasureparam, para, envir = tmp.env)
          sum(eval(exlogdens, envir = tmp.env))
        }
        mydens <-function(para){
          exdens <- parse(text = sprintf("d%s", dist))
          assign(myjumpname, ures, envir = tmp.env)
          assign(mymeasureparam, para, envir = tmp.env)
          eval(exdens, envir = tmp.env)
        }
        intervaljump <- c(lowerjump[[1]], upperjump[[1]])
        esti <- optimize(logdens, interval = intervaljump, maximum = TRUE)
        return(list(sde=esort, meas=esti$maximum))
      }else{
        
        logdens <- function(para){
          exlogdens <- parse(text = sprintf("log(d%s)", dist))
          assign(oldyuima@model@jump.variable, ures, envir = tmp.env)
          for(i in c(1:length(oldyuima@model@parameter@measure)))
            assign(oldyuima@model@parameter@measure[i], para[[oldyuima@model@parameter@measure[i]]], envir = tmp.env)
  
          -sum(eval(exlogdens, envir = tmp.env),na.rm = T)
        }
        
        mydens1 <-function(par,x){
          exdens <- parse(text = sprintf("d%s", dist))
          assign(myjumpname, x, envir = tmp.env)
          assign(mymeasureparam, par, envir = tmp.env)
          eval(exdens, envir = tmp.env)
        }
  
        esti <- optim(fn=logdens, lower = lowerjump, upper = upperjump, par = startjump,
                      method = "L-BFGS-B")
        
        HessianEta <- optimHess(par=esti$par, fn=logdens)
        res@coef<-c(res@coef,esti$par)
        res@fullcoef[lev.names]<-esti$par
        if(length(oldyuima@data@zoo.data)==1){
          if(!aggregation){
            Ter <- yuima@sampling@Terminal
            ures <- numeric(floor(Ter))
            for(l in 1:floor(Ter)){
              ures[l] <- sum(resi[(floor((l-1)/h)):(floor(l/h)-1)]) 
            } 
          }
          mypar<-res@coef[oldyuima@model@parameter@measure]
          #mypar[oldyuima@model@measure$df@time.var]<-1
          fdataeta<-mydens1(par=mypar,x=ures)# f(eps, eta)
          del <- 10^-3
          fdatadeltaeta<- mydens1(par=mypar,x=ures+del) # f(eps + delta, eta)
          dummy <- t(rep(1,length(mypar[lev.names])))
          myeta<- as.matrix(mypar[lev.names])%*%dummy
          myetapert <- myeta+del*diag(rep(1,dim(myeta)[1])) 
          fdataetadelta<-sapply(X=1:dim(myeta)[1],FUN = function(X){
            par<-myetapert[,X]
            #par[oldyuima@model@measure$df@time.var]<-1
            mydens1(par=par,x=ures)
          }
          )# f(eps, eta+delta)
          fdatadeltaetadelta<-sapply(X=1:dim(myeta)[1],FUN = function(X){
            par<-myetapert[,X]
            #par[oldyuima@model@measure$df@time.var]<-1
            mydens1(par=par,x=ures+del)
          }
          ) # f(eps +deta, eta+delta) 
          term1<-1/(fdataeta)
          term2 <- fdatadeltaeta*term1%*%dummy
          
          term2 <- fdatadeltaetadelta-term2*fdataetadelta
          mixpartial<-t((as.matrix(term1)%*%dummy)/del^2*term2/oldyuima@sampling@Terminal)
          
          
          # construction of b_i 
          # DiffJumpCoeff, DriftDerCoeff, jump.term length(resi)
          # step1 <- t(DiffJumpCoeff)%*%DiffJumpCoeff
          DerMeta <- 1/del*(fdataetadelta - as.matrix(fdataeta)%*%rep(1,dim(fdataetadelta)[2]))*(as.matrix(term1)%*%rep(1,dim(fdataetadelta)[2]))
          SigmaEta <- t(DerMeta)%*%DerMeta/oldyuima@sampling@Terminal
          
          #SigmaEtaAlpha<- 1/oldyuima@sampling@n*DriftDerCoeff%*%(t(DiffJumpCoeff)/(as.matrix(jump.term[-length(jump.term)]^2)%*%rep(1,dim(DiffJumpCoeff)[1])) )
          SigmaEtaAlpha<- 1/oldyuima@sampling@n*DriftDerCoeff%*%(t(DiffJumpCoeff)/(as.matrix(jump.term^2)%*%rep(1,dim(DiffJumpCoeff)[1])) )
          SigmaEtaAlpha <- SigmaEtaAlpha*sum(resi^3)/oldyuima@sampling@delta
          
          b_i <- matrix(0,floor(Ter),length(c(oldyuima@model@parameter@drift, oldyuima@model@parameter@jump)))
          Coef1 <- matrix(0,floor(Ter),length( oldyuima@model@parameter@jump))
          Coef2 <- matrix(0,floor(Ter),length( oldyuima@model@parameter@drift))
          
          for(l in 1:floor(Ter)){
            pos <- (floor((l-1)/h)):(floor(l/h)-1)
            if(length(oldyuima@model@parameter@jump)==1){
              b_i[l,1:length(oldyuima@model@parameter@jump)] <- sum(-DiffJumpCoeff[,pos]*resi[pos]/jump.term[pos])
              Coef1[l,]<-sum(DiffJumpCoeff[,pos]/jump.term[pos]*(resi[pos]^2-h))
            }else{
              interm <- as.matrix(resi[pos]/jump.term[pos])
              b_i[l,1:length(oldyuima@model@parameter@jump)] <-  -t(DiffJumpCoeff[,pos]%*%interm)
              interm2 <- as.matrix((resi[pos]^2-h)/jump.term[pos])
              Coef1[l,] <-DiffJumpCoeff[,pos]%*%interm2
            }
            if(length(oldyuima@model@parameter@drift)==1){
              b_i[l,1:length(oldyuima@model@parameter@drift)+length(oldyuima@model@parameter@jump)] <- -sum(h* DriftDerCoeff[,pos]%*%jump.term[pos])
              Coef2[l,]<-sum(DriftDerCoeff[,pos]/jump.term[pos]*(resi[pos]))
            }else{
              b_i[l,1:length(oldyuima@model@parameter@drift)+length(oldyuima@model@parameter@jump)] <--h* DriftDerCoeff[,pos]%*%jump.term[pos]
              interm3 <- as.matrix((resi[pos])/jump.term[pos])
              Coef2[l,] <-DriftDerCoeff[,pos]%*%interm3
            }
          }
          MatrUnder <- mixpartial%*%b_i
          I_n <- cbind(rbind(myGamhat,MatrUnder),rbind(matrix(0, dim(myGamhat)[1],dim(HessianEta)[2]),HessianEta/oldyuima@sampling@Terminal))
          SigmaGammaEta <- t(DerMeta)%*%Coef1/Ter
          SigmaAlphaEta <- t(DerMeta)%*%Coef2/Ter
          # dim(fdataetadelta), length(fdataeta),  length(term1)
          dum <- cbind(SigmaGammaEta , SigmaAlphaEta)
          MatSigmaHat1<- rbind(cbind(MatSigmaHat,t(dum)),cbind(dum,SigmaEta)) 
          InvIn<- solve(I_n)
          res@vcov <-InvIn%*%MatSigmaHat1%*%t(InvIn)/Ter
        }else{
          res@vcov<-cbind(res@vcov,matrix(NA,ncol=length(esti$par),nrow=dim(res@vcov)[1]))
        }
        colnames(res@vcov)<-names(res@fullcoef)
        #res@vcov<-rbind(res@vcov,matrix(NA,nrow=length(esti$par),ncol=dim(res@vcov)[2]))
        rownames(res@vcov)<-names(res@fullcoef)
        res@min<-c(res@min,esti$value)
        res@nobs<-c(res@nobs,length(Incr.Lev@zoo.data[[1]]))
        
        result<- new("yuima.qmleLevy.incr",Incr.Lev=Incr.Lev,
                     minusloglLevy = logdens,logL.Incr=-esti$value,
                     Data = yuima@data,  yuima=res, Levydetails=esti)
        
        return(result)
        
        
        #return(list(sde=esort, meas=esti$par))
      }
    }
  
}


