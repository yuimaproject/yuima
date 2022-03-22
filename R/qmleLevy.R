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
  #if(class(sdeModel@measure$df)!="yuima.law"){
  if(!inherits(sdeModel@measure$df, "yuima.law")){ # fixed by YK
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
      }else{
        myGamhat <- matrix(0,length(coef),length(coef))
        myGamhat[1:dim(DiffHessian)[1],1:dim(DiffHessian)[2]]<-DiffHessian*oldyuima@sampling@Terminal/oldyuima@sampling@n
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
    
    if(length(oldyuima@data@zoo.data)==1){  
      result<- new("yuima.qmleLevy.incr",Incr.Lev=Incr.Lev,
                 Data = yuima@data,  yuima=res)
    }else{
      result<- new("yuima.qmleLevy.incr",Incr.Lev=Incr.Lev,
                   Data = yuima@data,  yuima=res)
      
      vcovLevyNoMeas <- function(myres,  Gammahat0, sq=TRUE){
        #myres <- res.VG2
        DeltaX <- apply(myres@Data@original.data,2,diff)
        myData<- myres@Data@original.data
        CmatExpr<- myres@model@jump.coeff
        ncolC <-length(myres@model@jump.coeff[[1]])
        param <- myres@coef[myres@model@parameter@jump]
        aexpr<- myres@model@drift
        namedrift<-myres@model@parameter@drift
        pardrif <- myres@coef[namedrift]
        avect_exp<-myres@model@drift
        
        Jac_Drift <- function(aexpr,namedrift,nobs=length(aexpr)){
          lapply(X=c(1:nobs),FUN=function(X,namedrift,aexpr){
            return(deriv(aexpr[X],namedrift))
          },namedrift=namedrift,aexpr=aexpr)
        }
        
        Jac_DriftExp <- Jac_Drift(aexpr=avect_exp,namedrift)  
        
        FUNDum<-function(foo,myenv) sapply(foo, function(x,env) eval(x,envir = env), env= myenv)
        
        h<-diff(time(myres@Data@zoo.data[[1]]))[1]
        
        dummyf1n<- function(DeltaX, Cmat, h){
          C2<-t(Cmat)%*%Cmat
          dec <- chol(C2) # Eventualy Add a tryCatch
          tmp <- t(DeltaX)%*%solve(C2)%*%DeltaX
          logretval <- -h*sum(log(diag(dec)))  - 0.5 * as.numeric(tmp)
          return(logretval)
        }
        
        dummyInc<- function(DeltaX,Cmat, avect,h,sq=TRUE){
          if(sq){
            Incr <- solve(Cmat)%*%(DeltaX-avect*h)
          }else{
            C2<-t(Cmat)%*%Cmat
            Incr <- solve(C2)%*%Cmat%*%(DeltaX-avect*h)
          }
          return(Incr)
        }
        
        dummyf2n<- function(DeltaX, avect, Cmat, h){
          C2 <- t(Cmat)%*%Cmat
          Incr <- DeltaX-h*avect
          tmp <- t(Incr)%*%solve(C2)%*%Incr
          logretval <-   - 0.5/h * as.numeric(tmp)
          return(logretval)
        }
        
        dumGrad_f2n<-function(DeltaX, avect, Jac_a,Cmat, h){
          C2 <- t(Cmat)%*%Cmat
          Incr <- DeltaX-h*avect
          Grad<-t(Jac_a)%*%solve(C2)%*%Incr
          return(Grad)
        }
        
        f1n_j<-function(x, myObsDelta, CmatExpr,h, myData, myres){
          names(myData)<-myres@model@solve.variable
          newenvJumpCoef <- list2env(as.list(c(x,myData)))
          Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenvJumpCoef )
          return(dummyf1n(DeltaX=myObsDelta, Cmat=Cmat, h=h))
        }
        
        Incr_Func <- function(x, myObsDelta, CmatExpr,avect_exp, h, myData, myres, sq=TRUE){
          names(myData)<-myres@model@solve.variable
          newenv <- list2env(as.list(c(x,myData)))
          Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenv )
          dd<-length(avect_exp)
          avect<-numeric(length=dd)
          for(j in c(1:dd)){
            avect[j]<-eval(avect_exp[j],envir=newenv)
          }
          return(dummyInc(DeltaX= myObsDelta,Cmat, avect,h,sq=sq))
        }
        
        
        f2n_j <- function(x, parDiff, myObsDelta, CmatExpr, avect_exp, h, myData, myres){
          names(myData)<-myres@model@solve.variable
          newenvJumpCoef <- list2env(as.list(c(parDiff,myData)))
          Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenvJumpCoef)
          newenvDriftCoef <- list2env(as.list(c(x,myData)))
          dd<-length(avect_exp)
          avect<-numeric(length=dd)
          for(j in c(1:dd)){
            avect[j]<-eval(avect_exp[j],envir=newenvDriftCoef)
          }
          return(dummyf2n(DeltaX=myObsDelta, avect=avect, Cmat=Cmat, h=h))
        }
        
        Gradf2n_j <-function(x, parDiff, myObsDelta, CmatExpr, 
                             avect_exp, Jac_DriftExp, h, myData, myres){
          names(myData)<-myres@model@solve.variable
          newenvJumpCoef <- list2env(as.list(c(parDiff,myData)))
          Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenvJumpCoef)
          newenvDriftCoef <- list2env(as.list(c(x,myData)))
          dd<-length(avect_exp)
          avect<-numeric(length=dd)
          Jac_a<- matrix(0,dd,length(x))
          for(j in c(1:dd)){
            avect[j]<-eval(avect_exp[j],envir=newenvDriftCoef)
            Jac_a[j,]<-attr(eval(Jac_DriftExp[[j]],envir=newenvDriftCoef),"gradient")
          }
          return(dumGrad_f2n(DeltaX=myObsDelta, avect=avect, Jac_a=Jac_a, Cmat=Cmat, h=h))
        }
        
        # f1n_j(x=param, myObsDelta=DeltaX[1,], CmatExpr=CmatExpr, h=h)
        #i=1
        # debug(dummyInc)
        # Incr_Func(x=c(param,pardrif), myObsDelta=DeltaX[i,], CmatExpr,avect_exp, h, myData=myData[i,], myres, sq=TRUE)
        # f1n_j(x=param, myObsDelta=DeltaX[i,], CmatExpr=CmatExpr, h=h,  myData=myData[i,], myres=myres)
        # f2n_j(x, parDiff=param, myObsDelta=DeltaX[i,], CmatExpr=CmatExpr, avect_exp=avect_exp, 
        #       h=h, myData=myData[i,], myres=myres)
        # Gradf2n_j(x, parDiff=param, myObsDelta=DeltaX[i,], CmatExpr, 
        #     avect_exp, Jac_DriftExp, h, myData[i,], myres)
        del <- 10^-3
        dummy <- t(rep(1,length(param)))
        myeta<- as.matrix(param)%*%dummy
        dummyG <- t(rep(1,length(c(param,pardrif))))
        globalmyeta <- as.matrix(c(param,pardrif))%*%dummyG
        Incr <- del*diag(rep(1,dim(myeta)[1]))
        GlobIncr <- del*diag(rep(1,dim(globalmyeta)[1]))
        Sigma_gamma <- matrix(0 ,length(param),length(param))
        Sigma_alpha <- matrix(0 ,length(pardrif),length(pardrif))
        Sigma_algam <- matrix(0 ,length(pardrif),length(param))
        histDeltaf1n_j_delta <-matrix(0 , dim(DeltaX)[1],length(param))
        histDeltaf2n_j_delta <- matrix(0 , dim(DeltaX)[1],length(pardrif))
        #  histb_incr <- array(0, c(dim(DeltaX)[2],dim(DeltaX)[1],length(c(param,pardrif))))
        for(i in c(1:dim(DeltaX)[1])){
          Deltaf1n_j_delta<-sapply(X=1:dim(myeta)[1],
                                   FUN = function(X,myObsDelta, CmatExpr, h, myData, myres,del){
                                     par1<-myeta[,X]+Incr[,X]
                                     #par[oldyuima@model@measure$df@time.var]<-1
                                     f1<-f1n_j(x=par1, myObsDelta, CmatExpr, h, myData, myres)
                                     par2<-myeta[,X]-Incr[,X]
                                     f2<-f1n_j(x=par2, myObsDelta, CmatExpr, h, myData, myres)
                                     return((f1-f2)/(2*del))
                                   },
                                   myObsDelta=DeltaX[i,], CmatExpr=CmatExpr, h=h, myData=myData[i,],
                                   myres=myres, del=del)
          
          # DeltaInc <- sapply(X=1:dim(GlobIncr)[1],
          #                    FUN = function(X, myObsDelta, CmatExpr,avect_exp, h, 
          #                                   myData, myres, del, sq){
          #                      par1<-globalmyeta[,X]+GlobIncr[,X]
          #                      
          #                      # f1<-f1n_j(x=par1, myObsDelta, CmatExpr, h, myData, myres)
          #                      f1<-Incr_Func(x=par1, myObsDelta, CmatExpr,avect_exp, h, myData, myres, sq)
          #                      par2<-globalmyeta[,X]-GlobIncr[,X]
          #                      f2<-Incr_Func(x=par2, myObsDelta, CmatExpr,avect_exp, h, myData, myres, sq)
          #                      #f2<-f1n_j(x=par2, myObsDelta, CmatExpr, h, myData, myres)
          #                      return((f1-f2)/(2*del))
          #                    },
          #                    myObsDelta=DeltaX[i,], CmatExpr,avect_exp, h, myData=myData[i,], myres, del=del,sq=sq
          # )
          # 
          # histDeltaf1n_j_delta[i, ]<-Deltaf1n_j_delta
          
          #histb_incr[ , i, ] <- DeltaInc
          Sigma_gamma <- Sigma_gamma + as.matrix(Deltaf1n_j_delta)%*%t(Deltaf1n_j_delta)
          
          Deltaf2n_j_delta<- Gradf2n_j(x=pardrif, parDiff=param, myObsDelta=DeltaX[i,], CmatExpr, 
                                       avect_exp, Jac_DriftExp, h, myData[i,], myres)
          histDeltaf2n_j_delta[i,]<-Deltaf2n_j_delta
          
          Sigma_alpha <- Sigma_alpha + Deltaf2n_j_delta%*%t(Deltaf2n_j_delta)
          Sigma_algam <- Sigma_algam + Deltaf2n_j_delta%*%t(Deltaf1n_j_delta)
          #cat("\n",i)
        }
        Tn<-tail(time(myres@Data@original.data),1L)
        Sigma_gamma0<-Sigma_gamma/Tn
        Sigma_alpha0 <- Sigma_alpha/Tn     
        Sigma_algam0 <- Sigma_algam/Tn
        Sigma0 <- cbind(rbind(Sigma_gamma0,Sigma_algam0),rbind(t(Sigma_algam0),Sigma_alpha0))
        
        InvGammaHAT <- solve(Gammahat0)
        vcov <- InvGammaHAT %*% Sigma0 %*% InvGammaHAT/Tn
        
        myres@vcov<- vcov
        return(myres)
      }
      
      if(length(result@model@jump.coeff)==length(result@model@jump.coeff[[1]])){
        sq<-TRUE
      }else{
        sq<-FALSE
      }
      
      result<-vcovLevyNoMeas(myres=result,  Gammahat0=myGamhat, sq=sq)
    }
    return(result)
  }
  cat("\nEstimation Levy parameters ... \n")
  
  #if(class(mylaw)=="yuima.law"){
  if(inherits(mylaw, "yuima.law")){ # YK, Mar. 22, 2022
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
      colnames(res@vcov)<-names(res@fullcoef)
      #res@vcov<-rbind(res@vcov,matrix(NA,nrow=length(esti$par),ncol=dim(res@vcov)[2]))
      rownames(res@vcov)<-names(res@fullcoef)
      res@min<-c(res@min,esti$value)
      res@nobs<-c(res@nobs,length(Incr.Lev@zoo.data[[1]]))
      result<- new("yuima.qmleLevy.incr",Incr.Lev=Incr.Lev,
                   minusloglLevy = minusloglik,logL.Incr=-esti$value,
                   Data = yuima@data,  yuima=res, Levydetails=esti)
    }else{
      Tn<-yuima@sampling@Terminal
      HessianEta_divTn <- HessianEta/Tn
      if(!aggregation){
        Ter <- min(floor(yuima@sampling@Terminal))
        ures<-t(sapply(1:Ter,function(i) colSums(resi[(floor((i-1)/h)):(floor(i/h)-1),]) ))
        #yuima.stop("da fare")
        vcovLevy1 <- function(myres,   HessianEta_divTn,  Gammahat0, ures, sq=TRUE){
          #myres <- res.VG2
          DeltaX <- apply(myres@Data@original.data,2,diff)
          myData<- myres@Data@original.data
          CmatExpr<- myres@model@jump.coeff
          ncolC <-length(myres@model@jump.coeff[[1]])
          param <- myres@coef[myres@model@parameter@jump]
          aexpr<- myres@model@drift
          namedrift<-myres@model@parameter@drift
          pardrif <- myres@coef[namedrift]
          avect_exp<-myres@model@drift
          
          Jac_Drift <- function(aexpr,namedrift,nobs=length(aexpr)){
            lapply(X=c(1:nobs),FUN=function(X,namedrift,aexpr){
              return(deriv(aexpr[X],namedrift))
            },namedrift=namedrift,aexpr=aexpr)
          }
          
          Jac_DriftExp <- Jac_Drift(aexpr=avect_exp,namedrift)  
          
          FUNDum<-function(foo,myenv) sapply(foo, function(x,env) eval(x,envir = env), env= myenv)
          
          h<-diff(time(myres@Data@zoo.data[[1]]))[1]
          
          dummyf1n<- function(DeltaX, Cmat, h){
            C2<-t(Cmat)%*%Cmat
            dec <- chol(C2) # Eventualy Add a tryCatch
            tmp <- t(DeltaX)%*%solve(C2)%*%DeltaX
            logretval <- -h*sum(log(diag(dec)))  - 0.5 * as.numeric(tmp)
            return(logretval)
          }
          dummyInc<- function(DeltaX,Cmat, avect,h,sq=TRUE){
            if(sq){
              Incr <- solve(Cmat)%*%(DeltaX-avect*h)
            }else{
              C2<-t(Cmat)%*%Cmat
              Incr <- solve(C2)%*%Cmat%*%(DeltaX-avect*h)
            }
            return(Incr)
          }
          
          dummyf2n<- function(DeltaX, avect, Cmat, h){
            C2 <- t(Cmat)%*%Cmat
            Incr <- DeltaX-h*avect
            tmp <- t(Incr)%*%solve(C2)%*%Incr
            logretval <-   - 0.5/h * as.numeric(tmp)
            return(logretval)
          }
          dumGrad_f2n<-function(DeltaX, avect, Jac_a,Cmat, h){
            C2 <- t(Cmat)%*%Cmat
            Incr <- DeltaX-h*avect
            Grad<-t(Jac_a)%*%solve(C2)%*%Incr
            return(Grad)
          }
          
          f1n_j<-function(x, myObsDelta, CmatExpr,h, myData, myres){
            names(myData)<-myres@model@solve.variable
            newenvJumpCoef <- list2env(as.list(c(x,myData)))
            Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenvJumpCoef )
            return(dummyf1n(DeltaX=myObsDelta, Cmat=Cmat, h=h))
          }
          
          Incr_Func <- function(x, myObsDelta, CmatExpr,avect_exp, h, myData, myres, sq=TRUE){
            names(myData)<-myres@model@solve.variable
            newenv <- list2env(as.list(c(x,myData)))
            Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenv )
            dd<-length(avect_exp)
            avect<-numeric(length=dd)
            for(j in c(1:dd)){
              avect[j]<-eval(avect_exp[j],envir=newenv)
            }
            return(dummyInc(DeltaX= myObsDelta,Cmat, avect,h,sq=sq))
          }
          
          
          f2n_j <- function(x, parDiff, myObsDelta, CmatExpr, avect_exp, h, myData, myres){
            names(myData)<-myres@model@solve.variable
            newenvJumpCoef <- list2env(as.list(c(parDiff,myData)))
            Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenvJumpCoef)
            newenvDriftCoef <- list2env(as.list(c(x,myData)))
            dd<-length(avect_exp)
            avect<-numeric(length=dd)
            for(j in c(1:dd)){
              avect[j]<-eval(avect_exp[j],envir=newenvDriftCoef)
            }
            return(dummyf2n(DeltaX=myObsDelta, avect=avect, Cmat=Cmat, h=h))
          }
          
          Gradf2n_j <-function(x, parDiff, myObsDelta, CmatExpr, 
                               avect_exp, Jac_DriftExp, h, myData, myres){
            names(myData)<-myres@model@solve.variable
            newenvJumpCoef <- list2env(as.list(c(parDiff,myData)))
            Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenvJumpCoef)
            newenvDriftCoef <- list2env(as.list(c(x,myData)))
            dd<-length(avect_exp)
            avect<-numeric(length=dd)
            Jac_a<- matrix(0,dd,length(x))
            for(j in c(1:dd)){
              avect[j]<-eval(avect_exp[j],envir=newenvDriftCoef)
              Jac_a[j,]<-attr(eval(Jac_DriftExp[[j]],envir=newenvDriftCoef),"gradient")
            }
            return(dumGrad_f2n(DeltaX=myObsDelta, avect=avect, Jac_a=Jac_a, Cmat=Cmat, h=h))
          }
          
          # f1n_j(x=param, myObsDelta=DeltaX[1,], CmatExpr=CmatExpr, h=h)
          #i=1
          # debug(dummyInc)
          # Incr_Func(x=c(param,pardrif), myObsDelta=DeltaX[i,], CmatExpr,avect_exp, h, myData=myData[i,], myres, sq=TRUE)
          # f1n_j(x=param, myObsDelta=DeltaX[i,], CmatExpr=CmatExpr, h=h,  myData=myData[i,], myres=myres)
          # f2n_j(x, parDiff=param, myObsDelta=DeltaX[i,], CmatExpr=CmatExpr, avect_exp=avect_exp, 
          #       h=h, myData=myData[i,], myres=myres)
          # Gradf2n_j(x, parDiff=param, myObsDelta=DeltaX[i,], CmatExpr, 
          #     avect_exp, Jac_DriftExp, h, myData[i,], myres)
          del <- 10^-3
          dummy <- t(rep(1,length(param)))
          myeta<- as.matrix(param)%*%dummy
          dummyG <- t(rep(1,length(c(param,pardrif))))
          globalmyeta <- as.matrix(c(param,pardrif))%*%dummyG
          Incr <- del*diag(rep(1,dim(myeta)[1]))
          GlobIncr <- del*diag(rep(1,dim(globalmyeta)[1]))
          Sigma_gamma <- matrix(0 ,length(param),length(param))
          Sigma_alpha <- matrix(0 ,length(pardrif),length(pardrif))
          Sigma_algam <- matrix(0 ,length(pardrif),length(param))
          histDeltaf1n_j_delta <-matrix(0 , dim(DeltaX)[1],length(param))
          histDeltaf2n_j_delta <- matrix(0 , dim(DeltaX)[1],length(pardrif))
          histb_incr <- array(0, c(dim(DeltaX)[2],dim(DeltaX)[1],length(c(param,pardrif))))
          for(i in c(1:dim(DeltaX)[1])){
            Deltaf1n_j_delta<-sapply(X=1:dim(myeta)[1],
                                     FUN = function(X,myObsDelta, CmatExpr, h, myData, myres,del){
                                       par1<-myeta[,X]+Incr[,X]
                                       #par[oldyuima@model@measure$df@time.var]<-1
                                       f1<-f1n_j(x=par1, myObsDelta, CmatExpr, h, myData, myres)
                                       par2<-myeta[,X]-Incr[,X]
                                       f2<-f1n_j(x=par2, myObsDelta, CmatExpr, h, myData, myres)
                                       return((f1-f2)/(2*del))
                                     },
                                     myObsDelta=DeltaX[i,], CmatExpr=CmatExpr, h=h, myData=myData[i,],
                                     myres=myres, del=del)
            
            DeltaInc <- sapply(X=1:dim(GlobIncr)[1],
                               FUN = function(X, myObsDelta, CmatExpr,avect_exp, h, 
                                              myData, myres, del, sq){
                                 par1<-globalmyeta[,X]+GlobIncr[,X]
                                 
                                 # f1<-f1n_j(x=par1, myObsDelta, CmatExpr, h, myData, myres)
                                 f1<-Incr_Func(x=par1, myObsDelta, CmatExpr,avect_exp, h, myData, myres, sq)
                                 par2<-globalmyeta[,X]-GlobIncr[,X]
                                 f2<-Incr_Func(x=par2, myObsDelta, CmatExpr,avect_exp, h, myData, myres, sq)
                                 #f2<-f1n_j(x=par2, myObsDelta, CmatExpr, h, myData, myres)
                                 return((f1-f2)/(2*del))
                               },
                               myObsDelta=DeltaX[i,], CmatExpr,avect_exp, h, myData=myData[i,], myres, del=del,sq=sq
            )
            
            histDeltaf1n_j_delta[i, ]<-Deltaf1n_j_delta
            
            histb_incr[ , i, ] <- DeltaInc
            Sigma_gamma <- Sigma_gamma + as.matrix(Deltaf1n_j_delta)%*%t(Deltaf1n_j_delta)
            
            Deltaf2n_j_delta<- Gradf2n_j(x=pardrif, parDiff=param, myObsDelta=DeltaX[i,], CmatExpr, 
                                         avect_exp, Jac_DriftExp, h, myData[i,], myres)
            histDeltaf2n_j_delta[i,]<-Deltaf2n_j_delta
            
            Sigma_alpha <- Sigma_alpha + Deltaf2n_j_delta%*%t(Deltaf2n_j_delta)
            Sigma_algam <- Sigma_algam + Deltaf2n_j_delta%*%t(Deltaf1n_j_delta)
            #cat("\n",i)
          }
          Tn<-tail(time(myres@Data@original.data),1L)
          Sigma_gamma0<-Sigma_gamma/Tn
          Sigma_alpha0 <- Sigma_alpha/Tn     
          Sigma_algam0 <- Sigma_algam/Tn
          Sigma0 <- cbind(rbind(Sigma_gamma0,Sigma_algam0),rbind(t(Sigma_algam0),Sigma_alpha0))
          
          myLogDens <- function(x, myObsDelta, CmatExpr,h,myData, myres, pos){
            return(f1n_j(x, myObsDelta[pos, ], CmatExpr,h,myData=myData[pos,], myres=myres))
          }
          
          myLogDens2<- function(x, parDiff, myObsDelta, CmatExpr, avect_exp, h, myData, myres, pos){
            return(f2n_j(x, parDiff, myObsDelta[pos, ], CmatExpr, avect_exp, h, myData[pos,], myres))
            #return(f1n_j(x, myObsDelta[pos, ], CmatExpr,h,myData=myData[pos,], myres=myres))
          }
          
          
          VectLogDens<-Vectorize(FUN=myLogDens,vectorize.args = "pos")
          #VectLogDens2<-Vectorize(FUN=myLogDens2,vectorize.args = "pos")
          
          mylogLik<- function(par,DeltaX,CmatExpr, myData, myres, h, Tn, nobs=dim(DeltaX)[1]){
            term<-VectLogDens(x=par, 
                              myObsDelta=DeltaX, CmatExpr=CmatExpr,h=h, myData=myData, myres=myres,pos=c(1:nobs))
            return(sum(term)/Tn)
          }
          
          mylogLik2<- function(par, parDiff, myObsDelta, CmatExpr, avect_exp, h, Tn, myData, myres, nobs=dim(DeltaX)[1]){
            # term<-VectLogDens2(x=par, parDiff, myObsDelta, CmatExpr, avect_exp, h, myData, myres,
            #                    pos=c(1:nobs))
            term<-0
            for(j in  c(1:nobs)){
              term<-term+f2n_j(x=par, parDiff, myObsDelta[j, ], CmatExpr, avect_exp, h, myData[j,], myres)
            }
            return(term/Tn)
          }
          
          GradlogLik2<-function(par, parDiff, myObsDelta, CmatExpr, avect_exp, h, Tn, myData, myres, nobs=dim(DeltaX)[1]){
            term<-matrix(0, length(par),1)
            for(j in  c(1:nobs)){ 
              term <- term + Gradf2n_j(x=par, parDiff, myObsDelta=myObsDelta[j, ], CmatExpr, 
                                       avect_exp, Jac_DriftExp, h, myData=myData[j,], myres)
            }
            return(as.numeric(term)/Tn)
          }
          
          # Gammahat <- optimHess(par=param,fn=mylogLik,
          #                       DeltaX=DeltaX,CmatExpr=CmatExpr,
          #                       h=h, Tn=Tn, nobs=dim(DeltaX)[1], 
          #                       myData=myData, myres=myres)
          
          # GammaAlfahat <- optimHess(par=pardrif,fn=mylogLik2,gr=GradlogLik2,
          #                           parDiff=param, myObsDelta=DeltaX, 
          #                           CmatExpr=CmatExpr, avect_exp=avect_exp, h=h,
          #                           Tn=Tn, nobs=dim(DeltaX)[1], myData=myData, myres=myres)
          # 
          # Gammahat0 <- cbind(rbind(Gammahat,matrix(0,dim(GammaAlfahat)[1],dim(Gammahat)[2])),
          #                    rbind(matrix(0,dim(Gammahat)[2],dim(GammaAlfahat)[1]),GammaAlfahat))
          
          uhistDeltaf1n_j<- matrix(0,floor(Tn),length(param))
          uhistDeltaf2n_j<- matrix(0,floor(Tn),length(pardrif))
          aaa<-length(c(param,pardrif))
          uhistb_incr <- array(0, c(dim(DeltaX)[2],floor(Tn),aaa))
          for(l in 1:floor(Tn)){
            #ures[l] <- sum(resi[(floor((l-1)/h)):(floor(l/h)-1)]) 
            uhistDeltaf1n_j[l, ] <- colSums(histDeltaf1n_j_delta[(floor((l-1)/h)):(floor(l/h)-1),])
            uhistDeltaf2n_j[l, ] <- colSums(histDeltaf2n_j_delta[(floor((l-1)/h)):(floor(l/h)-1),])
            for(i in c(1:aaa)){
              dumIn <- histb_incr[,(floor((l-1)/h)):(floor(l/h)-1),i]
              uhistb_incr[,l,i]<- rowSums(dumIn)
            }
          }
          
          nMeaspar<-myres@model@parameter@measure
          mypar_meas0 <- myres@coef[nMeaspar]#res@coef[oldyuima@model@parameter@measure]
          mypar_meas <- c(mypar_meas0,1)
          # ures<-myres@Incr.Lev@original.data
          names(mypar_meas)<-c(nMeaspar,myres@model@measure$df@time.var)
          fdataeta<-dens(object=myres@model@measure$df,x=ures,param=mypar_meas)# f(eps, eta)
          del <- 10^-3
          dummy <- t(rep(1,length(mypar_meas[myres@model@measure$df@param.measure])))
          myeta<- as.matrix(mypar_meas[myres@model@measure$df@param.measure])%*%dummy
          myetapert <- myeta+del*diag(rep(1,dim(myeta)[1])) 
          fdataetadelta<-sapply(X=1:dim(myeta)[1],FUN = function(X){
            par<-myetapert[,X]
            par[myres@model@measure$df@time.var]<-1
            dens(object=myres@model@measure$df,x=ures,param=par)
          }
          )# f(eps, eta+delta)
          
          term1<-1/(fdataeta)
          
          dummyvecdel<-numeric(length=dim(ures)[2]) 
          fdatadeltaeta<- matrix(0,floor(Tn),dim(ures)[2])
          mixpartial <-array(0,c(length(nMeaspar),floor(Tn),dim(ures)[2]))
          for(j in c(1:dim(ures)[2])){
            dummyvecdel[j]<-del
            fdatadeltaeta[,j]<- dens(object=myres@model@measure$df,x=ures+rep(1,dim(ures)[1])%*%t(dummyvecdel),param=mypar_meas)
            
            fdatadeltaetadelta<-sapply(X=1:dim(myeta)[1],FUN = function(X){
              par<-myetapert[,X]
              par[myres@model@measure$df@time.var]<-1
              dens(object=myres@model@measure$df,x=ures+rep(1,dim(ures)[1])%*%t(dummyvecdel),param=par)
            }
            ) # f(eps +deta, eta+delta)
            dummyvecdel<-numeric(length=dim(ures)[2])
            term2 <- fdatadeltaeta[,j]*term1%*%dummy
            term2 <- fdatadeltaetadelta-term2*fdataetadelta
            mixpartial[,,j]<-t((as.matrix(term1)%*%dummy)/del^2*term2/Tn)
          }
          
          DerMeta <- 1/del*(fdataetadelta - as.matrix(fdataeta)%*%rep(1,dim(fdataetadelta)[2]))*(as.matrix(term1)%*%rep(1,dim(fdataetadelta)[2]))
          SigmaEta <- t(DerMeta)%*%DerMeta/Tn
          
          minusloglik <- function(para){
            para[length(para)+1]<-1
            names(para)[length(para)]<-myres@model@time.variable
            -sum(dens(object=myres@model@measure$df, x=ures, param = para, log = TRUE), 
                 na.rm = T)/Tn
          }
          #  HessianEta_divTn <- optimHess(par=mypar_meas0, fn=minusloglik)
          
          Gammaeta_theta <- matrix(0,length(mypar_meas0), aaa)
          for(t in c(1:floor(Tn))){
            Gammaeta_theta<-Gammaeta_theta+mixpartial[,t,]%*%uhistb_incr[,t,]
          }
          GammaEta_Theta<- -1/Tn*Gammaeta_theta
          
          Sigma_eta_theta<-t(DerMeta)%*%cbind(uhistDeltaf1n_j,uhistDeltaf2n_j)/Tn
          
          Sigma <- cbind(rbind(Sigma0,Sigma_eta_theta),rbind(t(Sigma_eta_theta),SigmaEta))
          
          GammaHAT<-cbind(rbind(Gammahat0, GammaEta_Theta), rbind(t(GammaEta_Theta), HessianEta_divTn))
          InvGammaHAT <- solve(GammaHAT)
          vcov <- InvGammaHAT %*% Sigma %*% InvGammaHAT/Tn
          
          myres@vcov<- vcov
          return(myres)
        }
      }
      
        res@vcov<-cbind(res@vcov,matrix(NA,ncol=length(esti$par),nrow=dim(res@vcov)[1]))
        res@vcov<-rbind(res@vcov,matrix(NA,nrow=length(esti$par),ncol=dim(res@vcov)[2]))
        colnames(res@vcov)<-names(res@fullcoef)
        #res@vcov<-rbind(res@vcov,matrix(NA,nrow=length(esti$par),ncol=dim(res@vcov)[2]))
        rownames(res@vcov)<-names(res@fullcoef)
        res@min<-c(res@min,esti$value)
        res@nobs<-c(res@nobs,length(Incr.Lev@zoo.data[[1]]))
      
        result<- new("yuima.qmleLevy.incr",Incr.Lev=Incr.Lev,
                     minusloglLevy = minusloglik,logL.Incr=-esti$value,
                     Data = yuima@data,  yuima=res, Levydetails=esti)
        
        if(length(result@model@jump.coeff)==length(result@model@jump.coeff[[1]])){
          sq<-TRUE
        }else{
          sq<-FALSE
        }
        if(!aggregation){
          result <- vcovLevy1(myres = result,   HessianEta_divTn = HessianEta_divTn,
                Gammahat0 = myGamhat, ures = ures, sq=TRUE)
        }else{
          vcovLevy <- function(myres,   HessianEta_divTn,  Gammahat0, sq=TRUE){
        #myres <- res.VG2
        DeltaX <- apply(myres@Data@original.data,2,diff)
        myData<- myres@Data@original.data
        CmatExpr<- myres@model@jump.coeff
        ncolC <-length(myres@model@jump.coeff[[1]])
        param <- myres@coef[myres@model@parameter@jump]
        aexpr<- myres@model@drift
        namedrift<-myres@model@parameter@drift
        pardrif <- myres@coef[namedrift]
        avect_exp<-myres@model@drift
        
        Jac_Drift <- function(aexpr,namedrift,nobs=length(aexpr)){
          lapply(X=c(1:nobs),FUN=function(X,namedrift,aexpr){
            return(deriv(aexpr[X],namedrift))
          },namedrift=namedrift,aexpr=aexpr)
        }
        
        Jac_DriftExp <- Jac_Drift(aexpr=avect_exp,namedrift)  
        
        FUNDum<-function(foo,myenv) sapply(foo, function(x,env) eval(x,envir = env), env= myenv)
        
        h<-diff(time(myres@Data@zoo.data[[1]]))[1]
        
        dummyf1n<- function(DeltaX, Cmat, h){
          C2<-t(Cmat)%*%Cmat
          dec <- chol(C2) # Eventualy Add a tryCatch
          tmp <- t(DeltaX)%*%solve(C2)%*%DeltaX
          logretval <- -h*sum(log(diag(dec)))  - 0.5 * as.numeric(tmp)
          return(logretval)
        }
        dummyInc<- function(DeltaX,Cmat, avect,h,sq=TRUE){
          if(sq){
            Incr <- solve(Cmat)%*%(DeltaX-avect*h)
          }else{
            C2<-t(Cmat)%*%Cmat
            Incr <- solve(C2)%*%Cmat%*%(DeltaX-avect*h)
          }
          return(Incr)
        }
        
        dummyf2n<- function(DeltaX, avect, Cmat, h){
          C2 <- t(Cmat)%*%Cmat
          Incr <- DeltaX-h*avect
          tmp <- t(Incr)%*%solve(C2)%*%Incr
          logretval <-   - 0.5/h * as.numeric(tmp)
          return(logretval)
        }
        dumGrad_f2n<-function(DeltaX, avect, Jac_a,Cmat, h){
          C2 <- t(Cmat)%*%Cmat
          Incr <- DeltaX-h*avect
          Grad<-t(Jac_a)%*%solve(C2)%*%Incr
          return(Grad)
        }
        
        f1n_j<-function(x, myObsDelta, CmatExpr,h, myData, myres){
          names(myData)<-myres@model@solve.variable
          newenvJumpCoef <- list2env(as.list(c(x,myData)))
          Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenvJumpCoef )
          return(dummyf1n(DeltaX=myObsDelta, Cmat=Cmat, h=h))
        }
        
        Incr_Func <- function(x, myObsDelta, CmatExpr,avect_exp, h, myData, myres, sq=TRUE){
          names(myData)<-myres@model@solve.variable
          newenv <- list2env(as.list(c(x,myData)))
          Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenv )
          dd<-length(avect_exp)
          avect<-numeric(length=dd)
          for(j in c(1:dd)){
            avect[j]<-eval(avect_exp[j],envir=newenv)
          }
          return(dummyInc(DeltaX= myObsDelta,Cmat, avect,h,sq=sq))
        }
        
        
        f2n_j <- function(x, parDiff, myObsDelta, CmatExpr, avect_exp, h, myData, myres){
          names(myData)<-myres@model@solve.variable
          newenvJumpCoef <- list2env(as.list(c(parDiff,myData)))
          Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenvJumpCoef)
          newenvDriftCoef <- list2env(as.list(c(x,myData)))
          dd<-length(avect_exp)
          avect<-numeric(length=dd)
          for(j in c(1:dd)){
            avect[j]<-eval(avect_exp[j],envir=newenvDriftCoef)
          }
          return(dummyf2n(DeltaX=myObsDelta, avect=avect, Cmat=Cmat, h=h))
        }
        
        Gradf2n_j <-function(x, parDiff, myObsDelta, CmatExpr, 
                             avect_exp, Jac_DriftExp, h, myData, myres){
          names(myData)<-myres@model@solve.variable
          newenvJumpCoef <- list2env(as.list(c(parDiff,myData)))
          Cmat<-sapply(X=CmatExpr, FUN=FUNDum, myenv=newenvJumpCoef)
          newenvDriftCoef <- list2env(as.list(c(x,myData)))
          dd<-length(avect_exp)
          avect<-numeric(length=dd)
          Jac_a<- matrix(0,dd,length(x))
          for(j in c(1:dd)){
            avect[j]<-eval(avect_exp[j],envir=newenvDriftCoef)
            Jac_a[j,]<-attr(eval(Jac_DriftExp[[j]],envir=newenvDriftCoef),"gradient")
          }
          return(dumGrad_f2n(DeltaX=myObsDelta, avect=avect, Jac_a=Jac_a, Cmat=Cmat, h=h))
        }
        
        # f1n_j(x=param, myObsDelta=DeltaX[1,], CmatExpr=CmatExpr, h=h)
        #i=1
        # debug(dummyInc)
        # Incr_Func(x=c(param,pardrif), myObsDelta=DeltaX[i,], CmatExpr,avect_exp, h, myData=myData[i,], myres, sq=TRUE)
        # f1n_j(x=param, myObsDelta=DeltaX[i,], CmatExpr=CmatExpr, h=h,  myData=myData[i,], myres=myres)
        # f2n_j(x, parDiff=param, myObsDelta=DeltaX[i,], CmatExpr=CmatExpr, avect_exp=avect_exp, 
        #       h=h, myData=myData[i,], myres=myres)
        # Gradf2n_j(x, parDiff=param, myObsDelta=DeltaX[i,], CmatExpr, 
        #     avect_exp, Jac_DriftExp, h, myData[i,], myres)
        del <- 10^-3
        dummy <- t(rep(1,length(param)))
        myeta<- as.matrix(param)%*%dummy
        dummyG <- t(rep(1,length(c(param,pardrif))))
        globalmyeta <- as.matrix(c(param,pardrif))%*%dummyG
        Incr <- del*diag(rep(1,dim(myeta)[1]))
        GlobIncr <- del*diag(rep(1,dim(globalmyeta)[1]))
        Sigma_gamma <- matrix(0 ,length(param),length(param))
        Sigma_alpha <- matrix(0 ,length(pardrif),length(pardrif))
        Sigma_algam <- matrix(0 ,length(pardrif),length(param))
        histDeltaf1n_j_delta <-matrix(0 , dim(DeltaX)[1],length(param))
        histDeltaf2n_j_delta <- matrix(0 , dim(DeltaX)[1],length(pardrif))
        histb_incr <- array(0, c(dim(DeltaX)[2],dim(DeltaX)[1],length(c(param,pardrif))))
        for(i in c(1:dim(DeltaX)[1])){
          Deltaf1n_j_delta<-sapply(X=1:dim(myeta)[1],
                                   FUN = function(X,myObsDelta, CmatExpr, h, myData, myres,del){
                                     par1<-myeta[,X]+Incr[,X]
                                     #par[oldyuima@model@measure$df@time.var]<-1
                                     f1<-f1n_j(x=par1, myObsDelta, CmatExpr, h, myData, myres)
                                     par2<-myeta[,X]-Incr[,X]
                                     f2<-f1n_j(x=par2, myObsDelta, CmatExpr, h, myData, myres)
                                     return((f1-f2)/(2*del))
                                   },
                                   myObsDelta=DeltaX[i,], CmatExpr=CmatExpr, h=h, myData=myData[i,],
                                   myres=myres, del=del)
          
          DeltaInc <- sapply(X=1:dim(GlobIncr)[1],
                             FUN = function(X, myObsDelta, CmatExpr,avect_exp, h, 
                                            myData, myres, del, sq){
                               par1<-globalmyeta[,X]+GlobIncr[,X]
                               
                               # f1<-f1n_j(x=par1, myObsDelta, CmatExpr, h, myData, myres)
                               f1<-Incr_Func(x=par1, myObsDelta, CmatExpr,avect_exp, h, myData, myres, sq)
                               par2<-globalmyeta[,X]-GlobIncr[,X]
                               f2<-Incr_Func(x=par2, myObsDelta, CmatExpr,avect_exp, h, myData, myres, sq)
                               #f2<-f1n_j(x=par2, myObsDelta, CmatExpr, h, myData, myres)
                               return((f1-f2)/(2*del))
                             },
                             myObsDelta=DeltaX[i,], CmatExpr,avect_exp, h, myData=myData[i,], myres, del=del,sq=sq
          )
          
          histDeltaf1n_j_delta[i, ]<-Deltaf1n_j_delta
          
          histb_incr[ , i, ] <- DeltaInc
          Sigma_gamma <- Sigma_gamma + as.matrix(Deltaf1n_j_delta)%*%t(Deltaf1n_j_delta)
          
          Deltaf2n_j_delta<- Gradf2n_j(x=pardrif, parDiff=param, myObsDelta=DeltaX[i,], CmatExpr, 
                                       avect_exp, Jac_DriftExp, h, myData[i,], myres)
          histDeltaf2n_j_delta[i,]<-Deltaf2n_j_delta
          
          Sigma_alpha <- Sigma_alpha + Deltaf2n_j_delta%*%t(Deltaf2n_j_delta)
          Sigma_algam <- Sigma_algam + Deltaf2n_j_delta%*%t(Deltaf1n_j_delta)
          #cat("\n",i)
        }
        Tn<-tail(time(myres@Data@original.data),1L)
        Sigma_gamma0<-Sigma_gamma/Tn
        Sigma_alpha0 <- Sigma_alpha/Tn     
        Sigma_algam0 <- Sigma_algam/Tn
        Sigma0 <- cbind(rbind(Sigma_gamma0,Sigma_algam0),rbind(t(Sigma_algam0),Sigma_alpha0))
        
        myLogDens <- function(x, myObsDelta, CmatExpr,h,myData, myres, pos){
          return(f1n_j(x, myObsDelta[pos, ], CmatExpr,h,myData=myData[pos,], myres=myres))
        }
        
        myLogDens2<- function(x, parDiff, myObsDelta, CmatExpr, avect_exp, h, myData, myres, pos){
          return(f2n_j(x, parDiff, myObsDelta[pos, ], CmatExpr, avect_exp, h, myData[pos,], myres))
          #return(f1n_j(x, myObsDelta[pos, ], CmatExpr,h,myData=myData[pos,], myres=myres))
        }
        
        
        VectLogDens<-Vectorize(FUN=myLogDens,vectorize.args = "pos")
        #VectLogDens2<-Vectorize(FUN=myLogDens2,vectorize.args = "pos")
        
        mylogLik<- function(par,DeltaX,CmatExpr, myData, myres, h, Tn, nobs=dim(DeltaX)[1]){
          term<-VectLogDens(x=par, 
                            myObsDelta=DeltaX, CmatExpr=CmatExpr,h=h, myData=myData, myres=myres,pos=c(1:nobs))
          return(sum(term)/Tn)
        }
        
        mylogLik2<- function(par, parDiff, myObsDelta, CmatExpr, avect_exp, h, Tn, myData, myres, nobs=dim(DeltaX)[1]){
          # term<-VectLogDens2(x=par, parDiff, myObsDelta, CmatExpr, avect_exp, h, myData, myres,
          #                    pos=c(1:nobs))
          term<-0
          for(j in  c(1:nobs)){
            term<-term+f2n_j(x=par, parDiff, myObsDelta[j, ], CmatExpr, avect_exp, h, myData[j,], myres)
          }
          return(term/Tn)
        }
        
        GradlogLik2<-function(par, parDiff, myObsDelta, CmatExpr, avect_exp, h, Tn, myData, myres, nobs=dim(DeltaX)[1]){
          term<-matrix(0, length(par),1)
          for(j in  c(1:nobs)){ 
            term <- term + Gradf2n_j(x=par, parDiff, myObsDelta=myObsDelta[j, ], CmatExpr, 
                                     avect_exp, Jac_DriftExp, h, myData=myData[j,], myres)
          }
          return(as.numeric(term)/Tn)
        }
        
        # Gammahat <- optimHess(par=param,fn=mylogLik,
        #                       DeltaX=DeltaX,CmatExpr=CmatExpr,
        #                       h=h, Tn=Tn, nobs=dim(DeltaX)[1], 
        #                       myData=myData, myres=myres)
        
        # GammaAlfahat <- optimHess(par=pardrif,fn=mylogLik2,gr=GradlogLik2,
        #                           parDiff=param, myObsDelta=DeltaX, 
        #                           CmatExpr=CmatExpr, avect_exp=avect_exp, h=h,
        #                           Tn=Tn, nobs=dim(DeltaX)[1], myData=myData, myres=myres)
        # 
        # Gammahat0 <- cbind(rbind(Gammahat,matrix(0,dim(GammaAlfahat)[1],dim(Gammahat)[2])),
        #                    rbind(matrix(0,dim(Gammahat)[2],dim(GammaAlfahat)[1]),GammaAlfahat))
        
        uhistDeltaf1n_j<- matrix(0,floor(Tn),length(param))
        uhistDeltaf2n_j<- matrix(0,floor(Tn),length(pardrif))
        aaa<-length(c(param,pardrif))
        uhistb_incr <- array(0, c(dim(DeltaX)[2],floor(Tn),aaa))
        for(l in 1:floor(Tn)){
          #ures[l] <- sum(resi[(floor((l-1)/h)):(floor(l/h)-1)]) 
          uhistDeltaf1n_j[l, ] <- colSums(histDeltaf1n_j_delta[(floor((l-1)/h)):(floor(l/h)-1),])
          uhistDeltaf2n_j[l, ] <- colSums(histDeltaf2n_j_delta[(floor((l-1)/h)):(floor(l/h)-1),])
          for(i in c(1:aaa)){
            dumIn <- histb_incr[,(floor((l-1)/h)):(floor(l/h)-1),i]
            uhistb_incr[,l,i]<- rowSums(dumIn)
          }
        }
        
        nMeaspar<-myres@model@parameter@measure
        mypar_meas0 <- myres@coef[nMeaspar]#res@coef[oldyuima@model@parameter@measure]
        mypar_meas <- c(mypar_meas0,1)
        ures<-myres@Incr.Lev@original.data
        names(mypar_meas)<-c(nMeaspar,myres@model@measure$df@time.var)
        fdataeta<-dens(object=myres@model@measure$df,x=ures,param=mypar_meas)# f(eps, eta)
        del <- 10^-3
        dummy <- t(rep(1,length(mypar_meas[myres@model@measure$df@param.measure])))
        myeta<- as.matrix(mypar_meas[myres@model@measure$df@param.measure])%*%dummy
        myetapert <- myeta+del*diag(rep(1,dim(myeta)[1])) 
        fdataetadelta<-sapply(X=1:dim(myeta)[1],FUN = function(X){
          par<-myetapert[,X]
          par[myres@model@measure$df@time.var]<-1
          dens(object=myres@model@measure$df,x=ures,param=par)
        }
        )# f(eps, eta+delta)
        
        term1<-1/(fdataeta)
        
        dummyvecdel<-numeric(length=dim(ures)[2]) 
        fdatadeltaeta<- matrix(0,floor(Tn),dim(ures)[2])
        mixpartial <-array(0,c(length(nMeaspar),floor(Tn),dim(ures)[2]))
        for(j in c(1:dim(ures)[2])){
          dummyvecdel[j]<-del
          fdatadeltaeta[,j]<- dens(object=myres@model@measure$df,x=ures+rep(1,dim(ures)[1])%*%t(dummyvecdel),param=mypar_meas)
          
          fdatadeltaetadelta<-sapply(X=1:dim(myeta)[1],FUN = function(X){
            par<-myetapert[,X]
            par[myres@model@measure$df@time.var]<-1
            dens(object=myres@model@measure$df,x=ures+rep(1,dim(ures)[1])%*%t(dummyvecdel),param=par)
          }
          ) # f(eps +deta, eta+delta)
          dummyvecdel<-numeric(length=dim(ures)[2])
          term2 <- fdatadeltaeta[,j]*term1%*%dummy
          term2 <- fdatadeltaetadelta-term2*fdataetadelta
          mixpartial[,,j]<-t((as.matrix(term1)%*%dummy)/del^2*term2/Tn)
        }
        
        DerMeta <- 1/del*(fdataetadelta - as.matrix(fdataeta)%*%rep(1,dim(fdataetadelta)[2]))*(as.matrix(term1)%*%rep(1,dim(fdataetadelta)[2]))
        SigmaEta <- t(DerMeta)%*%DerMeta/Tn
        
        minusloglik <- function(para){
          para[length(para)+1]<-1
          names(para)[length(para)]<-myres@model@time.variable
          -sum(dens(object=myres@model@measure$df, x=ures, param = para, log = TRUE), 
               na.rm = T)/Tn
        }
        #  HessianEta_divTn <- optimHess(par=mypar_meas0, fn=minusloglik)
        
        Gammaeta_theta <- matrix(0,length(mypar_meas0), aaa)
        for(t in c(1:floor(Tn))){
          Gammaeta_theta<-Gammaeta_theta+mixpartial[,t,]%*%uhistb_incr[,t,]
        }
        GammaEta_Theta<- -1/Tn*Gammaeta_theta
        
        Sigma_eta_theta<-t(DerMeta)%*%cbind(uhistDeltaf1n_j,uhistDeltaf2n_j)/Tn
        
        Sigma <- cbind(rbind(Sigma0,Sigma_eta_theta),rbind(t(Sigma_eta_theta),SigmaEta))
        
        GammaHAT<-cbind(rbind(Gammahat0, GammaEta_Theta), rbind(t(GammaEta_Theta), HessianEta_divTn))
        InvGammaHAT <- solve(GammaHAT)
        vcov <- InvGammaHAT %*% Sigma %*% InvGammaHAT/Tn
        
        myres@vcov<- vcov
        return(myres)
          }
          result <- vcovLevy(result,HessianEta_divTn,myGamhat,sq=sq)
        }
        
      
        #result <- vcovLevy(result,HessianEta_divTn,myGamhat,sq=sq)
      # colnames(result@vcov)<-names(res@fullcoef)
      # rownames(result@vcov)<-names(res@fullcoef)
    }
    # colnames(res@vcov)<-names(res@fullcoef)
    # #res@vcov<-rbind(res@vcov,matrix(NA,nrow=length(esti$par),ncol=dim(res@vcov)[2]))
    # rownames(res@vcov)<-names(res@fullcoef)
    # res@min<-c(res@min,esti$value)
    # res@nobs<-c(res@nobs,length(Incr.Lev@zoo.data[[1]]))
    
    
    
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


