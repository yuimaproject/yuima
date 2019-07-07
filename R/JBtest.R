JBtest<-function(yuima,start,lower,upper,alpha,skewness=TRUE,kurtosis=TRUE,withdrift=FALSE){
  
  ## The new options "skewness" and "kurtosis" are added (3/27).
  ## Multivariate Mardia's skewness is tentativity written, not work now (5/30), workable (6/6)!.
  
  if(yuima@model@equation.number == 1){
    data <- get.zoo.data(yuima)
    s.size<-yuima@sampling@n 
    modelstate<-yuima@model@solve.variable
    modeltime<-yuima@model@time.variable
    DRIFT<-yuima@model@drift
    DIFFUSION<-yuima@model@diffusion
    PARLENGS<-length(yuima@model@parameter@drift)+length(yuima@model@parameter@diffusion)
    XINIT<-yuima@model@xinit
    X<-as.numeric(data[[1]])
    dX <- abs(diff(X))
    odX <- sort(dX, decreasing=TRUE) 
    plot(yuima,main="Original path")
    abline(0,0,lty=5)
    pX<-X[1:(s.size-1)]
    sY<-as.numeric(data[[1]])
    derDIFFUSION<-D(DIFFUSION[[1]],modelstate)
    tmp.env<-new.env()
    INTENSITY<-eval(yuima@model@measure$intensity)
    
    inc<-double(s.size-1)
    inc<-X[2:(s.size)]-pX
    preservedinc<-inc
    ainc<-abs(inc)
    oinc<-sort(ainc,decreasing=TRUE)
    if(length(yuima@model@parameter@measure)!=0){
      extp <- match(yuima@model@parameter@measure, yuima@model@parameter@all)
      extplow <- match(yuima@model@parameter@measure, names(lower))
      extpupp <- match(yuima@model@parameter@measure, names(upper))
      extpsta <- match(yuima@model@parameter@measure, names(start))
      lower <- lower[-extplow]
      upper <- upper[-extpupp]
      start <- start[-extpsta]
      yuima@model@parameter@all <- yuima@model@parameter@all[-extp]
      yuima@model@parameter@measure <- as.character("abcdefg")
      yuima@model@measure$df$expr <- expression(dunif(z,abcdefg,10))
      lower <- c(lower,abcdefg=2-10^(-2))
      upper <- c(upper,abcdefg=2+10^(-2))
      start <- c(start,abcdefg=2)
      yuima@model@parameter@all <- c(yuima@model@parameter@all,as.character("abcdefg"))
    }else{
      yuima@model@parameter@measure <- as.character("abcdefg")
      yuima@model@measure$df$expr <- expression(dunif(z,abcdefg,10))
      lower <- c(lower,abcdefg=2-10^(-2))
      upper <- c(upper,abcdefg=2+10^(-2))
      start <- c(start,abcdefg=2)
      yuima@model@parameter@all <- c(yuima@model@parameter@all,as.character("abcdefg"))
    }
    #yuima@model@drift<-expression((0))
    qmle<-qmle(yuima, start = start, lower = lower, upper = upper,threshold = oinc[1]+1) # initial estimation
    parameter<-yuima@model@parameter@all
    mp<-match(names(qmle@coef),parameter)
    esort <- qmle@coef[order(mp)]
    
    for(i in 1:length(parameter))
    {
      assign(parameter[i],esort[[i]],envir=tmp.env)
    }
    
    resi<-double(s.size-1) 
    derv<-double(s.size-1)
    total<-c()
    j.size<-c()
    
    assign(modeltime,yuima@sampling@delta,envir=tmp.env)
    h<-yuima@sampling@delta
    
    assign(modelstate,pX,envir=tmp.env)
    diff.term<-eval(DIFFUSION[[1]],envir=tmp.env)
    drif.term<-eval(DRIFT,envir=tmp.env)
    if(length(diff.term)==1){
      diff.term <- rep(diff.term, s.size)
    }
    if(length(drif.term)==1){
      drif.term <- rep(drif.term, s.size)
    } # vectorization (note. if an expression type object does not include state.variable, the length of the item after "eval" operation is 1.)
    for(s in 1:(s.size-1)){
      nova<-sqrt((diff.term)^2) # normalized variance
      resi[s]<-(1/(nova[s]*sqrt(h)))*(inc[s]-h*withdrift*drif.term[s])
    }    
    assign(modelstate,pX,envir=tmp.env)
    derv<-eval(derDIFFUSION,envir=tmp.env)
    if(length(derv)==1){
      derv<-rep(derv,s.size) # vectorization
    }
    
    mresi<-mean(resi)
    vresi<-mean((resi-mresi)^2)
    snr<-(resi-mresi)/sqrt(vresi)
    isnr<-snr
    
    if(skewness==TRUE&kurtosis==FALSE){
      bias<-3*sqrt(h)*sum(derv)
      JB<-1/(6*s.size)*(sum(snr^3)-bias)^2
      rqmle<-list() # jump removed qmle
      statis<-c() # the value of JB statistics
      limJB<-5*INTENSITY*yuima@sampling@Terminal
      # the limit of JB repetition 
      klarinc<-numeric(floor(limJB)) 
      k<-1
      if(JB<=qchisq(alpha,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)){
        print("There is no jump.")
        result<-list(OGQMLE=qmle@coef[1:PARLENGS])
      }else{
        while(JB>qchisq(alpha,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE))
        { 
          klarinc[k]<-which(ainc==oinc[k]) ## detect the largest increment
          j.size<-append(j.size,round(preservedinc[klarinc[k]],digits=3))
          rqmle[[k]]<-qmle(yuima, start = start, lower = lower, upper = upper, threshold = oinc[k]-(oinc[k]-oinc[k+1])/2)@coef
          ## calculate the modified qmle
          mp<-match(names(rqmle[[k]]),parameter)
          resort <- rqmle[[k]][order(mp)]
          
          for(i in 1:length(parameter))
          {
            assign(parameter[i],resort[[i]],envir=tmp.env)
          }
          
          diff.term<-eval(DIFFUSION[[1]],envir=tmp.env)
          drif.term<-eval(DRIFT,envir=tmp.env)
          if(length(diff.term)==1){
            diff.term <- rep(diff.term, s.size-1)
          }
          if(length(drif.term)==1){
            drif.term <- rep(drif.term, s.size-1)
          } # vectorization 
          for(s in 1:(s.size-1)){
            nova<-sqrt((diff.term)^2) # normalized variance
            resi[s]<-(1/(nova[s]*sqrt(h)))*(inc[s]-h*withdrift*drif.term[s])
          }          ## rebuild Euler-residuals
          
          derv<-eval(derDIFFUSION,envir=tmp.env)
          if(length(derv)==1){
            derv<-rep(derv,s.size) # vectorization
          }
          ## rebuild derivative vector
          
          resi[klarinc]<-0
          derv[klarinc]<-0
          mresi<-mean(resi)
          vresi<-mean((resi-mresi)^2)
          snr<-(resi-mresi)/sqrt(vresi) ## rebuild self-normalized residuals
          snr[klarinc]<-0
          
          bias<-3*sqrt(h)*sum(derv)
          JB<-1/(6*(s.size-k))*(sum(snr^3)-bias)^2
          jump.point<-klarinc[k]
          points(jump.point*h,sY[klarinc[k]],col="red", type="h", lty=3, lwd=2)
          points(jump.point*h,sY[klarinc[k]],col="red", pch="+",lwd=2)
          points(jump.point*h,0,col="red", pch=17,cex=1.5)
          total<-append(total,as.character(round(jump.point*h,digits=3)))
          statis<-append(statis,JB)
          k<-k+1
          if(k == floor(limJB))
          { 
            warning("Removed jumps seems to be too many. Change intensity for more removal.")
            break
          }
        }
        par(mfrow=c(2,2)) 
        # plot(isnr,main="Raw self-normalized residual")
        par(new=T)
        abline(0,0,lty=5,col="green")
        snr <- snr[-(klarinc)]
        # h <- dpih(snr)      
        # bins <- seq(min(snr)-0.1, max(snr)+0.1+h, by=h)  
        # hist(snr, prob=1, breaks=bins,xlim=c(-3,3),ylim=c(0,1))
        hist(snr, freq = FALSE, breaks = "freedman-diaconis", xlim=c(min(snr),max(snr)),ylim=c(0,1))
        xval <- seq(min(snr),max(snr),0.01) # to add plot the corresponding density
        lines(xval,dnorm(xval,0,1),xlim=c(-3,3),ylim=c(0,1),col="red",main="Jump removed self-normalized residual")
        if(k > 2){
          estimat<-matrix(0,PARLENGS,k-1)
          for(i in 1:(k-1)){
            for(j in 1:PARLENGS){
              estimat[j,i]<-rqmle[[i]][j]
            }
          }
          oo<-match(parameter,names(estimat[1,]))
          parameter<-parameter[order(oo)]
          for(i in 1:PARLENGS){
            plot(estimat[i,], main=paste("Transition of",parameter[i]),
                 xlab="The number of trial",ylab = "Estimated value",type="l",xlim = c(1,k-1),xaxp = c(1,k-1,k-2))
          }
          plot(log(statis),main="Transition of the logJB statistics",
               xlab="The number of trial",ylab="log JB"
               ,type="l",xlim = c(1,k-1),xaxp = c(1,k-1,k-2))
          par(new=T)
          abline(log(qchisq(alpha,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)),0,lty=5,col="red")
        }
        names(j.size)<-total
        removed<-j.size
        resname<-c("Jump size")
        names(resname)<-"Jump time"
        removed<-append(resname,removed)
        result<-list(Removed=removed,OGQMLE=qmle@coef[1:PARLENGS],JRGQMLE=rqmle[[k-1]][1:PARLENGS])
      }
    }else if(skewness==FALSE&kurtosis==TRUE){
      JB<-1/(24*s.size)*(sum(snr^4-3))^2
      fbias<-rep(3,s.size-1)
      rqmle<-list() # jump removed qmle
      statis<-c() # the value of JB statistics
      limJB<-5*INTENSITY*yuima@sampling@Terminal
      # the limit of JB repetition 
      klarinc<-numeric(floor(limJB)) 
      k<-1
      
      if(JB<=qchisq(alpha,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)){
        print("There is no jump.")
        result<-list(OGQMLE=qmle@coef[1:PARLENGS])
      }else{
        while(JB>qchisq(alpha,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE))
        { 
          klarinc[k]<-which(ainc==oinc[k]) ## detect the largest increment
          j.size<-append(j.size,round(preservedinc[klarinc[k]],digits=3))
          rqmle[[k]]<-qmle(yuima, start = start, lower = lower, upper = upper, threshold = oinc[k]-(oinc[k]-oinc[k+1])/2)@coef
          ## calculate the modified qmle
          mp<-match(names(rqmle[[k]]),parameter)
          resort <- rqmle[[k]][order(mp)]
          
          for(i in 1:length(parameter))
          {
            assign(parameter[i],resort[[i]],envir=tmp.env)
          }
          
          diff.term<-eval(DIFFUSION[[1]],envir=tmp.env)
          drif.term<-eval(DRIFT,envir=tmp.env)
          if(length(diff.term)==1){
            diff.term <- rep(diff.term, s.size-1)
          }
          if(length(drif.term)==1){
            drif.term <- rep(drif.term, s.size-1)
          } # vectorization 
          for(s in 1:(s.size-1)){
            nova<-sqrt((diff.term)^2) # normalized variance
            resi[s]<-(1/(nova[s]*sqrt(h)))*(inc[s]-h*withdrift*drif.term[s])
          }          ## rebuild Euler-residuals
          
          derv<-eval(derDIFFUSION,envir=tmp.env)
          if(length(derv)==1){
            derv<-rep(derv,s.size) # vectorization
          }
          
          resi[klarinc]<-0
          derv[klarinc]<-0
          mresi<-mean(resi)
          vresi<-mean((resi-mresi)^2)
          snr<-(resi-mresi)/sqrt(vresi) ## rebuild self-normalized residuals
          snr[klarinc]<-0
          fbias[klarinc]<-0
          
          JB<-1/(24*(s.size-k))*(sum(snr^4-fbias))^2
          jump.point<-klarinc[k]
          points(jump.point*h,sY[klarinc[k]],col="red", type="h", lty=3, lwd=2)
          points(jump.point*h,sY[klarinc[k]],col="red", pch="+",lwd=2)
          points(jump.point*h,0,col="red", pch=17,cex=1.5)
          total<-append(total,as.character(round(jump.point*h,digits=3)))
          statis<-append(statis,JB)
          k<-k+1
          if(k == floor(limJB))
          { 
            warning("Removed jumps seems to be too many. Change intensity for more removal.")
            break
          }
        }
        par(mfrow=c(2,2)) 
        # plot(isnr,main="Raw self-normalized residual")
        par(new=T)
        abline(0,0,lty=5,col="green")
        snr <- snr[-(klarinc)]
        # h <- dpih(snr)      
        # bins <- seq(min(snr)-0.1, max(snr)+0.1+h, by=h)  
        # hist(snr, prob=1, breaks=bins,xlim=c(-3,3),ylim=c(0,1))
        hist(snr, freq = FALSE, breaks = "freedman-diaconis", xlim=c(min(snr),max(snr)),ylim=c(0,1))
        xval <- seq(min(snr),max(snr),0.01) # to add plot the corresponding density
        lines(xval,dnorm(xval,0,1),xlim=c(-3,3),ylim=c(0,1),col="red",main="Jump removed self-normalized residual")
        if(k > 2){
          estimat<-matrix(0,PARLENGS,k-1)
          for(i in 1:(k-1)){
            for(j in 1:PARLENGS){
              estimat[j,i]<-rqmle[[i]][j]
            }
          }
          oo<-match(parameter,names(estimat[1,]))
          parameter<-parameter[order(oo)]
          for(i in 1:PARLENGS){
            plot(estimat[i,], main=paste("Transition of",parameter[i]),
                 xlab="The number of trial",ylab = "Estimated value",type="l",xlim = c(1,k-1),xaxp = c(1,k-1,k-2))
          }
          plot(log(statis),main="Transition of the logJB statistics",
               xlab="The number of trial",ylab="log JB"
               ,type="l",xlim = c(1,k-1),xaxp = c(1,k-1,k-2))
          par(new=T)
          abline(log(qchisq(alpha,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)),0,lty=5,col="red")
        }
        names(j.size)<-total
        removed<-j.size
        resname<-c("Jump size")
        names(resname)<-"Jump time"
        removed<-append(resname,removed)
        result<-list(Removed=removed,OGQMLE=qmle@coef[1:PARLENGS],JRGQMLE=rqmle[[k-1]][1:PARLENGS])
      }
    }else{
      fbias<-rep(3,s.size-1)
      bias<-3*sqrt(h)*sum(derv)
      JB<-1/(6*s.size)*(sum(snr^3)-bias)^2+1/(24*s.size)*(sum(snr^4-fbias))^2
      rqmle<-list() # jump removed qmle
      statis<-c() # the value of JB statistics
      limJB<-5*INTENSITY*yuima@sampling@Terminal
      # the limit of JB repetition 
      klarinc<-numeric(floor(limJB)) 
      k<-1
      
      
      if(JB<=qchisq(alpha,df=2,ncp=0,lower.tail=FALSE,log.p=FALSE)){
        print("There is no jump.")
        result<-list(OGQMLE=qmle@coef[1:PARLENGS])
      }else{
        while(JB>qchisq(alpha,df=2,ncp=0,lower.tail=FALSE,log.p=FALSE))
        { 
          klarinc[k]<-which(ainc==oinc[k]) ## detect the largest increment
          j.size<-append(j.size,round(preservedinc[klarinc[k]],digits=3))
          rqmle[[k]]<-qmle(yuima, start = start, lower = lower, upper = upper, threshold = oinc[k]-(oinc[k]-oinc[k+1])/2)@coef
          ## calculate the modified qmle
          mp<-match(names(rqmle[[k]]),parameter)
          resort <- rqmle[[k]][order(mp)]
          
          for(i in 1:length(parameter))
          {
            assign(parameter[i],resort[[i]],envir=tmp.env)
          }
          
          diff.term<-eval(DIFFUSION[[1]],envir=tmp.env)
          drif.term<-eval(DRIFT,envir=tmp.env)
          if(length(diff.term)==1){
            diff.term <- rep(diff.term, s.size-1)
          }
          if(length(drif.term)==1){
            drif.term <- rep(drif.term, s.size-1)
          } # vectorization 
          for(s in 1:(s.size-1)){
            nova<-sqrt((diff.term)^2) # normalized variance
            resi[s]<-(1/(nova[s]*sqrt(h)))*(inc[s]-h*withdrift*drif.term[s])
          }          ## rebuild Euler-residuals
          
          
          derv<-eval(derDIFFUSION,envir=tmp.env)
          if(length(derv)==1){
            derv<-rep(derv,s.size) # vectorization
          }
          ## rebuild derivative vector
          
          resi[klarinc]<-0
          derv[klarinc]<-0
          mresi<-mean(resi)
          vresi<-mean((resi-mresi)^2)
          snr<-(resi-mresi)/sqrt(vresi) ## rebuild self-normalized residuals
          snr[klarinc]<-0
          
          bias<-3*sqrt(h)*sum(derv)
          fbias[klarinc]<-0
          JB<-1/(6*(s.size-k))*(sum(snr^3)-bias)^2+1/(24*(s.size-k))*(sum(snr^4-fbias))^2
          jump.point<-klarinc[k]
          points(jump.point*h,sY[klarinc[k]],col="red", type="h", lty=3, lwd=2)
          points(jump.point*h,sY[klarinc[k]],col="red", pch="+",lwd=2)
          points(jump.point*h,0,col="red", pch=17,cex=1.5)
          total<-append(total,as.character(round(jump.point*h,digits=3)))
          statis<-append(statis,JB)
          k<-k+1
          if(k == floor(limJB))
          { 
            warning("Removed jumps seems to be too many. Change intensity for more removal.")
            break
          }
        }
        par(mfrow=c(2,2)) 
        # plot(isnr,main="Raw snr")
        par(new=T)
        abline(0,0,lty=5,col="green")
        snr <- snr[-(klarinc)]
        # h <- dpih(snr)      
        # bins <- seq(min(snr)-0.1, max(snr)+0.1+h, by=h)  
        # hist(snr, prob=1, breaks=bins,xlim=c(-3,3),ylim=c(0,1))
        hist(snr, freq = FALSE, breaks = "freedman-diaconis", xlim=c(min(snr),max(snr)),ylim=c(0,1))
        xval <- seq(min(snr),max(snr),0.01) # to add plot the corresponding density
        lines(xval,dnorm(xval,0,1),xlim=c(-3,3),ylim=c(0,1),col="red",main="Jump removed self-normalized residual")
        if(k > 2){
          estimat<-matrix(0,PARLENGS,k-1)
          for(i in 1:(k-1)){
            for(j in 1:PARLENGS){
              estimat[j,i]<-rqmle[[i]][j]
            }
          }
          oo<-match(parameter,names(estimat[1,]))
          parameter<-parameter[order(oo)]
          for(i in 1:PARLENGS){
            plot(estimat[i,], main=paste("Transition of",parameter[i]),
                 xlab="The number of trial",ylab = "Estimated value",type="l",xlim = c(1,k-1),xaxp = c(1,k-1,k-2))
          }
          plot(log(statis),main="Transition of the logJB statistics",
               xlab="The number of trial",ylab="log JB"
               ,type="l",xlim = c(1,k-1),xaxp = c(1,k-1,k-2))
          par(new=T)
          abline(log(qchisq(alpha,df=2,ncp=0,lower.tail=FALSE,log.p=FALSE)),0,lty=5,col="red")
        }
        names(j.size)<-total
        removed<-j.size
        resname<-c("Jump size")
        names(resname)<-"Jump time"
        removed<-append(resname,removed)
        result<-list(Removed=removed,OGQMLE=qmle@coef[1:PARLENGS],JRGQMLE=rqmle[[k-1]][1:PARLENGS])
      }
    }
    # }else{
    #   ## Multivariate version based on Mardia's kurtosis, not work now !
    #   DRIFT <- yuima@model@drift
    #   DIFFUSION <- yuima@model@diffusion
    #   d.size <- yuima@model@equation.number
    #   data <- matrix(0,length(yuima@data@zoo.data[[1]]),d.size)
    #   for(i in 1:d.size) data[,i] <- as.numeric(yuima@data@zoo.data[[i]])
    #   dx_set <- as.matrix((data-rbind(numeric(d.size),as.matrix(data[-length(data[,1]),])))[-1,])
    #   preservedinc <- dx_set
    #   ainc <- numeric(length(yuima@data@zoo.data[[1]])-1)
    #   for(i in 1:length(yuima@data@zoo.data[[1]])-1){
    #     ainc[i] <- sqrt(t(dx_set[i,])%*%dx_set[i,])
    #   }
    #   oinc<-sort(ainc,decreasing=TRUE)
    #   qmle<-qmle(yuima, start = start, lower = lower, upper = upper,threshold = oinc[1]+1) # initial estimation
    #   
    #   parameter<-yuima@model@parameter@all
    #   mp<-match(names(qmle@coef),parameter)
    #   esort <- qmle@coef[order(mp)]
    #   
    #   tmp.env <- new.env()
    #   for(i in 1:length(parameter))
    #   {
    #     assign(parameter[i],esort[[i]],envir=tmp.env)
    #   }
    #   
    #   noise_number <- yuima@model@noise.number
    #   
    #   for(i in 1:d.size) assign(yuima@model@state.variable[i], data[-length(data[,1]),i], envir=tmp.env)
    #   
    #   d_b <- NULL
    #   for(i in 1:d.size){
    #     if(length(eval(DRIFT[[i]],envir=tmp.env))==(length(data[,1])-1)){
    #       d_b[[i]] <- DRIFT[[i]] #this part of model includes "x"(state.variable)
    #     }
    #     else{
    #       if(is.na(c(DRIFT[[i]][2]))){ #ex. yuima@model@drift=expression(0) (we hope "expression((0))")
    #         DRIFT[[i]] <- parse(text=paste(sprintf("(%s)", DRIFT[[i]])))[[1]]
    #       }
    #       d_b[[i]] <- parse(text=paste("(",DRIFT[[i]][2],")*rep(1,length(data[,1])-1)",sep=""))
    #       #vectorization
    #     }
    #   }
    #   
    #   v_a<-matrix(list(NULL),d.size,noise_number)
    #   for(i in 1:d.size){
    #     for(j in 1:noise_number){
    #       if(length(eval(DIFFUSION[[i]][[j]],envir=tmp.env))==(length(data[,1])-1)){
    #         v_a[[i,j]] <- DIFFUSION[[i]][[j]] #this part of model includes "x"(state.variable)
    #       }
    #       else{
    #         if(is.na(c(DIFFUSION[[i]][[j]][2]))){
    #           DIFFUSION[[i]][[j]] <- parse(text=paste(sprintf("(%s)", DIFFUSION[[i]][[j]])))[[1]]
    #         }
    #         v_a[[i,j]] <- parse(text=paste("(",DIFFUSION[[i]][[j]][2],")*rep(1,length(data[,1])-1)",sep=""))
    #         #vectorization
    #       }
    #     }
    #   }
    #   
    #   resi<-NULL
    #   nvari<-matrix(0,d.size,noise_number)
    #   drivec<-numeric(d.size)
    #   for(k in 1:(length(data[,1])-1)){
    #     for(i in 1:d.size)
    #     {
    #       for(j in 1:noise_number){
    #         nvari[i,j]<-eval(v_a[[i,j]],envir=tmp.env)[k]
    #       }
    #       drivec[i]<-eval(d_b[[i]],envir=tmp.env)[k]
    #     }
    #     resi[[k]]<-1/sqrt(yuima@sampling@delta[[1]])*solve(t(nvari)%*%nvari)%*%t(nvari)%*%(dx_set[k,]-withdrift*drivec)
    #   }
    #   
    #   mresi<-numeric(d.size)
    #   for(k in 1:(length(data[,1])-1)){
    #     mresi<-mresi+resi[[k]]
    #   }
    #   mresi<-mresi/(length(data[,1])-1)
    #   
    #   svresi<-matrix(0,d.size,d.size)
    #   for(k in 1:(length(data[,1])-1)){
    #     svresi<-svresi+(resi[[k]]-mresi)%*%t(resi[[k]]-mresi)
    #   }
    #   svresi<-svresi/(length(data[,1])-1)
    #   
    #   snr<-NULL
    #   invsvresi<-solve(svresi)
    #   U <- svd(invsvresi)$u
    #   V <- svd(invsvresi)$v
    #   D <- diag(sqrt(svd(invsvresi)$d))
    #   sqinvsvresi <- U %*% D %*% t(V) 
    #   for(k in 1:(length(data[,1])-1)){
    #     snr[[k]]<-sqinvsvresi%*%(resi[[k]]-mresi)
    #   }
    #   fourmom<-numeric(1)
    #   for(k in 1:(length(data[,1])-1)){
    #     fourmom<-fourmom+(t(snr[[k]])%*%snr[[k]])^2
    #   }
    #   
    #   Mard<-1/((length(data[,1])-1)*8*d.size*(d.size+2))*(fourmom-(length(data[,1])-1)*d.size*(d.size+2))^2
    #   if(Mard<=qchisq(alpha,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)){
    #     print("There is no jump.")
    #   }else{ 
    #     INTENSITY<-eval(yuima@model@measure$intensity)
    #     limMard<-5*INTENSITY*yuima@sampling@Terminal[1]
    #     # the limit of the Mardia's kurtosis based opration repetition 
    #     klarinc<-numeric(floor(limMard)) 
    #     statis<-c()
    #     trial<-1
    #     jumpdec<-NULL
    #     rqmle<<-NULL
    #     result<-NULL
    #     while((Mard>qchisq(alpha,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE))){
    #       klarinc[trial]<-which(ainc==oinc[trial]) ## detect the largest increment
    #       jumpdec[[trial]]<-preservedinc[klarinc[trial],]
    #       rqmle[[trial]]<-qmle(yuima, start = start, lower = lower, upper = upper, threshold = oinc[trial]-(oinc[trial]-oinc[trial+1])/2)@coef
    #       ## calculate the modified qmle
    #       mp<-match(names(rqmle[[trial]]),parameter)
    #       resort <- rqmle[[trial]][order(mp)]
    #       
    #       for(i in 1:length(parameter))
    #       {
    #         assign(parameter[i],resort[[i]],envir=tmp.env)
    #       }
    #       
    #       
    #       d_b <- NULL
    #       for(i in 1:d.size){
    #         if(length(eval(DRIFT[[i]],envir=tmp.env))==(length(data[,1])-1)){
    #           d_b[[i]] <- DRIFT[[i]] #this part of model includes "x"(state.variable)
    #         }
    #         else{
    #           if(is.na(c(DRIFT[[i]][2]))){ #ex. yuima@model@drift=expression(0) (we hope "expression((0))")
    #             DRIFT[[i]] <- parse(text=paste(sprintf("(%s)", DRIFT[[i]])))[[1]]
    #           }
    #           d_b[[i]] <- parse(text=paste("(",DRIFT[[i]][2],")*rep(1,length(data[,1])-1)",sep=""))
    #           #vectorization
    #         }
    #       }
    #       
    #       v_a<-matrix(list(NULL),d.size,noise_number)
    #       for(i in 1:d.size){
    #         for(j in 1:noise_number){
    #           if(length(eval(DIFFUSION[[i]][[j]],envir=tmp.env))==(length(data[,1])-1)){
    #             v_a[[i,j]] <- DIFFUSION[[i]][[j]] #this part of model includes "x"(state.variable)
    #           }
    #           else{
    #             if(is.na(c(DIFFUSION[[i]][[j]][2]))){
    #               DIFFUSION[[i]][[j]] <- parse(text=paste(sprintf("(%s)", DIFFUSION[[i]][[j]])))[[1]]
    #             }
    #             v_a[[i,j]] <- parse(text=paste("(",DIFFUSION[[i]][[j]][2],")*rep(1,length(data[,1])-1)",sep=""))
    #             #vectorization
    #           }
    #         }
    #       }
    #       
    #       resi<-NULL
    #       nvari<-matrix(0,d.size,noise_number)
    #       drivec<-numeric(d.size)
    #       
    #       for(k in 1:(length(data[,1])-1)){
    #         for(i in 1:d.size)
    #         {
    #           for(j in 1:noise_number){
    #             nvari[i,j]<-eval(v_a[[i,j]],envir=tmp.env)[k]
    #           }
    #           drivec[i]<-eval(d_b[[i]],envir=tmp.env)[k]
    #         }
    #         if(sum(k == klarinc)==0) resi[[k]] <- 1/sqrt(yuima@sampling@delta[[1]])*solve(t(nvari)%*%nvari)%*%t(nvari)%*%(dx_set[k,]--withdrift*drivec)
    #         else resi[[k]] <- numeric(d.size)
    #       }
    #       
    #       mresi<-numeric(d.size)
    #       for(k in 1:(length(data[,1])-1)){
    #         mresi<-mresi+resi[[k]]
    #       }
    #       mresi<-mresi/(length(data[,1])-1)
    #       
    #       svresi<-matrix(0,d.size,d.size)
    #       for(k in 1:(length(data[,1])-1)){
    #         svresi<-svresi+(resi[[k]]-mresi)%*%t(resi[[k]]-mresi)
    #       }
    #       svresi<-svresi/(length(data[,1])-1)
    #       
    #       snr<-NULL
    #       invsvresi<-solve(svresi)
    #       U <- svd(invsvresi)$u
    #       V <- svd(invsvresi)$v
    #       D <- diag(sqrt(svd(invsvresi)$d))
    #       sqinvsvresi <- U %*% D %*% t(V) 
    #       for(k in 1:(length(data[,1])-1)){
    #         if(sum(k == klarinc)==0) snr[[k]]<-sqinvsvresi%*%(resi[[k]]-mresi)
    #         else snr[[k]] <- numeric(d.size)
    #       }
    #       fourmom<-numeric(1)
    #       for(k in 1:(length(data[,1])-1)){
    #         fourmom<-fourmom+(t(snr[[k]])%*%snr[[k]])^2
    #       }
    #       
    #       Mard<-1/((length(data[,1])-1-trial)*8*d.size*(d.size+2))*(fourmom-(length(data[,1])-1)*d.size*(d.size+2))^2
    #       statis[trial]<-Mard
    #       trial<-trial+1
    #       if(trial == floor(limMard))
    #       { 
    #         warning("Removed jumps seems to be too many. Change intensity for more removal.")
    #         result[[1]]<-jumpdec
    #         result[[2]]<-rqmle[[trial-1]]
    #         break
    #       }
    #       result[[1]]<-jumpdec
    #       result[[2]]<-rqmle[[trial-1]]     
    #     }
    #   }
    # }
    plot(odX, type="p", main="Ordered absolute-value increments with threshold", ylab="increment sizes")
    abline(odX[k-1],0,lty=5,col="red")
    # arrows(2*s.size/3,odX[k-1],2*s.size/3,odX[k-1]+2.5, col="red")
    text(4*s.size/5,odX[k-1],labels="Threshold")
  }
  result
}