# library(KernSmooth)


setMethod("show", "yuima.snr", 
          function (object){
            cat("\nCall:\n")
            print(object@call)
            cat("\nCoefficients:\n")
            print(object@coef)
          }
)


snr<-function(yuima,start,lower,upper,withdrift=FALSE){

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

    inc<-double(s.size-1)
    inc<-X[2:(s.size)]-pX
    preservedinc<-inc
    ainc<-abs(inc)
    oinc<-sort(ainc,decreasing=TRUE)
    yuima@model@jump.coeff <- list()
    yuima@model@measure <- list()
    yuima@model@jump.variable <- character(0)
    yuima@model@measure.type <- character(0)
    yuima@model@parameter@jump <- character(0)
    yuima@model@parameter@measure <- character(0)

    #yuima@model@drift<-expression((0))
    qmle<-qmle(yuima, start = start, lower = lower, upper = upper) # initial estimation
    parameter<-yuima@model@parameter@all
    mp<-match(names(qmle@coef),parameter)
    esort <- qmle@coef[order(mp)]
    
    for(i in 1:length(qmle@coef))
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
    
    
    # assign(modelstate,pX,envir=tmp.env)
    # derv<-eval(derDIFFUSION,envir=tmp.env)
    # if(length(derv)==1){
    #   derv<-rep(derv,s.size) # vectorization
    # }
    
    mresi<-mean(resi)
    vresi<-mean((resi-mresi)^2)
    snr<-(resi-mresi)/sqrt(vresi)
   
    plot(snr,main="Raw self-normalized residual")
    # h <- dpih(snr)
    # bins <- seq(min(snr)-0.1, max(snr)+0.1+h, by=h)
    # hist(snr, prob=1, breaks=bins,xlim=c(-3,3),ylim=c(0,1))
    hist(snr, freq = FALSE, breaks = "freedman-diaconis", xlim=c(min(snr),max(snr)),ylim=c(0,1))
    xval <- seq(min(snr),max(snr),0.01) # to add plot the corresponding density
    lines(xval,dnorm(xval,0,1),xlim=c(min(snr),max(snr)),ylim=c(0,1),col="red",main="Jump removed self-normalized residual")
  
    
    call <- match.call()
    final_res <- new("yuima.snr", call = call, coef = qmle@coef[order(mp)], snr = snr, model = yuima@model)
    
    
    return(final_res)
  }else{
    yuima.stop("This function currently works only for one-dimensional case.")
  }
}



