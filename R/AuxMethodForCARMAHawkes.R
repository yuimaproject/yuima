is.CarmaHawkes<-function(obj){
  if(is(obj,"yuima"))
    return(is(obj@model, "yuima.carmaHawkes"))
  if(is(obj,"yuima.model"))
    return(is(obj, "yuima.carmaHawkes"))
  return(FALSE)
}
# CARMA Hawkes Utilities
Atilde<-function(a,b){
  # a <- [a_1,...a_p]
  # b <- [b_0,...,b_{p-1}]
  d <- length(a)
  btrue<-numeric(length =d)
  btrue[1:length(b)]<-b 
  lastrow <- btrue-rev(a)
  if(d>1){
    A <- cbind(0,diag(1,d-1,d-1))
    Atilde <-rbind(A,lastrow)
  }else{
    Atilde <- as.matrix(lastrow)
  }
  
  rownames(Atilde) <- rep(" ",d)
  return(Atilde)
}
myebold<-function(d){
  res <- matrix(0, d,1) 
  res[d]<-1
  return(res)
}
bbold<-function(b,d){
  btrue<-numeric(length =d)
  btrue[1:length(b)]<-b
  return(btrue)
}
Atildetilde<-function(a,b,p){
  res <- matrix(0, p*(p+1)/2, p*(p+1)/2)
  aaa<- length(b)
  btrue<-numeric(length =p)
  btrue[1:aaa]<-b 
  b<-btrue
  for(j in c(1:p)){
    for(i in c(1:p)){
      dimrow <- p-j+1
      dimcol <- p-i+1
      if(dimrow==dimcol){
        #DMat<- Matr(0,dimrow,dimcol)
        posdf<- j*p-sum(c(0:(j-1)))
        if(dimrow!=1){
          id<-matrix(0,dimrow,dimcol)
          id[1,2]<-1
          DMat<-Atilde(a=a[1:dimrow],b[j:aaa]) + id 
        }else{
          
          DMat<- as.matrix(2*tail(btrue-rev(a),1L))
        }
        res[(posdf-dimrow+1):posdf,(posdf-dimrow+1):posdf]<-DMat
      }
      if(dimrow<dimcol){
        posdfcol<- i*p-sum(c(0:(i-1)))
        posdfrow<- j*p-sum(c(0:(j-1)))
        LMat<- matrix(0,dimrow,dimcol)
        
        if(dimrow!=1){
          LMat[dim(LMat)[1],j-i+1]<-b[i]-a[p-i+1]
        }else{
          LMat[dim(LMat)[1],j-i+1]<-2*(b[i]-a[p-i+1])
        }        
        res[(posdfrow-dimrow+1):posdfrow,(posdfcol-dimcol+1):posdfcol]<-LMat
      }
      if(dimrow==dimcol+1){
        posdfcol<- i*p-sum(c(0:(i-1)))
        posdfrow<- j*p-sum(c(0:(j-1)))
        UMat<- rbind(0,diag(dimcol))
        res[(posdfrow-dimrow+1):posdfrow,(posdfcol-dimcol+1):posdfcol]<-UMat
      }
    }
  }
  return(res)
}
BboldFunct<-function(b,p){
  aaa<- length(b)
  btrue<-numeric(length =p)
  btrue[1:aaa]<-b 
  b<-btrue
  
  #Bbold[1,]<-b
  if(p>1){
    Bbold<- diag(rep(b[1],p))
    Bbold[1,]<-b
    for(i in c(2:p)){
      A0<-matrix(0,i-1,p-i+1)
      if(p-i+1!=1){
        Ai<-diag(rep(b[i],p-i+1))
      }else{Ai <- as.matrix(b[i])}
      Ai[1,]<-b[i:p]
      dumA <- rbind(A0,Ai)
      Bbold<- cbind(Bbold,dumA)
    }
  }else{
    Bbold<-as.matrix(b)
  } 
  return(Bbold)
}
Ctilde <- function(mu,b,p){
  if(p>1){  
    Ctilde <- matrix(0,p,p)
    Ctilde[p,1] <- mu
    for(j in c(1:(p-1))){
      dum<- matrix(0,p-j,p)
      dum[p-j,j+1]<- mu
      Ctilde<- rbind(Ctilde,dum)
    }
    aaa<- length(b)
    btrue<-numeric(length =p)
    btrue[1:aaa]<-b 
    dumVec <- numeric(length =p)
    dumVec[p]<- 2*mu
    vect <- btrue+dumVec
    Ctilde[p*(p+1)/2,]<- vect
  }else{
    Ctilde <- as.matrix(b+2*mu)
  }
  return(Ctilde)
}

ExpIncr<-function(mu,b,a,X0,t0,ti1,ti2){
  At<-Atilde(a,b)
  Ainv<-solve(At)
  d <- length(a)
  ebold <-myebold(d)
  bb<-bbold(b,d)
  term1 <-t(bb)%*%Ainv
  res <- mu*(1-term1%*%ebold)*(ti2-ti1)
  res <- res+term1%*%(expm(At*(ti2-t0))-expm(At*(ti1-t0)))%*%(X0+Ainv%*%ebold*mu)
  return(res)
}

tripletExpInt <- function(A11,A12,A22,A23,A33,T1){
  # dA <- dim(A)
  # dB <- dim(B)
  # bigA<- cbind(A,C)
  # dummy <-matrix(0,dB[1],dA[2])
  # big1<- cbind(dummy,B)
  # bigA<-rbind(bigA,big1)
  dA11 <- dim(A11)
  dA12 <- dim(A12)
  dA22<-dim(A22)
  dA23 <- dim(A23)
  dA33 <- dim(A33)
  dummy <-matrix(0,dA11[1],dA23[2])
  bigA<- cbind(A11,A12,dummy)
  dummy <-matrix(0,dA22[1],dA23[2])
  big1 <- cbind(dummy,A22,A23)
  bigA<-rbind(bigA,big1)
  dummy<- matrix(0,dA33[1],dA11[2]+dA22[2])
  big1 <- cbind(dummy,A33)
  bigA<-rbind(bigA,big1)
  return(bigA)
  #res <- expm(bigA*T1)[1:dA[1], 1:dB[2]+dA[2]]
}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

doubleExpInt <- function(A,C,B,T1){
  dA <- dim(A)
  dB <- dim(B)
  bigA<- cbind(A,C)
  dummy <-matrix(0,dB[1],dA[2])
  big1<- cbind(dummy,B)
  bigA<-rbind(bigA,big1)
  res <- expm(bigA*T1)[1:dA[1], 1:dB[2]+dA[2]]
}


VtrillX0 <- function(X0, Att, At, Ct, mu, p, t0=0, T1){
  e<- matrix(0,p,1)
  e[p,1]<-1
  et <- matrix(0,p*(p+1)/2,1)
  et[p*(p+1)/2,1]<-1
  XX0<- as.matrix(X0)%*%t(as.matrix(X0))
  vtXX0<-as.matrix(XX0[lower.tri(XX0,T)])
  AttInv<-solve(Att)
  AtInv <-solve(At)
  term1 <- expm(Att*(T1-t0))%*%(vtXX0+mu*AttInv%*%et-mu*AttInv%*%Ct%*%AtInv%*%e )# ->0 as T->+\infty
  term2 <- mu*AttInv%*%Ct%*%AtInv%*%e - mu*AttInv%*%et
  
  B12<-doubleExpInt(A=Att,C=Ct,B=At,T1=T1)
  term3 <- B12%*%expm(-At*t0)%*%(X0+AtInv%*%e*mu)
  return(list(res =term1+term2+term3, term1=term1, term2=term2, term3=term3))
}

gfunct <- function(a,b,mu){
  p<-length(a)
  At<-Atilde(a,b)
  eb<-myebold(p)
  
  AtInv <- solve(At) 
  Att <- Atildetilde(a,b,p)
  AttInv <- solve(Att)
  etilde <- myebold(p*(p+1)/2)
  B <- BboldFunct(b,p)
  Ct <- Ctilde(mu,b,p)
  
  aaa <- length(b)
  btrue<-numeric(length =p)
  btrue[1:aaa]<-b 
  b<-btrue
  
  const <- AtInv*mu
  Term1 <- as.numeric(t(b)%*%AtInv%*%eb)
  res <- eb*Term1-eb+AtInv%*%eb*mu*Term1+B%*%AttInv%*%(etilde-Ct%*%AtInv%*%eb)
  return(mu*AtInv%*%res)
}

# Function For simulation of the trajectories
S_Kfunction <- function(S_old, T_k, T_old, A, Id){
  res<- expm(A*(T_k-T_old))%*%(S_old + Id)
  return(res)
}

S_KfunctionSim <- function(S_old, T_k, T_old, A, Id){
  res<- expm(A*(T_k-T_old))%*%(S_old ) + Id
  return(res)
}

ConditionForSimul <-function(x,U,S_k, T_k, A, Ainv, Id, bbold1, ebold1, mu, t0=0, X0){
  res0 <- mu*(x-T_k)+t(bbold1)%*%expm(A*(T_k-t0))%*%(Ainv%*%(expm(A*(x-T_k))-Id)%*%X0)
  res0 <-res0+t(bbold1)%*%S_k%*%(Ainv%*%(expm(A*(x-T_k))-Id)%*%ebold1)
  res<-log(U)+as.numeric(res0)
  return(res)
}

simulateCarmaHawkes <-function(mu,b,a,X0,t0=0,FinalTime){
  p<- length(a)
  A <- Atilde(a, numeric(length =p))
  Ainv <- solve(A) 
  bbold1 <- bbold(b,p)
  ebold1 <- myebold(p)
  JumpTime <- NULL
  U<-runif(1)
  JumpTime <- -log(U)/mu
  if(JumpTime<FinalTime){
    Cond<-TRUE  
  }else{
    return(NULL)
  }
  S_old <- diag(1,p,p)
  Id <- diag(1,p,p)
  T_old<-JumpTime
  S_k <- S_old 
  T_k <- JumpTime
  while(Cond){
    U<-runif(1)
    # S_old <-S_Kfunction(S_k, T_k, T_old, A, Id)
    S_k <-S_KfunctionSim(S_k, T_k, T_old, A, Id)
    T_old <- T_k
    res<-uniroot(f=ConditionForSimul,interval=c(T_k, 10*FinalTime),
                 U=U,S_k=S_k, T_k=T_k, 
                 A=A, Ainv=Ainv, Id=Id, 
                 bbold1=bbold1, ebold1=ebold1, mu=mu,
                 X0=X0)
    T_k<-res$root
    
    
    
    if(T_k<=FinalTime){
      JumpTime<-c(JumpTime,T_k)
      #cat("\n",T_k)
      f_u<-ConditionForSimul(x=10*FinalTime,U,S_k, T_k, A, Ainv, Id, bbold1, ebold1, mu, t0=0, X0)
      f_l<-ConditionForSimul(x=T_k,U,S_k, T_k, A, Ainv, Id, bbold1, ebold1, mu, t0=0, X0)
      if(f_u*f_l>=0){
        Cond=FALSE
       # cat("\n vedi ", T_k)
      }
    }else{
      Cond <- FALSE
    }
    
  }
  #CountProc <- 1:length(JumpTime)
  JumpTime<- c(0,JumpTime)
  return(JumpTime)
}

aux.simulateCarmaHawkes<- function(object, true.parameter){
  model <- object@model
  ar.par <- true.parameter[model@info@ar.par]
  ma.par <- true.parameter[model@info@ma.par]
  if(!model@info@XinExpr){
    X0 <- rep(0,model@info@p)
  }else{
    yuima.stop("X0 will be available as soon as possible")
  }
  mu <- true.parameter[model@info@base.Int]
  t0 <- object@sampling@Initial[1]
  FinalTime <- object@sampling@Terminal[1]
  res <- simulateCarmaHawkes(mu = mu,b = ma.par,a = ar.par,
           X0 = X0,t0 = t0, FinalTime = FinalTime)
  numb<-length(res)
  
  if(length(object@model@measure$df@param.measure)==0){
    param <- 1
    names(param)<- object@model@measure$df@time.var
    Nt<-cumsum(c(0,rand(object@model@measure$df, n=numb, param=param)))
  }else{
    yuima.stop("Jump size different from 1 will be available as soon as possible")
  }
  Nt<-matrix(Nt)
  colnames(Nt)<-object@model@info@Counting.Process
  res1<-zoo(x=Nt, order.by=res)
  mydata <- setData(original.data = res1)
  obsgrid <- na.approx(res1,xout=object@sampling@grid[[1]], method ="constant")
  mydata@zoo.data[[1]]<-obsgrid 
  object@data <- mydata
  object@model@solve.variable <- object@model@info@Counting.Process # Check Again
  return(object)
}

SimThinAlg <- function(x,FinalT, p, q){
  mu0<-x[1]
  a0<-x[1:p+1]
  b0<-x[1:(q+1)+p+1]
  b0 <- bbold(b0,p)
  A <- Atilde(a0, numeric(length =p))
  eigenvalues <- sort(eigen(A)$values,decreasing = TRUE)
  dumexp <- 1:p-1
  S <- kronecker(t(eigenvalues),as.matrix(dumexp),"^")
  #S_new%*%diag(eigenvalues)%*%solve(S_new)
  # 
  coef<-as.matrix(sqrt(sum(abs(t(b0)%*%S)^2))*sqrt(sum(abs(solve(S)%*%myebold(p))^2)))
  expcoef <- as.matrix(-max(Re(eigenvalues)))
  # Initialize
  JumpTime <- NA
  s <- 0
  n <- 0
  U<-runif(1)
  JumpTime[1] <- -log(U)/mu0
  S_k <- diag(1,p,p)
  Id <- diag(1,p,p)
  S_k_H <- diag(1)
  Id_H <- diag(1)
  s <- s+ JumpTime[1]
  T_old <- 0
  T_k <-  JumpTime[1]
  S_k <-S_KfunctionSim(S_k, T_k, T_old, A, Id)
  S_k_H <-S_KfunctionSim(S_k_H, T_k, T_old, -expcoef, Id_H)
  bbold1 <- bbold(b0,p)
  ebold1 <- myebold(p)
  n<-1
  #i=0
  while(s<FinalT){
    U<-runif(1)
    # S_old <-S_Kfunction(S_k, T_k, T_old, A, Id)
    lambda_k_bar <- as.numeric(mu0+exp(-expcoef*(s-T_k))*coef*S_k_H)
    omega <- -log(U)/lambda_k_bar
    s <- s+omega
    D <- runif(1)
    lambda_k <- mu0+ as.numeric(t(bbold1)%*%expm(A*(s-T_k))%*%S_k%*%ebold1)
    if(lambda_k_bar*D<lambda_k){
      n<-n+1
      T_old <- T_k
      T_k <- s
      JumpTime <- c(JumpTime, T_k)
      S_k <-S_KfunctionSim(S_k, T_k, T_old, A, Id)
      S_k_H <-S_KfunctionSim(S_k_H, T_k, T_old, -expcoef, Id_H)
    }
    #i<-i+1
    # cat("\n",i)
  }
  if(T_k <=FinalT){
    res <- zoo(0:n,order.by = c(0,JumpTime))
  }else{
    JumpTime <- JumpTime[-n]
    res <- zoo(1:n-1, order.by =c(0,JumpTime))
  }
  return(res)
}


aux.simulateCarmaHawkes_thin<- function(object, true.parameter){
  model <- object@model
  p <- length(model@info@ar.par)
  q <- length(model@info@ma.par)
  if(!model@info@XinExpr){
    X0 <- rep(0,p)
  }else{
    yuima.stop("X0 will be available as soon as possible")
  }
#  mu <- true.parameter[model@info@base.Int]
 # t0 <- object@sampling@Initial[1]
  true.parameter <- true.parameter[c(model@info@base.Int,model@info@ar.par,model@info@ma.par)]
  FinalTime <- object@sampling@Terminal[1]
  
  res <- SimThinAlg(x=true.parameter,FinalT=FinalTime, p = p, q = q-1)
  numb<-length(res)
  
  if(length(object@model@measure$df@param.measure)==0){
    param <- 1
    names(param)<- object@model@measure$df@time.var
    Nt<-cumsum(c(0,rand(object@model@measure$df, n=numb, param=param)))
  }else{
    yuima.stop("Jump size different from 1 will be available as soon as possible")
  }
  Nt<-matrix(Nt)
  colnames(Nt)<-object@model@info@Counting.Process
  res1<-zoo(x=Nt, order.by=index(res))
  mydata <- setData(original.data = res1)
  obsgrid <- na.approx(res1,xout=object@sampling@grid[[1]], method ="constant")
  mydata@zoo.data[[1]]<-obsgrid 
  object@data <- mydata
  object@model@solve.variable <- object@model@info@Counting.Process # Check Again
  return(object)
}
Int_LikelihoodCarmaHawkes <- function(x,p,q,JumpTime,T_k,
                                        numbOfJump,
                                        t0 = 0, Id = diag(p),
                                        X0 = matrix(0,p),
                                        display = FALSE){
  mu<-x[1]
  a<-x[1:p+1]
  b<-x[1:(q+1)+p+1]
  A <- Atilde(a,rep(0,p))
  InvA <- solve(A)
  bbold1 <- bbold(b,p)
  ebold1 <- myebold(p)
  res0<-mu*(T_k-t0)
  dum <- t(bbold1)%*%InvA
  res1 <- dum%*%(expm(A*(T_k-t0))-Id)%*%X0
  lambda<-numeric(length=numbOfJump)
  S<-Id*0 #S_Kfunction(Id, T_k=JumpTime[1], JumpTime[1], A, Id) hawkes
  #S<-S_Kfunction(Id, T_k=JumpTime[1], JumpTime[1], A, Id)
  # interstep <- expm(A*(JumpTime[1]-t0))%*%X0+S%*%ebold1
  lambda[1]<-mu# + as.numeric(t(bbold1)%*%interstep)
  for(i in c(2:numbOfJump)){
    S<-S_Kfunction(S, T_k=JumpTime[i], JumpTime[i-1], A, Id)
    interstep <- expm(A*(JumpTime[i]-t0))%*%X0+S%*%ebold1
    lambda[i]<-mu+as.numeric(t(bbold1)%*%interstep)
  }
  
  if(any(is.infinite(lambda))){
    res <--10^6
  }else{
    if(any(is.nan(lambda))){
      res <--10^6
    }else{
      if(all(lambda>0)){
        res1 <- res1+dum%*%S%*%ebold1-numbOfJump*dum%*%ebold1
        res <--res0-as.numeric(res1)+sum(log(lambda[-1]))
        if(display){cat("\n", c(res,x))}
      }else{
        res<--10^6
      }
    }  
    
  }
  return(-res/numbOfJump)
}


EstimCarmaHawkes <- function(yuima, start, est.method = "qmle", method = "BFGS",
        lower = NULL, upper = NULL, lags = NULL, display = FALSE){
  model <- yuima@model
  names.par <- c(model@info@base.Int, model@info@ar.par, model@info@ma.par) 
  cond0 <- all(names(start) %in% names.par) #& all(names.par %in% names(start))
  if(!cond0){
    yuima.stop(paste("names in start are ", 
                     paste(names(start), collapse = ", "), " while param names in the model are ", 
                     paste(names.par, collapse = ", "), collapse =""))
  }
  true.par <- start[names.par]
  p <- model@info@p
  q <- model@info@q
  t0 <- yuima@sampling@Initial[1]
  if(est.method=="qmle"){
    JumpTime <- time(yuima@data@original.data)
    T_k <- tail(JumpTime,1L)
    #t0 <- JumpTime[1]
    if(is.null(lower) & is.null(upper)){
        res <- optim(par=true.par,fn = Int_LikelihoodCarmaHawkes,
                 p = p,q=q, JumpTime = JumpTime,
                 T_k= T_k,
                 numbOfJump=tail(as.numeric(yuima@data@original.data),1L),
                 t0=t0,Id=diag(1,p,p), display = display,
                 method = method)
    }else{
      yuima.stop("constraints will be available as soon as possible.")
    }
    res$value <- res$value*T_k
  }else{
    yuima.warn("We estimate the Carma(p,q)-Hawkes using the empirical autocorrelation function")
    yuima.stop("This method is not available yet!!! It will be implemented as soon as possible")
  }
  return(res)
}