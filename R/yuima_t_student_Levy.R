# Algorithms for the density of the increment of a L\'evy t - student process

## Internal functions:

### internal characteristic function 
phi_int <- function(t,lambda,v){
  if(t==0){
    return(1)
  }else{
    t=1i*t
    my<-(sqrt(v)*abs(t))^lambda*besselK(sqrt(v)*abs(t),lambda)/(gamma(lambda)*2^(lambda-1))
    return(my)
  }
  
}
phi <- Vectorize(phi_int, vectorize.args = "t")

### FFT internal function

Inter_dens <-function (n, nu, h, alim, blim) {
  v=nu
  lambda = nu/2
  i <- 0:(n - 1)
  dx <- (blim - alim)/n
  x <- alim + i * dx
  dt <- 2 * pi/(n * dx)
  c <- -n/2 * dt
  d <- n/2 * dt
  t <- c + i * dt
  phi_t <-  phi(t,lambda,v)^h
  X <- exp(-(0 + (0+1i)) * i * dt * alim) * phi_t
  Y <- fft(X)
  density <- dt/(2 * pi) * exp(-(0 + (0+1i)) * c * x) * Y
  data.frame(i = i, t = t, characteristic_function = phi_t, 
             x = x, density = Re(density))
}

Internal_FTT_dens<-function (x, nu, h, N = 2^10) {
  # if (length(x) == 1) {
  #   alim <- min((-abs(x) - 0.5), setInf)
  #   blim <- max((abs(x) + 0.5), setSup)
  # }
  # else {
  #   xdummy <- na.omit(x[is.finite(x)])
  #   alim <- min(min(xdummy) - 0.1, setInf)
  #   blim <- max(max(xdummy) + 0.1, setSup)
  # }
  alim <- min(x)
  blim <- max(x)
  invFFT <-Inter_dens(n=N, nu=nu, h=h, alim = alim, blim = blim)
  
  dens <- spline(invFFT$x, invFFT$density, xout=x)
  return(dens$y)
  
}

### COS internal function

inter_fun_cos1 <- function(x,FK,k,low,Dt,w){
  C_fun <- cos(k*pi*(x-low)/Dt)
  f   <- sum(w*FK*C_fun)
  return(f)
}
inter_fun_cos_vect <- Vectorize(FUN=inter_fun_cos1,vectorize.args ="x")

cos_method_dens_t <-function(pts, nu, h, N=100){
  Dt  <- diff(range(pts))
  k   <- c(0:N)
  w   <- c(0.5,rep(1,N))
  char_h <- phi(t=k*pi/Dt,lambda=nu/2,v=nu)^h
  low <- min(pts) 
  FK  <- 2/Dt*Re(char_h*exp(-1i*k*pi*low/Dt))
  dens <- inter_fun_cos_vect(x=pts,FK=FK,k=k, low=low,Dt=Dt,w=w)
  return(dens)
} 

### Laguerre 

Internal_LAG_dens <-function(pts,nu,h=1, n=40){
  # without tail approx
  #if(is.null(regular_par)){
  x<- pts/sqrt(nu)
  alpha <- nu*h/2
  out <- gauss.quad(n,"laguerre",alpha=alpha)
  nodes<-out$nodes
  myw_1<-out$weights*exp(out$nodes)*(besselK(out$nodes,nu/2))^h
  arg <- cos(as.matrix(nodes)%*%t(x))
  dens <-((2^(1-nu/2)/gamma(nu/2))^h)/(pi)*t(myw_1)%*%arg
  dens <- as.numeric(dens)/sqrt(nu)
  return(dens)
}

intern_dens <- function(mydens,x,nu,sigma=1,h=1, n=40){
  dens<-mydens
  pos_low<-tail(which(dens<0 & x<=0),1L)+1
  pos_up<-which(dens<0 & x>0)[1]-1
  if(is.na(pos_low)||is.infinite(pos_low)||length(pos_low)==0){
    pos_low<-0
  }
  if(is.na(pos_up)){
    pos_up <- +Inf
  }
  B_nu <- h*gamma(0.5*nu+0.5)/(sqrt(pi)*gamma(0.5*nu))
  expnu <- -1/2*(nu+1)
  if(pos_low>1&pos_up<length(x)){
    dens_low <- B_nu*(1+x[1:pos_low]^2)^expnu #beta*sin(beta*pi/(2))*gamma(beta)/(abs(stddata[1:pos_low])^(1+beta)*pi)
    dens_up <- B_nu*(1+x[pos_up:length(x)]^2)^expnu # beta*sin(beta*pi/(2))*gamma(beta)/(abs(stddata[pos_up:length(stddata)])^(1+beta)*pi)
    dens <- c(dens_low,dens[(pos_low+1):(pos_up-1)],dens_up)
  }else{
    if(pos_low>1&!pos_up<length(x)){
      dens_low <- B_nu*(1+x[1:pos_low]^2)^expnu #beta*sin(beta*pi/(2))*gamma(beta)/(abs(stddata[1:pos_low])^(1+beta)*pi)
      dens <- c(dens_low,dens[(pos_low+1):length(x)])
    }
    if(!pos_low>1&pos_up<length(x)){
      dens_up <- B_nu*(1+x[pos_up:length(x)]^2)^expnu
      dens <- c(dens[1:(pos_up-1)],dens_up)
    }
    if(!pos_low>1&!pos_up<length(x)){
      dens<-dens
    }
  }
  return(dens)
}

mydtlevy_1vers <-function(x,nu,sigma=1,h=1, n=40){
  # without tail approx
  #if(is.null(regular_par)){
  alpha <- nu*h/2
  out <- gauss.quad(n,"laguerre",alpha=alpha)
  nodes<-out$nodes
  myw_1<-out$weights*exp(out$nodes)*(besselK(out$nodes,nu/2))^h
  int <- as.matrix(nodes)%*%t(x/sigma)
  #arg <- matrix(0,dim(int)[1],dim(int)[2])
  a <- cos(int)
  dens <-((2^(1-nu/2)/gamma(nu/2))^h)/(pi*sigma)*t(myw_1)%*%a
  dens <- as.numeric(dens)
  return(dens)
}


mydtlevy_4vers<- function(x,nu,sigma,h, n, regular_par,ImprovGibbs){
  if(is.null(regular_par)){
    L<-0
  }else{
    L<-regular_par
  }
  if(length(x)>1){
    oldx <- x
    x<-sort(unique(x))
    dens <- as.numeric(mydtlevy_1vers(x,nu,sigma,h, n))
    
    pos_low<-max(tail(which(diff(dens)<0 & x[-1]<=0),1L), tail(which(dens<0 & x<=0),1L))+1
    pos_up<-min(which(diff(dens)>0 & x[-1]>0)[1],which(dens<0 & x>0)[1])-1
    if(is.na(pos_low)){
      pos_low<-0
    }
    if(is.na(pos_up)){
      pos_up <- +Inf
    }
    
    if(pos_low>1&pos_up<length(x)){
      x_dow <- x[1:pos_low]
      dens_low <- (as.numeric(mydtlevy_1vers(x_dow-L,nu,sigma,h, n))+as.numeric(mydtlevy_1vers(x_dow+L,nu,sigma,h, n)))/2 #beta*sin(beta*pi/(2))*gamma(beta)/(abs(stddata[1:pos_low])^(1+beta)*pi)
      x_up <- x[pos_up:length(x)]
      dens_up <- (as.numeric(mydtlevy_1vers(x_up-L,nu,sigma,h, n))+as.numeric(mydtlevy_1vers(x_up+L,nu,sigma,h, n)))/2
      dens <- c(dens_low,dens[(pos_low+1):(pos_up-1)],dens_up)
      if(ImprovGibbs){
        dens<-intern_dens(dens,x,nu,sigma,h, n)
      }
    }else{
      if(pos_low>1&!pos_up<length(x)){
        x_dow <- x[1:pos_low]
        dens_low <- (as.numeric(mydtlevy_1vers(x_dow-L,nu,sigma,h, n))+as.numeric(mydtlevy_1vers(x_dow+L,nu,sigma,h, n)))/2 #beta*sin(beta*pi/(2))*gamma(beta)/(abs(stddata[1:pos_low])^(1+beta)*pi) 
        dens <- c(dens_low,dens[(pos_low+1):length(x)])
        if(ImprovGibbs){
          dens<-intern_dens(dens,x,nu,sigma,h, n)
        }
      }
      if(!pos_low>1&pos_up<length(x)){
        x_up <- x[pos_up:length(x)]
        dens_up <- (as.numeric(mydtlevy_1vers(x_up-L,nu,sigma,h, n))+as.numeric(mydtlevy_1vers(x_up+L,nu,sigma,h, n)))/2
        dens <- c(dens[1:(pos_up-1)],dens_up)
        if(ImprovGibbs){
          dens<-intern_dens(dens,x,nu,sigma,h, n)
        }
      }
      if(!pos_low>1&!pos_up<length(x)){
        dens<-dens
      }
    }
    
    pos <- match(oldx,x)
    dens <- dens[pos]
    
  }else{
    dens <- as.numeric(mydtlevy_1vers(x,nu,sigma,h, n))
    if(dens<0){
      B_nu <- h*gamma(0.5*nu+0.5)/(sqrt(pi)*gamma(0.5*nu))
      expnu <- -1/2*(nu+1)
      dens <- B_nu*(1+x^2)^expnu
    }
  }
}


rt_Levy<-function(n,nu=1,h=1,n_laguerre=180,up=7,low=-7,
                  N_grid=500001,
                  regular_par=NULL,ImprovGibbs=TRUE){
  xin<- seq(low,up,length.out=N_grid)
  dens_g <- mydtlevy_4vers(xin/sqrt(nu),nu=nu,sigma=1,h=h, n=n_laguerre, 
                            regular_par=regular_par,
                           ImprovGibbs=ImprovGibbs)/sqrt(nu)
  cond <- xin>=low & xin<=up
  cdf_g <- cumsum(dens_g*diff(xin)[1])
  cdf_g0<-cdf_g[cond]
  cond_cdf<-(cdf_g0-cdf_g0[1])/(tail(cdf_g0,1L)-cdf_g0[1])
  U=runif(n)
  # res<-as.numeric(approx(x=cond_cdf,
  #                        y=xin[cond],xout=U,method="constant",f=1)$y)
  res<-as.numeric(spline(x=cond_cdf,
                         y=xin[cond],xout=U)$y)
  return(res)
}


## NO_Approx
tdens_noapprox <-  function(x, nu, h, method="COS", N = 2^7){
  #pts <- sort(unique(c(x,seq(low,up,length.out=N_grid))))
  if(method=="COS"){
     dens0<- cos_method_dens_t(pts=x, nu=nu, h=h, N=N)
  }else{
    if(method=="FFT"){
      dens0 <- Internal_FTT_dens(x=x, nu=nu, h=h, N = N) 
      }else{
        if(method=="LAG"){
          dens0 <- Internal_LAG_dens(pts=x,nu=nu,h=h, n=N)
        }else{
          yuima.stop("method must be either 'COS' or 'FFT' or 'LAG')")
        }
      }
  }
  #pos <- match(x,pts)
  #mydens <- dens0[pos]
  return(dens0)
}


## Tail approx

approx_dens1 <- function(mydens,x,nu,h){
  dens<-mydens
  pos_low<-tail(which(dens<0 & x<=0),1L)+1
  pos_up<-which(dens<0 & x>0)[1]-1
  # pos_low<-max(tail(which(diff(dens)<0 & x[-1]<=0),1L), tail(which(dens<0 & x<=0),1L))+1
  # pos_up<-min(which(diff(dens)>0 & x[-1]>0)[1],which(dens<0 & x>0)[1])-1
  # 
  if(is.na(pos_low)||is.infinite(pos_low)||length(pos_low)==0){
    pos_low<-0
  }
  if(is.na(pos_up)){
    pos_up <- +Inf
  }
  B_nu <- h*gamma(0.5*nu+0.5)/(sqrt(pi)*gamma(0.5*nu))
  expnu <- -1/2*(nu+1)
  if(pos_low>1&pos_up<length(x)){
    dens_low <- B_nu*(1+x[1:pos_low]^2)^expnu #beta*sin(beta*pi/(2))*gamma(beta)/(abs(stddata[1:pos_low])^(1+beta)*pi)
    dens_up <- B_nu*(1+x[pos_up:length(x)]^2)^expnu # beta*sin(beta*pi/(2))*gamma(beta)/(abs(stddata[pos_up:length(stddata)])^(1+beta)*pi)
    dens <- c(dens_low,dens[(pos_low+1):(pos_up-1)],dens_up)
  }else{
    if(pos_low>1&!pos_up<length(x)){
      dens_low <- B_nu*(1+x[1:pos_low]^2)^expnu #beta*sin(beta*pi/(2))*gamma(beta)/(abs(stddata[1:pos_low])^(1+beta)*pi)
      dens <- c(dens_low,dens[(pos_low+1):length(x)])
    }
    if(!pos_low>1&pos_up<length(x)){
      dens_up <- B_nu*(1+x[pos_up:length(x)]^2)^expnu
      dens <- c(dens[1:(pos_up-1)],dens_up)
    }
    if(!pos_low>1&!pos_up<length(x)){
      dens<-dens
    }
  }
  return(dens)
}
# Inter_mytail2<-function(dens,x,pos, pos2,C_nu){
#   t_low <- 1/(1+x[pos]^2)
#   t_low_inc <- 1/(1+x[pos2]^2)
#   f<- dens[pos]
#   delta_t_low <- t_low-t_low_inc
#   delta_f_low <- (f-dens[pos2])/(t_low-t_low_inc)
#   C_nu_low <- C_nu*t_low^(0.5*nu+.5)
#   Mat_low <- matrix(c(2*t_low,-1,-t_low^2,t_low),2,2)/(C_nu_low*t_low^2)
#   Coeff_low <- Mat_low%*%(c(f,delta_f_low)-c(C_nu_low,1/2*(nu+1)*t_low^(-1)*f))
#   return(Coeff_low)
# }
# 
# 
# Inter_mytail<-function(dens,x,pos,C_nu){
#   t_low <- 1/(1+x[pos]^2)
#   f<- dens[pos]
#   C_nu_low <- C_nu*t_low^(0.5*nu+.5)
#   Coeff_low <- (f/C_nu_low-1)/t_low
#   return(Coeff_low)
# }
# 
# 
# approx_dens2 <- function(mydens,x,nu,h){
#   dens<-mydens
#   pos_low<-tail(which(dens<0 & x<=0),1L)+1
#   pos_up<-which(dens<0 & x>0)[1]-1
#   # pos_low<-max(tail(which(diff(dens)<0 & x[-1]<=0),1L), tail(which(dens<0 & x<=0),1L))+1
#   # pos_up<-min(which(diff(dens)>0 & x[-1]>0)[1],which(dens<0 & x>0)[1])-1
#   
#   if(is.na(pos_low)||is.infinite(pos_low)||length(pos_low)==0){
#     pos_low<-0
#   }
#   if(is.na(pos_up)){
#     pos_up <- +Inf
#   }
#   C_nu<- h*gamma(0.5*nu+0.5)/(sqrt(pi)*gamma(0.5*nu))
#   if(pos_low>1&pos_up<length(x)){
#     #pos_low2<-pos_low+1
#     #Coeff_low <- Inter_mytail(dens=dens,x=x,pos=pos_low, pos2=pos_low2,C_nu=C_nu)
#     Coeff_low <- Inter_mytail(dens=dens,x=x,pos=pos_low,C_nu=C_nu)
#     t_low <- 1/(1+x[1:pos_low]^2)
#     #dens_low<- C_nu*t_low^(0.5*nu+0.5)*(1+Coeff_low[1,1]*t_low+Coeff_low[2,1]*t_low^2)
#     dens_low<- C_nu*t_low^(0.5*nu+0.5)*(1+Coeff_low*t_low)
#                                         
#     # pos_up2<-pos_up-1
#     # Coeff_up <- Inter_mytail(dens=dens,x=x,pos=pos_up, pos2=pos_up2,C_nu=C_nu)
#     Coeff_up <- Inter_mytail(dens=dens,x=x,pos=pos_up,C_nu=C_nu)
#     t_up <- 1/(1+x[pos_up:length(x)]^2)
#     #dens_up<- C_nu*t_up^(0.5*nu+0.5)*(1+Coeff_up[1,1]*t_up+Coeff_up[2,1]*t_up^2)
#     dens_up<- C_nu*t_up^(0.5*nu+0.5)*(1+Coeff_up*t_up)
#     dens <- c(dens_low,dens[(pos_low+1):(pos_up-1)],dens_up)
#   }else{
#     if(pos_low>1&!pos_up<length(x)){
#       #dens_low <- B_nu*(1+x[1:pos_low]^2)^expnu #beta*sin(beta*pi/(2))*gamma(beta)/(abs(stddata[1:pos_low])^(1+beta)*pi)
#       # pos_low2<-pos_low+1
#       # Coeff_low <- Inter_mytail(dens=dens,x=x,pos=pos_low, pos2=pos_low2,C_nu=C_nu)
#       Coeff_low <- Inter_mytail(dens=dens,x=x,pos=pos_low, C_nu=C_nu)
#       t_low <- 1/(1+x[1:pos_low]^2)
#       #dens_low<- C_nu*t_low^(0.5*nu+0.5)*(1+Coeff_low[1,1]*t_low+Coeff_low[2,1]*t_low^2)
#       dens_low<- C_nu*t_low^(0.5*nu+0.5)*(1+Coeff_low*t_low)
#       dens <- c(dens_low,dens[(pos_low+1):length(x)])
#     }
#     if(!pos_low>1&pos_up<length(x)){
#       #dens_up <- B_nu*(1+x[pos_up:length(x)]^2)^expnu
#       #pos_up2<-pos_up-1
#       #Coeff_up <- Inter_mytail(dens=dens,x=x,pos=pos_up, pos2=pos_up2,C_nu=C_nu)
#       Coeff_up <- Inter_mytail(dens=dens,x=x,pos=pos_up, C_nu=C_nu)
#       t_up <- 1/(1+x[pos_up:length(x)]^2)
#       #dens_up<- C_nu*t_up^(0.5*nu+0.5)*(1+Coeff_up[1,1]*t_up+Coeff_up[2,1]*t_up^2)
#       dens_up<- C_nu*t_up^(0.5*nu+0.5)*(1+Coeff_up*t_up)
#       dens <- c(dens[1:(pos_up-1)],dens_up)
#     }
#     if(!pos_low>1&!pos_up<length(x)){
#       dens<-dens
#     }
#   }
#   
# }


# yuima density function

dtLevy<-function(x, nu, h, method="COS",  up = 7, low = -7, N = 2^7, N_grid=1000, regular_par=NULL){
  pts <- sort(unique(c(x,seq(low,up,length.out=N_grid))))
  
  if(method == "COS" | method == "FFT"){
    approx_tail <- "NO"
  }else{
    approx_tail <- "Tail_Approx1"
  }
  
  if(method=="LAG" & approx_tail=="Tail_Approx1"){
    mydens <- mydtlevy_4vers(pts/sqrt(nu),nu=nu,sigma=1,h=h, n=N, 
                             regular_par=regular_par,
                             ImprovGibbs=TRUE)/sqrt(nu)
    pos <- match(x,pts)
    mydens1 <-mydens[pos]
    return(mydens1)
  }
  
  dens0<-tdens_noapprox(x=pts, nu, h, method, N)
  
  #if(approx_tail=="NO"){
    pos <- match(x,pts)
    mydens <- dens0[pos]
    return(mydens)
  #}else{
    # if(approx_tail=="Tail_Approx1"){
    #   mydens <- approx_dens1(mydens=dens0,x=pts,nu=nu,h=h)
    #   pos <- match(x,pts)
    #   mydens1 <-mydens[pos]
    #   return(mydens1)
    # }
    # if(approx_tail=="Tail_Approx2"){
    #   mydens <- approx_dens2(mydens=dens0,x=pts,nu=nu,h=h)
    #   pos <- match(x,pts)
    #   mydens1 <-mydens[pos]
    #   return(mydens1)
    # }
  #}
}

# yuima cdf

# cdf
ptLevy<-function(x, nu, h, method="COS",  up = 7, low = -7, N = 2^7, N_grid=1000,regular_par=NULL){
  pts <- sort(unique(c(x,seq(low,up,length.out=N_grid))))
  if(method == "COS" | method == "FFT"){
    approx_tail <- "NO"
  }else{
    approx_tail <- "Tail_Approx1"
  }
  
  if(method=="LAG" & approx_tail=="Tail_Approx1"){
    dens_g <- mydtlevy_4vers(pts/sqrt(nu),nu=nu,sigma=1,h=h, n=N, 
                             regular_par=regular_par,
                             ImprovGibbs=TRUE)/sqrt(nu)
  }else{
    dens_g <- dtLevy(x=pts, nu=nu, h=h, method=method, up=max(pts), low=min(pts), N=N, N_grid=2)
  }
  cdf_g <- cumsum(dens_g*diff(pts)[1])
  cond <- pts>=low & pts<=up
  cdf_g0<-cdf_g[cond]
  cond_cdf<-(cdf_g0-cdf_g0[1])/(tail(cdf_g0,1L)-cdf_g0[1])
  
  res<-as.numeric(spline(x=pts[cond],
                         y=cond_cdf,xout=x)$y)
  return(res)
}

# yuima quantile

# quantile
qtLevy<-function(p, nu, h, method="COS",  up = 7, low = -7, N = 2^7, N_grid=1000,regular_par=NULL){
  pts <- seq(low,up,length.out=N_grid)
  if(method == "COS" | method == "FFT"){
    approx_tail <- "NO"
  }else{
    approx_tail <- "Tail_Approx1"
  }
  if(method=="LAG" & approx_tail=="Tail_Approx1"){
    dens_g <- mydtlevy_4vers(pts/sqrt(nu),nu=nu,sigma=1,h=h, n=N, 
                             regular_par=regular_par,
                             ImprovGibbs=TRUE)/sqrt(nu)
  }else{
    dens_g<-dtLevy(x=pts, nu=nu, h=h, method=method,  up=max(pts), low=min(pts), N=N, N_grid=2)
  }
  cdf_g <- cumsum(dens_g*diff(pts)[1])
  cond <- pts>=low & pts<=up
  cdf_g0<-cdf_g[cond]
  cond_cdf<-(cdf_g0-cdf_g0[1])/(tail(cdf_g0,1L)-cdf_g0[1])
  res<-as.numeric(approx(x=cond_cdf,
                         y=pts[cond],xout=p,method="constant",f=1)$y)
  
  # res<-as.numeric(spline(x=cond_cdf,
  #                        y=pts[cond],xout=p)$y)
  return(res)
}

# yuima rng

# Random Number Generator
rtLevy<-function(n, nu, h, method="COS", up = 7, low = -7, N = 2^7, N_grid=1000,regular_par=NULL){
  if(method == "COS" | method == "FFT"){
    approx_tail <- "NO"
  }else{
    approx_tail <- "Tail_Approx1"
  }
  if(method=="LAG" & approx_tail=="Tail_Approx1"){
    res<-rt_Levy(n=n,nu=nu,h=h,n_laguerre=N,up=up,low=low,
                 N_grid=c(N_grid+1),
                 regular_par=regular_par, ImprovGibbs=TRUE)
    return(res)
  }
  
  pts <- seq(low,up,length.out=N_grid)
  dens_g<-dtLevy(x=pts, nu=nu, h=h, method=method, up=max(pts), low=min(pts), N=N, N_grid=2)
  cdf_g <- cumsum(dens_g*diff(pts)[1])
  cond <- pts>=low & pts<=up
  cdf_g0<-cdf_g[cond]
  cond_cdf<-(cdf_g0-cdf_g0[1])/(tail(cdf_g0,1L)-cdf_g0[1])
  U <- runif(n)
  # res<-as.numeric(approx(x=cond_cdf,
  #                        y=pts[cond],xout=U,method="constant",f=1, ties = mean)$y) # if we use constant we have more zeros
  res<-as.numeric(approx(x=cond_cdf,
                         y=pts[cond],xout=U,method="linear",f=1, ties = mean)$y) # if we use constant we have more zeros
  # res<-as.numeric(spline(x=cond_cdf,
  #                   y=pts[cond],xout=U)$y)
  return(res)
  # newpts <- seq(pts[which(cdf_g0>0.001)[1]], pts[which(cdf_g0>1-0.001)[1]-1],length.out=N_grid)
  # 
  # dens_g<-dtLevy(x=newpts, nu=nu, h=h, method=method, 
  #                up=max(newpts), low=min(newpts), N=N, N_grid=2)
  # cdf_g0 <- cumsum(dens_g*diff(pts)[1])
  #   
  # cond_cdf<-(cdf_g0-cdf_g0[1])/(tail(cdf_g0,1L)-cdf_g0[1])
  # U <- runif(n)
  # # res<-as.numeric(approx(x=cond_cdf,
  # #                        y=newpts,xout=U,method="constant",f=1)$y) # if we use constant we have more zeros
  # res<-as.numeric(spline(x=cond_cdf,
  #                        y=newpts,xout=U)$y)
  # # U <- runif(n)
  # # res<-as.numeric(approx(x=cdf_g,
  # #                        y=pts,xout=U,method="constant",f=1)$y)
  # 
  # 
  # return(res)

}

