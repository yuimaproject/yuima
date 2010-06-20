# funcF
# function to calculate F in (13.2)
funcF <- function(yuima,X,e){ 
            division <- yuima@sampling@n  #number of observed time
            F <- getF(yuima@functional)
            d.size <- yuima@model@equation.number
            k.size <- length(F)
            modelstate <- yuima@model@state.variable
            XT <- X[division,]   #X observed last. size:vector[d.size]
            resF <- numeric(k.size)   #values of F. size:vector[k.size]
            
            for(d in 1:d.size){
              assign(modelstate[d],XT[d])  #input XT in state("x1","x2") to use eval function to F
            }  
            for(k in 1:k.size){
              resF[k] <- eval(F[k])  #calculate F with XT
            }
            return(resF)
          }

# funcf
# function to calculate fa in (13.2)
funcf <- function(yuima,X,e){ 
            division <- yuima@sampling@n  #number of observed time
            F <- getF(yuima@functional)
            f <- getf(yuima@functional)
            r.size <- yuima@model@noise.number
            d.size <- yuima@model@equation.number
            k.size <- length(F)
            modelstate <- yuima@model@state.variable
            resf <- array(0,dim=c(k.size,division,(r.size+1))) #value of f0,f1,~,fr. size:array[k.size,division,r.size+1]
            
            for(r in 1:(r.size+1)){
              for(t in 1:division){
                Xt <- X[t,]     #Xt is data X observed on time t. size:vector[d.size] 
                for(d in 1:d.size){
                  assign(modelstate[d],Xt[d]) #input Xt in state to use eval function to f
                }  
                for(k in 1:k.size){
                  resf[k,t,r] <- eval(f[[r]][k]) #calculate k th expression of fr.
                }
              }
            }
            return(resf)
          }

# funcFe.
# function to calculate Fe in (13.2).
# core function of 'simFunctional'
funcFe. <- function(yuima,X,e){
            F <- getF(yuima@functional)
            r.size <- yuima@model@noise.number
            d.size <- yuima@model@equation.number
            k.size <- length(F)
            modelstate <- yuima@model@state.variable
            division <- yuima@sampling@n
            Terminal <- yuima@sampling@Terminal
            delta <- Terminal/division  #length between observed times
            dw <- matrix(0,division-1,r.size+1) #Wr(t) size:matrix[division,r.size+1]
            
            dw[,1] <- rep(Terminal/(division-1),length=division-1) #W0(t)=t
            
            for(r in 2:(r.size+1)){
              tmp <- rnorm(division,0,sqrt(delta)) #calculate Wr(t)
              dw[,r] <- tmp[2:division]-tmp[1:(division-1)] #calculate dWr(t)
            }
            
            resF <- funcF(yuima,X,e) #calculate F with X,e. size:vector[k.size]
            resf <- funcf(yuima,X,e) #calculate f with X,e. size:array[k.size,division,r.size+1]
            
            Fe <- numeric(k.size)
            for(k in 1:k.size){
              Fe[k] <- sum(resf[k,1:division-1,]*dw)+resF[k]  #calculate Fe using resF and resf as (13.2).
            }
            return(Fe)            
          }

# Get.X.t0
# function to calculate X(t=t0) using runge kutta method
Get.X.t0 <- function(yuima){
            r.size <- yuima@model@noise.number
            d.size <- yuima@model@equation.number
            k.size <- length(getF(yuima@functional))
            modelstate <- yuima@model@state.variable
            V0 <- yuima@model@drift
            V <- yuima@model@diffusion

            Terminal <- yuima@sampling@Terminal
            division <- yuima@sampling@n
            
            pars <- yuima@model@parameter@all[1]  #epsilon
            xinit <- getxinit(yuima@functional)
            ## init
            dt <- Terminal/division
            X <- xinit
            Xt <- xinit
            k <- numeric(d.size)
            k1 <- numeric(d.size)
            k2 <- numeric(d.size)
            k3 <- numeric(d.size)
            k4 <- numeric(d.size)
            assign(pars,0)  ## epsilon==0
            ## runge kutta
            for(t in 1:division){
              ## k1
              for(i in 1:d.size){
                assign(modelstate[i],Xt[i])
              }
              for(i in 1:d.size){
                k1[i] <- dt*eval(V0[i])
              }
              ## k2
              for(i in 1:d.size){
                assign(modelstate[i],Xt[i]+k1[i]/2)
              }
              for(i in 1:d.size){
                k2[i] <- dt*eval(V0[i])
              }
              ## k3
              for(i in 1:d.size){
                assign(modelstate[i],Xt[i]+k2[i]/2)
              }
              for(i in 1:d.size){
                k3[i] <- dt*eval(V0[i])
              }
              ## k4
              for(i in 1:d.size){
                assign(modelstate[i],Xt[i]+k3[i])
              }
              for(i in 1:d.size){
                k4[i] <- dt*eval(V0[i])
              }
              ## V0(X(t+dt))
              k <- (k1+k2+k2+k3+k3+k4)/6
              Xt <- Xt+k
              X <- rbind(X,Xt)
            }
            ## return matrix : (division+1)*d
            rownames(X) <- NULL
            colnames(X) <- modelstate
            return(ts(X,deltat=dt,start=0))
          }

# simFunctional
# public function to calculate Fe in (13.2).
setGeneric("simFunctional",
           function(yuima)
           standardGeneric("simFunctional")
           )
setMethod("simFunctional", signature(yuima="yuima"),
          function(yuima){
		    Xlen <- length(yuima@data)
			if(sum(Xlen != mean(Xlen)) != 0) {
			  yuima.warn("All length must be same yet.")
			  return(NULL)
			}
            if( (Xlen[1]-1) != yuima@sampling@n){
              yuima.warn("Length of time series and division do not much.")
              return(NULL)
            }
            
            e <- gete(yuima@functional)
            
            Fe <- funcFe.(yuima,as.matrix(onezoo(yuima)),e)
            return(Fe)
          })

# F0
# public function to calculate Fe at e=0
setGeneric("F0",
           function(yuima)
           standardGeneric("F0")
           )
setMethod("F0", signature(yuima="yuima"),
          function(yuima){
            
            X.t0 <- Get.X.t0(yuima)
            F0 <- funcFe.(yuima, X.t0, 0)
            return(F0)
          })

# Fnorm
# public function to calculate (Fe-F0)/e
setGeneric("Fnorm",
           function(yuima)
           standardGeneric("Fnorm")
           )
setMethod("Fnorm", signature(yuima="yuima"),
          function(yuima){
            e <- gete(yuima@functional)
            Fe <- simFunctional(yuima)
            F0 <- F0(yuima)
            return((Fe-F0)/e)
          })
