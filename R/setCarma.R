# Carma_Model<-setClass("Carma_Model",
#                       slots = c(Cogarch_Model_Log="logical",
#                                 Under_Lev="yuima.model"),
#                       prototype=list(Cogarch_Model_Log = FALSE,
#                                      Under_Lev=NULL),
#                       contains= "yuima")
#ma.par, ar.par, lin.par=NULL  is.null
#state Variable consistente con Yuima.????
# Inserire Solve Variable???
# ... che mi serve per passare gli stessi parametri di setModel
# per i ... guardare qmle
# call<-matchcall()
# mydots guardare
setCarma<-function(p,q,loc.par=NULL,ar.par="beta",ma.par="alpha",lin.par=NULL,Carma.var="v",Latent.var="x", ...){
# We use the same parametrization as in Brockwell[2000] 
  
# mydots$Carma.var= V
# mydots$Latent.var= X ?????  
  # The Carma model is defined by two equations:  
  # \begin{eqnarray} 
  # V_{t}&=&\alpha_{0}+a'Y_{t-}\\
  # dY_{t}&=& BY_{t-}dt+ e dL\\
  # \end{eqnarray}  
  # p is the number of the moving average parameters \alpha
  # q  is the number of the autoregressive coefficient \beta
                     
#Default parameters
  
  call <- match.call()
  quadratic_variation<-FALSE
    
  mydots <- as.list(call)[-(1:2)]
  

#   hurst<-0.5 
#   jump.coeff <- character()
#   measure <- list() 
#   measure.type <- character() 
#   state.variable <- "Y"
#   jump.variable <- "z" 
#   time.variable <- "t"
#   mydots$xinit<- NULL
  if (is.null(mydots$hurst)){
    mydots$hurst<-0.5
  }
  
  if(is.null(mydots$time.variable)){
    mydots$time.variable<-"t"
  }

  if(is.null(mydots$jump.variable)){
    mydots$jump.variable<-"z"
  }
  
#   if(is.null(mydots$Carma.var)){
#     Carma.var<-"V"
#   } else{
#     Carma.var<-mydots$Carma.var
#   }
  
#   if(is.null(mydots$Latent.var)){
#     Latent.var<-"X"
#   } else{
#     Latent.var<-mydots$Latent.var
#   }
  
  if(is.null(mydots$xinit)){ 
    if(is.null(mydots$XinExpr)){
      mydots$xinit<-as.character(0*c(1:p))
    }else{  
      if(mydots$XinExpr==TRUE){
        Int.Var<-paste(Latent.var,"0",sep="")
        mydots$xinit<-paste(Int.Var,c(0:(p-1)),sep="")
      }
    }
  } else{
    dummy<-as.character(mydots$xinit)
    mydots$xinit<-dummy[-1]
  }
  
  if(p<q){
    yuima.stop("order of AR must be larger than MA order")
  }
    
  beta_coeff0<-paste("-",ar.par,sep="")
  beta_coeff<-paste(beta_coeff0,p:1,sep="")
  coeff_alpha<-c(paste(ma.par,0:q,sep=""),as.character(matrix(0,1,p-q-1)))
  fin_alp<-length(coeff_alpha)
  # We built the drift condition   
  
  Y_coeff<-paste(Latent.var,0:(p-1),sep="")
  fin_Y<-length(Y_coeff)
  V1<-paste(coeff_alpha,Y_coeff,sep="*")
  V2<-paste(V1,collapse="+")
#   alpha0<-paste(ma.par,0,sep="")
  if(is.null(loc.par)){
    V<-paste("(",V2,")",collapse="")
  } else {
    Vt<-paste(loc.par,V2,sep="+") 
    V<-paste("(",Vt,")",collapse="")
  }
  drift_last_cond<-paste(paste(beta_coeff,Y_coeff,sep="*"),collapse="")
  # Drift condition for the dV_{t}   
  
  drift_first_cond_1<-c(paste(coeff_alpha[-fin_alp],Y_coeff[-1],sep="*"))
  drift_first_cond_2<-paste(drift_first_cond_1,collapse="+")
  drift_first_cond_a<-paste("(",drift_last_cond,")",sep="")
  drift_first_cond_b<-paste(coeff_alpha[fin_alp],drift_first_cond_a,sep="*")
  drift_first_cond<-paste(drift_first_cond_2,drift_first_cond_b,sep="+")
  
  if(length(Y_coeff)>1)
  {drift_Carma<-c(drift_first_cond,Y_coeff[2:length(Y_coeff)],drift_last_cond)}else 
  {drift_Carma<-c(drift_first_cond,drift_last_cond)}
  # We need to consider three different situations   
  
  if(is.null(mydots$jump.coeff) & is.null(mydots$measure)  &  
        is.null(mydots$measure.type) & quadratic_variation==FALSE){
    # The Carma model is driven by a Brwonian motion
    if (is.null(lin.par)){ 
      diffusion_Carma<-matrix(c(coeff_alpha[fin_alp],as.character(matrix(0,(p-1),1)),"1"),(p+1),1)
      #     Latent.var<-Y_coeff
      Model_Carma1<-setModel(drift=drift_Carma, 
                             diffusion=diffusion_Carma,
                             hurst=mydots$hurst, 
                             state.variable=c(Carma.var,Y_coeff),  
                             solve.variable=c(Carma.var,Y_coeff),
                             xinit=c(V,mydots$xinit))
      if(length(Model_Carma1)==0){
        stop("Yuima model was not built") 
      } else { 
        return(Model_Carma1)
      }
    } else{
      if(ma.par==lin.par){
        first_term<-paste(coeff_alpha[fin_alp],V,sep="*")
        diffusion_Carma<-matrix(c(first_term,as.character(matrix(0,(p-1),1)),V),(p+1),1)  
        Model_Carma1<-setModel(drift=drift_Carma, 
                               diffusion=diffusion_Carma,
                               hurst=mydots$hurst, 
                               state.variable=c(Carma.var,Y_coeff),  
                               solve.variable=c(Carma.var,Y_coeff),
                               xinit=c(V,mydots$xinit))
        return(Model_Carma1)
      }else{ 
#         coeff_gamma<-c(paste(lin.par,1:p,sep=""),as.character(matrix(0,1,p-q)))
        coeff_gamma<-c(paste(lin.par,1:p,sep=""))
        Gamma1<-paste(coeff_gamma,Y_coeff,sep="*")
        Gamma2<-paste(Gamma1,collapse="+")
        gamma0<-paste(lin.par,0,sep="")
        Gammat<-paste(gamma0,Gamma2,sep="+") 
        Gamma<-paste("(",Gammat,")",collapse="")
        first_term<-paste(coeff_alpha[fin_alp],Gamma,sep="*")
        
        diffusion_Carma<-matrix(c(first_term,as.character(matrix(0,(p-1),1)),Gamma),(p+1),1)  
        Model_Carma1<-setModel(drift=drift_Carma, 
                               diffusion=diffusion_Carma,
                               hurst=mydots$hurst, 
                               state.variable=c(Carma.var,Y_coeff),  
                               solve.variable=c(Carma.var,Y_coeff),
                               xinit=c(V,mydots$xinit))
        
        return(Model_Carma1)
      }
      
    }                      
  } else {
    if( is.null(mydots$jump.coeff) & is.null(mydots$measure)  &  
          is.null(mydots$measure.type) & is.null(lin.par) & 
          quadratic_variation==TRUE){
      
      stop("The CoGarch model requires a Carma process driven by the discrete part of the quadratic covariation:
           You Must specify the Levy Measure")
      
    } else {
      
      if(quadratic_variation==FALSE & is.null(lin.par)){
        #         warning("At the moment, we need specify the underlying L?vy directly")
        #         diffusion_Carma<-matrix(c(coeff_alpha[fin_alp],as.character(matrix(0,(q-1),1)),"1"),(q+1),1)
        #         Model_Carma1<-setModel(drift=drift_Carma, diffusion=diffusion_Carma,
        #                                hurst=hurst, state.variable=c(V,Y_coeff),  solve.variable=c(V,Y_coeff))
        #         under_Lev1<-setModel(drift="0",diffusion="0",jump.coeff="1" ,
        #                              measure=measure ,measure.type=measure.type , 
        #                              jump.variable=jump.variable , time.variable=time.variable)
        #         if(length(Model_Carma1)==0){
        #           stop("Yuima model was not built") 
        #         } else {
        #           Model_Carma<-Carma_Model()
        #           Model_Carma@model <- Model_Carma1
        #           Model_Carma@Cogarch_Model_Log <- Cogarch_Model 
        #           Model_Carma@Under_Lev <-under_Lev1
        #           return(Model_Carma)
        #         }
        
        # LM 27/09 We use a modified 
        # setModel that allows us to build a sde where the slot model@jump.coeff is an vector
        
        # jump_Carma<-matrix(c(coeff_alpha[fin_alp],as.character(matrix(0,(q-1),1)),"1"),(q+1),1)
        jump_Carma<-c(coeff_alpha[fin_alp],as.character(matrix(0,(p-1),1)),"1")
        Model_Carma<-setModel(drift=drift_Carma, 
                              diffusion = NULL, 
                              hurst=mydots$hurst, 
                              jump.coeff=jump_Carma,
                              measure=eval(mydots$measure),
                              measure.type=mydots$measure.type, 
                              jump.variable=mydots$jump.variable, 
                              time.variable=mydots$time.variable,
                              state.variable=c(Carma.var,Y_coeff),  
                              solve.variable=c(Carma.var,Y_coeff),
                              xinit=c(V,mydots$xinit))
        return(Model_Carma)
      } else {
        if (quadratic_variation==FALSE ){
        # Selecting Quadratic_Variation==FALSE and specifying the Heteroskedatic.param in the model, 
        # The coefficient of the error term is a vector in which the last element is an affine linear function 
        # of the vector space Y               
        
        # We have to consider two different cases:
        # a) The last component of the error term is  $V_{t-}=\alpha_{0}+a'Y_{t-}$. Usually 
        #    this condition is used for the definition of the COGARCH(p,q) introduced in Brockwell and Davis     and 
        # b) The last component of the error term is a linear function of the state variable $Y_{t}$
        #     different of the observable variable V.  
        if(ma.par==lin.par){
          jump_first_comp<-paste(coeff_alpha[fin_alp],V,sep="*")
          jump_Carma<-c(jump_first_comp,as.character(matrix(0,(p-1),1)),V)
        }else{
#           coeff_gamma<-c(paste(lin.par,1:p,sep=""),as.character(matrix(0,1,q-p)))
          coeff_gamma<-c(paste(lin.par,1:p,sep=""))
          Gamma1<-paste(coeff_gamma,Y_coeff,sep="*")
          Gamma2<-paste(Gamma1,collapse="+")
          gamma0<-paste(lin.par,0,sep="")
          Gammat<-paste(gamma0,Gamma2,sep="+") 
          Gamma<-paste("(",Gammat,")",collapse="")
          jump_first_comp<-paste(coeff_alpha[fin_alp],Gamma,sep="*")
          jump_Carma<-c(jump_first_comp,as.character(matrix(0,(p-1),1)),Gamma)
        }
        
#         Model_Carma<-setModel(drift=drift_Carma,
#                               diffusion =NULL,
#                               hurst=0.5,
#                               jump.coeff=jump_Carma,
#                               measure=eval(mydots$measure),
#                               measure.type=mydots$measure.type,
#                               jump.variable=mydots$jump.variable,
#                               time.variable=mydots$time.variable,
#                               state.variable=c(Carma.var,Y_coeff),  
#                               solve.variable=c(Carma.var,Y_coeff),
#                               c(V,mydots$xinit))
#         return(Model_Carma)
        }
        
        Model_Carma<-setModel(drift=drift_Carma,
                              diffusion =NULL,
                              hurst=mydots$hurst,
                              jump.coeff=jump_Carma,
                              measure=eval(mydots$measure),
                              measure.type=mydots$measure.type,
                              jump.variable=mydots$jump.variable,
                              time.variable=mydots$time.variable,
                              state.variable=c(Carma.var,Y_coeff),  
                              solve.variable=c(Carma.var,Y_coeff),
                              c(V,mydots$xinit))
        return(Model_Carma)
         if(quadratic_variation==TRUE){
#           
                stop("Work in Progress: Implementation of CARMA model for CoGarch. 
                     We need the increments of Quadratic Variation")
#           
#           diffusion_first_cond<-paste(coeff_alpha[fin_alp],V,sep="*")
#           diffusion_Carma<-matrix(c(diffusion_first_cond,as.character(matrix(0,(q-1),1)),Vt),(q+1),1)
#           #    At the present time, Yuima does not support Multi - Jumps 
#           Model_Carma1<-setModel(drift=drift_Carma, diffusion=diffusion_Carma,
#                                  hurst=hurst, state.variable=c(V,Y_coeff),  solve.variable=c(V,Y_coeff))
#           under_Lev1<-setModel(drift="0",diffusion="0",jump.coeff="1" ,
#                                measure=measure ,measure.type=measure.type , 
#                                jump.variable=jump.variable , time.variable=time.variable)
#           if(length(Model_Carma1)==0){
#             stop("Yuima model was not built") 
#           } else {
#             Model_Carma<-Carma_Model()
#             Model_Carma@model <- Model_Carma1
#             Model_Carma@Cogarch_Model_Log <- Cogarch_Model 
#             Model_Carma@Under_Lev <-under_Lev1
#             return(Model_Carma)
#           }
#           
         } 
      } 
    }
    
  }  
}
