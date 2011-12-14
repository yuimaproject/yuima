#Estimating local Hölder exponent with the constant.


qgv<- function(yuima, filter.type="Daubechies", order=2, with.constant.est=TRUE, a=NULL, ...){

	call <- match.call()
	
	if(missing(yuima)){
		yuima.stop("yuima object is missing.")
	}

#Test for order to be integer
	
	if (missing(a)){

	if (filter.type=="Daubechies"){
		if (order==2){
			
		a<-1/sqrt(2)*c(.4829629131445341,
						  -.8365163037378077,
						   .2241438680420134,
						   .1294095225512603)
		}else{
			yuima.warn("There is no such order Daubechies' filter implemented, order have been fixed to 2.")	
				
		a<-1/sqrt(2)*c(.4829629131445341,
						  -.8365163037378077,
						   .2241438680420134,
						   .1294095225512603)
		
		
		} 
	}else if (filter.type=="Classical"){
		mesh<-0:order
		a=(-1)^(1-mesh)/2^order*choose(order,mesh)
	}
	}
	

	L<-length(a)
	a2<-rep(0,2*L)
	a2[seq(1,2*L,by=2)]<-a
	
	process<-yuima@data@zoo.data[[1]]
	N<-length(process)
	
	
	#Computation of the generalized quadratic variations
	
    V1<-sum(filter(process,a)^2,na.rm=TRUE)
	V2<-sum(filter(process,a2)^2,na.rm=TRUE)
	
	H<-1/2*log2(V2/V1)
	
	if (H>1){H<-0.999} #Over-estimation of H
	if (H<0){H<-0.001} 
	
	
	#Compute the estimation of the constant C.
	
	delta<-yuima@sampling@delta
	
    constfilt<-sum(a%*%t(a)*abs(matrix(0:(L-1),L,L)-matrix(0:(L-1),L,L,byrow=TRUE))^(2*H))
	
	C<- sqrt(-2*V1/(N-L)/(constfilt*delta^(2*H)))

	return(list(C=C,H=H))		 
}

