sampling2grid <- function(samp){
	
#Check if there is already a grid	
	if (length(samp@grid)>0){
		
#It is the general case where there is a grid 
#including irregular deterministic case
		return(samp@grid)
	}
	
#Check if this grid is random	
	if (samp@random==FALSE){			
#Regular deterministic grid
		return(seq(samp@Initial,samp@Terminal,length=samp@division))
	}
	
	return(1) #To do: Random grid		
}

