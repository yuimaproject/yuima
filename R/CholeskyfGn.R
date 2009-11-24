CholeskyfGn <-
function(mesh, H)
{

##--------------------------------------------------------------------------
## Author : Alexandre Brouste
## Project: Yuima and Fieldsim
##--------------------------------------------------------------------------

##--------------------------------------------------------------------------
## Input        :         mesh : mesh grid where the fBm is evaluated
##                        H    : self-similarity parameter          
##      		 
##
## Output       :         simulation of a standard fractional Brownian noise
##                        for the mesh grid evaluation by Choleki s 
##			  decomposition of the covariance matrix of the fGn.
##--------------------------------------------------------------------------

##---------------------------------------------------------
## Complexity O(N^3) via method chol
## Taille memoire N^2
## -------------------------------------------------------------------------
                     
	N<-length(mesh)-1	 #j'ai retirer 1    # N+1 is the size of the fGn sample to be simulated
matcov <- matrix(0,N+1,N+1) # Covariance matrix of the fGn
H2 <- 2*H
	
	for (i in (1:(N+1))) {
		j <- i:(N+1)
		matcov[i, j]<- 0.5 * (abs(mesh[i+1]-mesh[j])^H2 + abs(mesh[i]-mesh[j+1])^H2 - abs(mesh[i] - mesh[j])^H2-abs(mesh[i+1] - mesh[j+1])^H2)
		matcov[j, i]<- matcov[i, j]
	}
	L <- chol(matcov)
	Z <- rnorm(N+1)
	fGn <- t(L) %*% Z
	return(fGn)
}