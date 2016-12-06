#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

/*
// [[Rcpp::export]]
NumericMatrix Smake(NumericVector b,int d){   //tcrossprod(matrix(b,d,r))
int r = b.length()/d;
NumericMatrix S(d,d);

for(int i=0;i<d;i++){
for(int j=0;j<d;j++){
for(int k=0;k<r;k++){
S(i,j) += b(d*k+i)*b(d*k+j);
}
}
}
return S;
}
*/




// [[Rcpp::export]]
double W1(NumericMatrix crossdx,NumericMatrix b,NumericMatrix A,double h){
  int n = b.nrow();
  int m = b.ncol();
  int r = A.ncol()/m;
  double tmp1 = 0;
  NumericMatrix S(m,m);
  
  for(int i=0;i<n;i++){
    //S = Smake(A(i,_),m);
    for(int j=0;j<m;j++){
      for(int k=0;k<m;k++){
        for(int l=0;l<r;l++){
          S(j,k) += A(i,m*l+j)*A(i,m*l+k);
        }
        tmp1 += (crossdx(i,m*j+k)-h*S(j,k))*(crossdx(i,m*j+k)-h*S(j,k));
        S(j,k) = 0;
      }
    }
  }
  return tmp1;  //in "minusquasilogl_W1",QL <- W1()*(-0.5*h*h)
}

// [[Rcpp::export]]
double W2(NumericMatrix dx,NumericMatrix b,double h){
  int n = dx.nrow();
  int m = dx.ncol();
  double tmp1 = 0;
  
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      tmp1 += (dx(i,j)-h*b(i,j))*(dx(i,j)-h*b(i,j));
    }
  }
  return tmp1; //in "minusquasilogl_W2",QL <- W2()*(-0.5*h)
}
