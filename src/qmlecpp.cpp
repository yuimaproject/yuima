#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
double detcpp(NumericMatrix A){   //det(A)
    int n = A.nrow();
    double det = 1.0,buf;
    NumericMatrix B = clone(A);
    
    for(int i=0;i<n;i++){
        buf = 1/B(i,i);
        for(int j=i+1;j<n;j++){
            for(int k=i+1;k<n;k++){
                B(j,k)-=B(i,k)*B(j,i)*buf;
            }
        }
        det*=B(i,i);
    }
    return det;
}


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
// [[Rcpp::export]]
NumericMatrix solvecpp(NumericMatrix A){ //solve(A)
    int n = A.ncol();
    double buf;
    NumericMatrix B = clone(A);
    NumericMatrix inv(n,n);
    
    for(int i=0;i<n;i++){
        inv(i,i) = 1.0;
        buf = 1.0/B(i,i);
        for(int j=0;j<n;j++){
            if(j<i+1) inv(i,j)*=buf;
            else B(i,j)*=buf;
        }
        for(int j=0;j<n;j++){
            if(i!=j){
                buf=B(j,i);
                for(int k=0;k<n;k++){
                    if(k<i+1) inv(j,k)-=inv(i,k)*buf;
                    else B(j,k)-=B(i,k)*buf;
                }
            }
        }
    }
    return inv;
}

// [[Rcpp::export]]
double sub_f(NumericMatrix S,NumericVector b){   // tr(S%*%tcrossprod(b))
    int n = S.nrow();
    double tr = 0;
    
    for(int i=0;i<n;i++){
        for(int k=0;k<n;k++){
            tr += S(i,k)*b(k)*b(i);
        }
    }
    return tr;
}
// [[Rcpp::export]]
double likndim(NumericMatrix dx,NumericMatrix b,NumericMatrix A,double h){
    int n = dx.nrow();
    int m = dx.ncol();
    double tmp1 = 0;double tmp2 = 0;
    NumericMatrix S(m,m);
    
    for(int i=0;i<n;i++){
        S = Smake(A(i,_),m);
        tmp1 += log(detcpp(S));
        tmp2 += sub_f(solvecpp(S),dx(i,_)-h*b(i,_));
    }
    return tmp1+tmp2/h;
}
