#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector residualCpp(NumericVector dx, NumericVector a, NumericVector b, double w, double h){
  int n = dx.length();
  NumericVector residual (n);
  for(int i=0; i<n; i++){
    residual(i) += (dx(i)-w*a(i))/(sqrt(h)*b(i));
  }
  return residual;
}






