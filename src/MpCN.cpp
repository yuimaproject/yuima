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
double sqnorm(NumericVector x){
  double y=0;
  int d = x.length();
  for(int i=0;i<d;i++) y += x[i]*x[i];
  if(y==0){   //save NA risk of rgamma  
    y=0.00000000000001;
  }
  return y;
}
// [[Rcpp::export]]
NumericVector makeprop(NumericVector mu,NumericVector sample,
                       NumericVector low,NumericVector up){ 
  int d = mu.length();
  NumericVector prop(d);
  NumericVector tmp(d);
  double tmp2;double rho = 0.8;
  
  tmp = mu+sqrt(rho)*(sample-mu);
  tmp2 = 2.0/sqnorm(sample-mu);
  while(1){
    prop = tmp+rnorm(d)*sqrt((1.0-rho)/rgamma(1,0.5*d,tmp2)(0));
    if((sum(low>prop)+sum(up<prop))==0) break;
  }
  return prop;
}
