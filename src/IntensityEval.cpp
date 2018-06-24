#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector evalKernelCpp(StringMatrix Integrand2, 
                             ExpressionVector Integrand2expr, 
                             Environment myenvd1, Environment myenvd2,
                             LogicalVector ExistdN, LogicalVector ExistdX,
                             NumericVector gridTime, IntegerVector dimCol,
                             StringVector NameCol, StringVector JumpTimeName){
  int FinT = gridTime.size();
  NumericVector GlobalKernel(FinT);
  for(int t=0;t<FinT;t++){
      if(ExistdN[0]){
        String ttime = myenvd1.get("t.time");
        myenvd1.assign(ttime, gridTime[t]);
      }
      if(ExistdX[0]){
        String ttimedX = myenvd2.get("t.time");
        myenvd2.assign(ttimedX,gridTime[t]);
      }
      double IntegralKernel = 0;
      for(int i=0;i<dimCol[0];i++){
        if(ExistdN[0]){
          StringVector namedJumpTimeX1 = myenvd1.get("namedJumpTimeX"); 
          StringVector namedX1 = myenvd1.get("namedX");
          //  printf(namedJumpTimeX1[0]);
          for(int j=0;j<namedJumpTimeX1.size();j++){
            if(JumpTimeName[i]==namedJumpTimeX1[j]){
              String dummyString = as<std::string>(JumpTimeName[i]);
              NumericVector JumpTimeDiff = myenvd1.get(dummyString);
              String vartime =myenvd1.get("var.time");
              myenvd1.assign(vartime,JumpTimeDiff);
            }
          }
          for(int hh=0;hh<NameCol.size();hh++){
            if(namedX1[i]==NameCol[hh]){
              Language eval_call( "eval", Integrand2expr[hh], myenvd1);
              SEXP X = Rf_eval( eval_call, myenvd1);
              NumericVector XX(X);
              NumericVector dumY = na_omit(XX);
              IntegralKernel = IntegralKernel+sum(dumY);
            }
          }
/*
          if(any(cond)){
          assign(my.envd1$var.time,my.envd1[[my.envd1$namedJumpTimeX[cond]]],envir=my.envd1)
          condpos <- my.envd1$namedX %in% NameCol[i]  
          if(any(condpos)){
          IntegralKernelDum<- sum(eval(Integrand2expr[condpos], envir=my.envd1))
          IntegralKernel<-IntegralKernel+IntegralKernelDum
          }
          }
          */ 
        }
        if(ExistdX[0]){
          StringVector namedJumpTimeX2 = myenvd2.get("namedJumpTimeX");
          StringVector namedX2 = myenvd2.get("namedX");
          for(int j=0;j<namedJumpTimeX2.size();j++){
            if(JumpTimeName[i]==namedJumpTimeX2[j]){
              String dummyString2 = as<std::string>(JumpTimeName[i]);
              NumericVector JumpTimeDiff2 = myenvd2.get(dummyString2);
              String vartime2 =myenvd2.get("var.time");
              myenvd2.assign(vartime2,JumpTimeDiff2);
            }
          }
          for(int hh=0;hh<NameCol.size();hh++){
            if(namedX2[i]==NameCol[hh]){
              Language eval_call( "eval", Integrand2expr[hh], myenvd2);
              SEXP X2 = Rf_eval( eval_call, myenvd2);
              NumericVector XX2(X2);
              NumericVector dumY2 = na_omit(XX2);
              IntegralKernel = IntegralKernel+sum(dumY2);
            }
          }
          
        }
/*
        if(ExistdX){  
          cond <- my.envd2$namedJumpTimeX %in% JumpTimeName[i]
          if(any(cond)){
            assign(my.envd2$var.time,my.envd2[[my.envd2$namedJumpTimeX[cond]]],envir=my.envd2)
            condpos <- my.envd2$namedX %in% NameCol[i]
            if(any(condpos)){
              IntegralKernelDum<- sum(eval(Integrand2expr[condpos], envir=my.envd2))
              IntegralKernel<-IntegralKernel+IntegralKernelDum
            }
          }
        }
*/
    }
      GlobalKernel[t] = IntegralKernel;
  }
  return GlobalKernel;
}