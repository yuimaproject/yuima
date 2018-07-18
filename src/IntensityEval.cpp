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
              IntegralKernel = IntegralKernel+sum(XX);
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
// [[Rcpp::export]]
NumericVector evalKernelCpp2(StringMatrix Integrand2, 
                             ExpressionVector Integrand2expr, 
                             Environment myenvd1, Environment myenvd2, LogicalVector CondIntensity,
                             StringVector NameCountingVar,
                             StringVector Namecovariates,
                             LogicalVector ExistdN, LogicalVector ExistdX,
                             List gridTime, IntegerVector dimCol,
                             StringVector NameCol, StringVector JumpTimeName){
  int FinT = gridTime.size();
  int lengtcov= Namecovariates.size();
  NumericVector GlobalKernel(FinT);
  for(int t=0;t<FinT;t++){
    double IntegralKernel = 0;
    if(ExistdN[0]){
      String ttime = myenvd1.get("t.time");
      NumericVector PosInTimedN = gridTime[t];
      myenvd1.assign(ttime, PosInTimedN[0]);
    }
    if(CondIntensity[0]){
      StringVector PosListCountingVariableC = myenvd1.get("PosListCountingVariable");
      int dimCountingVar = NameCountingVar.size();
      for(int i=0;i<dimCountingVar;i++){
        String dummyStringListCount = as<std::string>(PosListCountingVariableC[i]);
        List DummyListCount = myenvd1.get(dummyStringListCount);
        String dummyStringCountN = as<std::string>(NameCountingVar[i]);
        myenvd1.assign(dummyStringCountN,DummyListCount[t]);
      }
    }
    if(lengtcov>0){
      StringVector PosListCovariatesC = myenvd1.get("PosListCovariates");
      for(int h=0;h<lengtcov;h++){
        String dummyStringListcovariates = as<std::string>(PosListCovariatesC[h]);
        List DummyListcovariates = myenvd1.get(dummyStringListcovariates);
        String dummyStringcovariatesN = as<std::string>(Namecovariates[h]);
        myenvd1.assign(dummyStringcovariatesN,DummyListcovariates[t]);
      }
    }
    
    for(int i=0;i<dimCol[0];i++){  
      if(ExistdN[0]){
        StringVector namedJumpTimeX1 = myenvd1.get("namedJumpTimeX"); 
        StringVector namedX1 = myenvd1.get("namedX");
        for(int j=0;j<namedJumpTimeX1.size();j++){
          if(JumpTimeName[i]==namedJumpTimeX1[j]){
            String dummyString = as<std::string>(JumpTimeName[i]);
            List JumpTimeDiff = myenvd1.get(dummyString);
            String vartime =myenvd1.get("var.time");
            myenvd1.assign(vartime,JumpTimeDiff[t]);
            String dumdN = NameCol[i];
            List dNAll = myenvd1.get(dumdN);
            NumericVector dN = dNAll[t];
            if(dN.size()==0){
              NumericVector dN[0] = {0};
            }
            Language eval_call( "eval", Integrand2expr[i], myenvd1);
            SEXP X = Rf_eval( eval_call, myenvd1);
            NumericVector XX(X);
            NumericVector dumY = na_omit(XX);
            double dumSumKer = 0;
            for(int l=0;l<dumY.size();l++){
              dumSumKer = dumSumKer + dumY[l]*dN[l];
            }
            IntegralKernel = IntegralKernel+dumSumKer;            
          }
        }
      }
      if(ExistdX[0]){
        StringVector namedJumpTimeX2 = myenvd2.get("namedJumpTimeX");
        for(int j=0;j<namedJumpTimeX2.size();j++){
          if(JumpTimeName[i]==namedJumpTimeX2[j]){
            String ttimedX = myenvd2.get("t.time");
            NumericVector PosInTimedX = gridTime[t];
            myenvd2.assign(ttimedX,PosInTimedX[0]);
            String dummyString2 = as<std::string>(JumpTimeName[i]);
            NumericVector JumpTime2 = myenvd2.get(dummyString2);
            String vartime2 =myenvd2.get("var.time");
            IntegerVector idx = seq(0,t);
            myenvd2.assign(vartime2,JumpTime2[idx]);
            Language eval_call( "eval", Integrand2expr[i], myenvd2);
            SEXP X2 = Rf_eval( eval_call, myenvd2);
            NumericVector XX2(X2);
            
            String dummyString3 = as<std::string>(NameCol[i]);
            NumericVector dXdummy = myenvd2.get(dummyString3);
            NumericVector dumY1 = XX2[idx]*dXdummy[idx];
            NumericVector dumY2 = na_omit(dumY1);
            IntegralKernel = IntegralKernel+sum(dumY2);
            
          }
        }
      }
    }
    GlobalKernel[t] = IntegralKernel;
  }
  return GlobalKernel;
}