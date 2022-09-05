#include <RcppArmadillo.h>

using namespace Rcpp;

/*
Inputs
int lengthObs
arma::mat B, arma::mat Btilde, arma::mat InvBtilde
double a0, double bq, double a1, double V, double PseudologLik
arma::rowvec ta

arma::colvec state, arma::colvec stateMean, arma::colvec e
arma::vec DeltaG2, arma::vec Deltat


double Irregular_PseudoLoglik_COG(int lengthObs,  arma::mat B, arma::mat Btilde, arma::mat InvBtilde,
                                  double a0, double bq, double a1, double V, double PseudologLik,arma::rowvec ta,
                                  arma::colvec state, arma::colvec stateMean, arma::colvec e,
                                  arma::vec DeltaG2, arma::vec Deltat) {
  double res=0;
  double VarDeltaG = 0;
  arma::mat I = arma::eye<arma::mat>(B.n_rows, B.n_rows);
  // for(i in c(1:(lengthObs))){
  for(int i=0;i<lengthObs;i++){
    arma::mat DeltatB1 = expmat(B * Deltat[i]);
    arma::mat DeltatB2 = expmat(Btilde * Deltat[i]);
    VarDeltaG = a0 * arma::as_scalar(Deltat[i]) * bq / (bq - a1); //+ arma::as_scalar(ta * DeltatB2 * InvBtilde * (I - DeltatB2) * (state - stateMean));
    state <- (I + DeltaG2[i] / V * e * ta) * DeltatB1 * state+a0 * DeltaG2[i] / V * e;
    if(VarDeltaG>0){
      res =res - 0.5*(arma::as_scalar((DeltaG2[i] / VarDeltaG )) + log(VarDeltaG) + log(2 * 3.141593));
    //  res = res - 1/2 * (arma::as_scalar(DeltaG2[i] / VarDeltaG) + log(VarDeltaG) + log(2 * 3.141593));
    }else{
      res = res - 1000000.6543;
    }
  }
  return PseudologLik;
}
*/

// [[Rcpp::export]]
double Irregular_PseudoLoglik_COG(int lengthObs,  arma::mat B, arma::mat Btilde, arma::mat InvBtilde,
             double a0, double bq, double a1, double V, double PseudologLik,arma::rowvec ta,
              arma::colvec state, arma::colvec stateMean, arma::colvec e,
              arma::vec DeltaG2, arma::vec Deltat) {
  double res = 0;
  double VarDeltaG = 0;
  double penal = 1000000000;
  arma::mat I = arma::eye<arma::mat>(B.n_rows, B.n_rows);
// for(i in c(1:(lengthObs))){
  for(int i=0;i<lengthObs;i++){
    arma::mat DeltatB1 = expmat(B * Deltat[i]);
    arma::mat DeltatB2 = expmat(Btilde * Deltat[i]);
    VarDeltaG = a0 * Deltat[i] * bq / (bq - a1) + arma::as_scalar(ta * DeltatB2 * InvBtilde * (I - DeltatB2) * (state - stateMean));
    state = (I + DeltaG2[i] / V * e * ta) * DeltatB1 * state+a0 * DeltaG2[i] / V * e;
    V = arma::as_scalar(a0 + ta*state);
  if(VarDeltaG>0){
    res = res - 0.5*(DeltaG2[i] / VarDeltaG  + log(VarDeltaG) + log(2 * 3.141593));
  }else{
    res = res-penal;
  }
}
return res;
}


/*
for(i in c(1:(lengthObs))){
DeltatB1 <- expm(B*Deltat[i])
DeltatB2 <- expm(Btilde*Deltat[i])
VarDeltaG <- as.numeric(a0*Deltat[i]*bq/(bq-a1)+ta%*%DeltatB2%*%InvBtilde%*%(I-DeltatB2)%*%(state-stateMean))
state <- (I+DeltaG2[i]/V*e%*%ta)%*%DeltatB1%*%state+a0*DeltaG2[i]/V*e
V <- as.numeric(a0+ta%*% state)
if(VarDeltaG>0){
PseudologLik<--1/2*(DeltaG2[i]/VarDeltaG+log(VarDeltaG)+log(2*pi))+PseudologLik
}else{
PseudologLik<-PseudologLik-10^6
}


}
*/

