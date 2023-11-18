#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix driftTermCpp(ExpressionVector drift, CharacterVector modelstate, NumericMatrix data, Environment env) {
    int dataNum = data.nrow();
    int stateNum = modelstate.length();
    int i, j;
    NumericMatrix res(dataNum, drift.length());
    for(i = 0; i < dataNum; i++) {
        for(j = 0; j < stateNum; j++) {
            env.assign(as<std::string>(modelstate[j]), data(i, j));
        }
        for(j = 0; j < drift.length(); j++) {
            SEXP evaluatedValue = Rf_eval(drift[j], env);
            res(i, j) = as<double>(evaluatedValue);
        }
    }
    return res;
}

// [[Rcpp::export]]
NumericVector diffusionTermCpp(List diffusion, CharacterVector modelstate, NumericMatrix data, Environment env) {
    int dataNum = data.nrow();
    int stateNum = modelstate.length();
    int i, j, k;
    int eqNum = diffusion.length();
    int wienerNum = (as<ExpressionVector>(diffusion[0])).length();
    NumericVector res(dataNum * eqNum * wienerNum);
    for(i = 0; i < dataNum; i++) {
        for(j = 0; j < stateNum; j++) {
            env.assign(as<std::string>(modelstate[j]), data(i, j));
        }
        for(j = 0; j < wienerNum; j++) {
            for(k = 0; k < eqNum; k++) {
                SEXP evaluatedValue = Rf_eval((as<ExpressionVector>(diffusion[k]))[j], env);
                res(i * eqNum * wienerNum + eqNum * j + k) = as<double>(evaluatedValue);
            }
        }
    }
    return res;
}

// [[Rcpp::export]]
NumericVector measureTermCpp(List measure, CharacterVector modelstate, NumericMatrix data, Environment env) {
    int dataNum = data.nrow();
    int stateNum = modelstate.length();
    int i, j, k;
    int eqNum = measure.length();
    int jumpNum = (as<ExpressionVector>(measure[0])).length();
    NumericVector res(dataNum * eqNum * jumpNum);
    for(i = 0; i < dataNum; i++) {
        for(j = 0; j < stateNum; j++) {
            env.assign(as<std::string>(modelstate[j]), data(i, j));
        }
        for(j = 0; j < jumpNum; j++) {
            for(k = 0; k < eqNum; k++) {
                SEXP evaluatedValue = Rf_eval((as<ExpressionVector>(measure[k]))[j], env);
                res(i * eqNum * jumpNum + eqNum * j + k) = as<double>(evaluatedValue);
            }
        }
    }
    return res;
}
