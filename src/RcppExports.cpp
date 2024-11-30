// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// evalKernelCpp
NumericVector evalKernelCpp(StringMatrix Integrand2, ExpressionVector Integrand2expr, Environment myenvd1, Environment myenvd2, LogicalVector ExistdN, LogicalVector ExistdX, NumericVector gridTime, IntegerVector dimCol, StringVector NameCol, StringVector JumpTimeName);
RcppExport SEXP _yuima_evalKernelCpp(SEXP Integrand2SEXP, SEXP Integrand2exprSEXP, SEXP myenvd1SEXP, SEXP myenvd2SEXP, SEXP ExistdNSEXP, SEXP ExistdXSEXP, SEXP gridTimeSEXP, SEXP dimColSEXP, SEXP NameColSEXP, SEXP JumpTimeNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringMatrix >::type Integrand2(Integrand2SEXP);
    Rcpp::traits::input_parameter< ExpressionVector >::type Integrand2expr(Integrand2exprSEXP);
    Rcpp::traits::input_parameter< Environment >::type myenvd1(myenvd1SEXP);
    Rcpp::traits::input_parameter< Environment >::type myenvd2(myenvd2SEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type ExistdN(ExistdNSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type ExistdX(ExistdXSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gridTime(gridTimeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dimCol(dimColSEXP);
    Rcpp::traits::input_parameter< StringVector >::type NameCol(NameColSEXP);
    Rcpp::traits::input_parameter< StringVector >::type JumpTimeName(JumpTimeNameSEXP);
    rcpp_result_gen = Rcpp::wrap(evalKernelCpp(Integrand2, Integrand2expr, myenvd1, myenvd2, ExistdN, ExistdX, gridTime, dimCol, NameCol, JumpTimeName));
    return rcpp_result_gen;
END_RCPP
}
// evalKernelCpp2
NumericVector evalKernelCpp2(StringMatrix Integrand2, ExpressionVector Integrand2expr, Environment myenvd1, Environment myenvd2, LogicalVector CondIntensity, StringVector NameCountingVar, StringVector Namecovariates, LogicalVector ExistdN, LogicalVector ExistdX, List gridTime, IntegerVector dimCol, StringVector NameCol, StringVector JumpTimeName);
RcppExport SEXP _yuima_evalKernelCpp2(SEXP Integrand2SEXP, SEXP Integrand2exprSEXP, SEXP myenvd1SEXP, SEXP myenvd2SEXP, SEXP CondIntensitySEXP, SEXP NameCountingVarSEXP, SEXP NamecovariatesSEXP, SEXP ExistdNSEXP, SEXP ExistdXSEXP, SEXP gridTimeSEXP, SEXP dimColSEXP, SEXP NameColSEXP, SEXP JumpTimeNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringMatrix >::type Integrand2(Integrand2SEXP);
    Rcpp::traits::input_parameter< ExpressionVector >::type Integrand2expr(Integrand2exprSEXP);
    Rcpp::traits::input_parameter< Environment >::type myenvd1(myenvd1SEXP);
    Rcpp::traits::input_parameter< Environment >::type myenvd2(myenvd2SEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type CondIntensity(CondIntensitySEXP);
    Rcpp::traits::input_parameter< StringVector >::type NameCountingVar(NameCountingVarSEXP);
    Rcpp::traits::input_parameter< StringVector >::type Namecovariates(NamecovariatesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type ExistdN(ExistdNSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type ExistdX(ExistdXSEXP);
    Rcpp::traits::input_parameter< List >::type gridTime(gridTimeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dimCol(dimColSEXP);
    Rcpp::traits::input_parameter< StringVector >::type NameCol(NameColSEXP);
    Rcpp::traits::input_parameter< StringVector >::type JumpTimeName(JumpTimeNameSEXP);
    rcpp_result_gen = Rcpp::wrap(evalKernelCpp2(Integrand2, Integrand2expr, myenvd1, myenvd2, CondIntensity, NameCountingVar, Namecovariates, ExistdN, ExistdX, gridTime, dimCol, NameCol, JumpTimeName));
    return rcpp_result_gen;
END_RCPP
}
// sqnorm
double sqnorm(NumericVector x);
RcppExport SEXP _yuima_sqnorm(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sqnorm(x));
    return rcpp_result_gen;
END_RCPP
}
// makeprop
NumericVector makeprop(NumericVector mu, NumericVector sample, NumericVector low, NumericVector up);
RcppExport SEXP _yuima_makeprop(SEXP muSEXP, SEXP sampleSEXP, SEXP lowSEXP, SEXP upSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type low(lowSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type up(upSEXP);
    rcpp_result_gen = Rcpp::wrap(makeprop(mu, sample, low, up));
    return rcpp_result_gen;
END_RCPP
}
// is_zero
bool is_zero(std::string const& x);
RcppExport SEXP _yuima_is_zero(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string const& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(is_zero(x));
    return rcpp_result_gen;
END_RCPP
}
// cpp_split
std::vector< std::vector<std::string> > cpp_split(std::vector<std::string> const& str, std::string const sep);
RcppExport SEXP _yuima_cpp_split(SEXP strSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> const& >::type str(strSEXP);
    Rcpp::traits::input_parameter< std::string const >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_split(str, sep));
    return rcpp_result_gen;
END_RCPP
}
// cpp_paste
std::vector<std::string> cpp_paste(std::vector<std::string> const& x, std::vector<std::string> const& y, std::string const sep);
RcppExport SEXP _yuima_cpp_paste(SEXP xSEXP, SEXP ySEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::string const >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_paste(x, y, sep));
    return rcpp_result_gen;
END_RCPP
}
// cpp_collapse
std::string cpp_collapse(std::vector<std::string> const& str, std::string const sep);
RcppExport SEXP _yuima_cpp_collapse(SEXP strSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> const& >::type str(strSEXP);
    Rcpp::traits::input_parameter< std::string const >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_collapse(str, sep));
    return rcpp_result_gen;
END_RCPP
}
// cpp_outer
std::vector<std::string> cpp_outer(std::vector<std::string> const& x, std::vector<std::string> const& y);
RcppExport SEXP _yuima_cpp_outer(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> const& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_outer(x, y));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ito_outer
std::vector<std::string> cpp_ito_outer(std::vector<std::string> const& x, std::vector<std::string> const& y);
RcppExport SEXP _yuima_cpp_ito_outer(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> const& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ito_outer(x, y));
    return rcpp_result_gen;
END_RCPP
}
// cpp_label
std::string cpp_label(std::vector<int> I);
RcppExport SEXP _yuima_cpp_label(SEXP ISEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type I(ISEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_label(I));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ito_product
std::vector<std::string> cpp_ito_product(std::vector<int> const& idx, List const& dZ, List const& Z_K, std::vector<int> const& K, int d, int a, int p, int q);
RcppExport SEXP _yuima_cpp_ito_product(SEXP idxSEXP, SEXP dZSEXP, SEXP Z_KSEXP, SEXP KSEXP, SEXP dSEXP, SEXP aSEXP, SEXP pSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> const& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< List const& >::type dZ(dZSEXP);
    Rcpp::traits::input_parameter< List const& >::type Z_K(Z_KSEXP);
    Rcpp::traits::input_parameter< std::vector<int> const& >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ito_product(idx, dZ, Z_K, K, d, a, p, q));
    return rcpp_result_gen;
END_RCPP
}
// cpp_E
std::vector<std::string> cpp_E(std::vector<std::string> str);
RcppExport SEXP _yuima_cpp_E(SEXP strSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type str(strSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_E(str));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ito
List cpp_ito(List const& K_set, List const& dZ, List const& Z_K, int d, int r);
RcppExport SEXP _yuima_cpp_ito(SEXP K_setSEXP, SEXP dZSEXP, SEXP Z_KSEXP, SEXP dSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List const& >::type K_set(K_setSEXP);
    Rcpp::traits::input_parameter< List const& >::type dZ(dZSEXP);
    Rcpp::traits::input_parameter< List const& >::type Z_K(Z_KSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ito(K_set, dZ, Z_K, d, r));
    return rcpp_result_gen;
END_RCPP
}
// calc_filter_vcov
arma::cube calc_filter_vcov(arma::cube un_dr_sl, arma::cube un_diff, arma::cube ob_dr_sl, arma::cube inv_sq_ob_diff, arma::mat init, double delta);
RcppExport SEXP _yuima_calc_filter_vcov(SEXP un_dr_slSEXP, SEXP un_diffSEXP, SEXP ob_dr_slSEXP, SEXP inv_sq_ob_diffSEXP, SEXP initSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type un_dr_sl(un_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type un_diff(un_diffSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ob_dr_sl(ob_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type inv_sq_ob_diff(inv_sq_ob_diffSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type init(initSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_filter_vcov(un_dr_sl, un_diff, ob_dr_sl, inv_sq_ob_diff, init, delta));
    return rcpp_result_gen;
END_RCPP
}
// calc_filter_vcov_time_homogeneous
arma::cube calc_filter_vcov_time_homogeneous(arma::mat un_dr_sl, arma::mat un_diff, arma::mat ob_dr_sl, arma::mat inv_sq_ob_diff, arma::mat init, double delta, int n);
RcppExport SEXP _yuima_calc_filter_vcov_time_homogeneous(SEXP un_dr_slSEXP, SEXP un_diffSEXP, SEXP ob_dr_slSEXP, SEXP inv_sq_ob_diffSEXP, SEXP initSEXP, SEXP deltaSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type un_dr_sl(un_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type un_diff(un_diffSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ob_dr_sl(ob_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sq_ob_diff(inv_sq_ob_diffSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type init(initSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_filter_vcov_time_homogeneous(un_dr_sl, un_diff, ob_dr_sl, inv_sq_ob_diff, init, delta, n));
    return rcpp_result_gen;
END_RCPP
}
// calc_filter_vcov_are
arma::mat calc_filter_vcov_are(arma::mat un_dr_sl, arma::mat un_diff, arma::mat ob_dr_sl, arma::mat inv_sq_ob_diff);
RcppExport SEXP _yuima_calc_filter_vcov_are(SEXP un_dr_slSEXP, SEXP un_diffSEXP, SEXP ob_dr_slSEXP, SEXP inv_sq_ob_diffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type un_dr_sl(un_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type un_diff(un_diffSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ob_dr_sl(ob_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sq_ob_diff(inv_sq_ob_diffSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_filter_vcov_are(un_dr_sl, un_diff, ob_dr_sl, inv_sq_ob_diff));
    return rcpp_result_gen;
END_RCPP
}
// calc_filter_mean
arma::mat calc_filter_mean(arma::cube un_dr_sl, arma::mat un_dr_in, arma::cube ob_dr_sl, arma::mat ob_dr_in, arma::cube inv_sq_ob_diff, arma::cube vcov, arma::vec init, double delta, arma::mat deltaY, int subsump_rate);
RcppExport SEXP _yuima_calc_filter_mean(SEXP un_dr_slSEXP, SEXP un_dr_inSEXP, SEXP ob_dr_slSEXP, SEXP ob_dr_inSEXP, SEXP inv_sq_ob_diffSEXP, SEXP vcovSEXP, SEXP initSEXP, SEXP deltaSEXP, SEXP deltaYSEXP, SEXP subsump_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type un_dr_sl(un_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type un_dr_in(un_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ob_dr_sl(ob_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ob_dr_in(ob_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type inv_sq_ob_diff(inv_sq_ob_diffSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type vcov(vcovSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type deltaY(deltaYSEXP);
    Rcpp::traits::input_parameter< int >::type subsump_rate(subsump_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_filter_mean(un_dr_sl, un_dr_in, ob_dr_sl, ob_dr_in, inv_sq_ob_diff, vcov, init, delta, deltaY, subsump_rate));
    return rcpp_result_gen;
END_RCPP
}
// calc_filter_mean_time_homogeneous_with_vcov_are
arma::mat calc_filter_mean_time_homogeneous_with_vcov_are(arma::mat un_dr_sl, arma::vec un_dr_in, arma::mat ob_dr_sl, arma::vec ob_dr_in, arma::mat inv_sq_ob_diff, arma::mat vcov, arma::vec init, double delta, arma::mat deltaY);
RcppExport SEXP _yuima_calc_filter_mean_time_homogeneous_with_vcov_are(SEXP un_dr_slSEXP, SEXP un_dr_inSEXP, SEXP ob_dr_slSEXP, SEXP ob_dr_inSEXP, SEXP inv_sq_ob_diffSEXP, SEXP vcovSEXP, SEXP initSEXP, SEXP deltaSEXP, SEXP deltaYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type un_dr_sl(un_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type un_dr_in(un_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ob_dr_sl(ob_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ob_dr_in(ob_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sq_ob_diff(inv_sq_ob_diffSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vcov(vcovSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type deltaY(deltaYSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_filter_mean_time_homogeneous_with_vcov_are(un_dr_sl, un_dr_in, ob_dr_sl, ob_dr_in, inv_sq_ob_diff, vcov, init, delta, deltaY));
    return rcpp_result_gen;
END_RCPP
}
// calc_filter_mean_time_homogeneous
arma::mat calc_filter_mean_time_homogeneous(arma::mat un_dr_sl, arma::vec un_dr_in, arma::mat ob_dr_sl, arma::vec ob_dr_in, arma::mat inv_sq_ob_diff, arma::cube vcov, arma::vec init, double delta, arma::mat deltaY, int subsump_rate);
RcppExport SEXP _yuima_calc_filter_mean_time_homogeneous(SEXP un_dr_slSEXP, SEXP un_dr_inSEXP, SEXP ob_dr_slSEXP, SEXP ob_dr_inSEXP, SEXP inv_sq_ob_diffSEXP, SEXP vcovSEXP, SEXP initSEXP, SEXP deltaSEXP, SEXP deltaYSEXP, SEXP subsump_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type un_dr_sl(un_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type un_dr_in(un_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ob_dr_sl(ob_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ob_dr_in(ob_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sq_ob_diff(inv_sq_ob_diffSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type vcov(vcovSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type deltaY(deltaYSEXP);
    Rcpp::traits::input_parameter< int >::type subsump_rate(subsump_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_filter_mean_time_homogeneous(un_dr_sl, un_dr_in, ob_dr_sl, ob_dr_in, inv_sq_ob_diff, vcov, init, delta, deltaY, subsump_rate));
    return rcpp_result_gen;
END_RCPP
}
// calc_filter_mean_explicit
arma::mat calc_filter_mean_explicit(arma::mat un_dr_sl, arma::vec un_dr_in, arma::mat ob_dr_sl, arma::vec ob_dr_in, arma::mat inv_sq_ob_diff, arma::mat vcov, arma::vec init, double delta, arma::mat deltaY);
RcppExport SEXP _yuima_calc_filter_mean_explicit(SEXP un_dr_slSEXP, SEXP un_dr_inSEXP, SEXP ob_dr_slSEXP, SEXP ob_dr_inSEXP, SEXP inv_sq_ob_diffSEXP, SEXP vcovSEXP, SEXP initSEXP, SEXP deltaSEXP, SEXP deltaYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type un_dr_sl(un_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type un_dr_in(un_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ob_dr_sl(ob_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ob_dr_in(ob_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sq_ob_diff(inv_sq_ob_diffSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vcov(vcovSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type init(initSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type deltaY(deltaYSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_filter_mean_explicit(un_dr_sl, un_dr_in, ob_dr_sl, ob_dr_in, inv_sq_ob_diff, vcov, init, delta, deltaY));
    return rcpp_result_gen;
END_RCPP
}
// calc_kalman_bucy_filter_cpp
Rcpp::List calc_kalman_bucy_filter_cpp(arma::cube un_dr_sl, arma::mat un_dr_in, arma::cube un_diff, arma::cube ob_dr_sl, arma::mat ob_dr_in, arma::cube inv_sq_ob_diff, arma::mat vcov_init, arma::vec mean_init, double delta, arma::mat deltaY, bool use_are, bool is_explicit, bool is_time_homogeneous, int subsump_rate);
RcppExport SEXP _yuima_calc_kalman_bucy_filter_cpp(SEXP un_dr_slSEXP, SEXP un_dr_inSEXP, SEXP un_diffSEXP, SEXP ob_dr_slSEXP, SEXP ob_dr_inSEXP, SEXP inv_sq_ob_diffSEXP, SEXP vcov_initSEXP, SEXP mean_initSEXP, SEXP deltaSEXP, SEXP deltaYSEXP, SEXP use_areSEXP, SEXP is_explicitSEXP, SEXP is_time_homogeneousSEXP, SEXP subsump_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type un_dr_sl(un_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type un_dr_in(un_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type un_diff(un_diffSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ob_dr_sl(ob_dr_slSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ob_dr_in(ob_dr_inSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type inv_sq_ob_diff(inv_sq_ob_diffSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vcov_init(vcov_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mean_init(mean_initSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type deltaY(deltaYSEXP);
    Rcpp::traits::input_parameter< bool >::type use_are(use_areSEXP);
    Rcpp::traits::input_parameter< bool >::type is_explicit(is_explicitSEXP);
    Rcpp::traits::input_parameter< bool >::type is_time_homogeneous(is_time_homogeneousSEXP);
    Rcpp::traits::input_parameter< int >::type subsump_rate(subsump_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_kalman_bucy_filter_cpp(un_dr_sl, un_dr_in, un_diff, ob_dr_sl, ob_dr_in, inv_sq_ob_diff, vcov_init, mean_init, delta, deltaY, use_are, is_explicit, is_time_homogeneous, subsump_rate));
    return rcpp_result_gen;
END_RCPP
}
// W1
double W1(NumericMatrix crossdx, NumericMatrix b, NumericMatrix A, double h);
RcppExport SEXP _yuima_W1(SEXP crossdxSEXP, SEXP bSEXP, SEXP ASEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type crossdx(crossdxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(W1(crossdx, b, A, h));
    return rcpp_result_gen;
END_RCPP
}
// W2
double W2(NumericMatrix dx, NumericMatrix b, double h);
RcppExport SEXP _yuima_W2(SEXP dxSEXP, SEXP bSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(W2(dx, b, h));
    return rcpp_result_gen;
END_RCPP
}
// Irregular_PseudoLoglik_COG
double Irregular_PseudoLoglik_COG(int lengthObs, arma::mat B, arma::mat Btilde, arma::mat InvBtilde, double a0, double bq, double a1, double V, double PseudologLik, arma::rowvec ta, arma::colvec state, arma::colvec stateMean, arma::colvec e, arma::vec DeltaG2, arma::vec Deltat);
RcppExport SEXP _yuima_Irregular_PseudoLoglik_COG(SEXP lengthObsSEXP, SEXP BSEXP, SEXP BtildeSEXP, SEXP InvBtildeSEXP, SEXP a0SEXP, SEXP bqSEXP, SEXP a1SEXP, SEXP VSEXP, SEXP PseudologLikSEXP, SEXP taSEXP, SEXP stateSEXP, SEXP stateMeanSEXP, SEXP eSEXP, SEXP DeltaG2SEXP, SEXP DeltatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type lengthObs(lengthObsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Btilde(BtildeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type InvBtilde(InvBtildeSEXP);
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type bq(bqSEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type PseudologLik(PseudologLikSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type ta(taSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type state(stateSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type stateMean(stateMeanSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type e(eSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type DeltaG2(DeltaG2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Deltat(DeltatSEXP);
    rcpp_result_gen = Rcpp::wrap(Irregular_PseudoLoglik_COG(lengthObs, B, Btilde, InvBtilde, a0, bq, a1, V, PseudologLik, ta, state, stateMean, e, DeltaG2, Deltat));
    return rcpp_result_gen;
END_RCPP
}
// minusloglcpp_linear_state_space_theta1
double minusloglcpp_linear_state_space_theta1(double logdet_sq_observed_diffusion, arma::mat inv_sq_observed_diffusion, arma::mat dx, double h, int drop_terms);
RcppExport SEXP _yuima_minusloglcpp_linear_state_space_theta1(SEXP logdet_sq_observed_diffusionSEXP, SEXP inv_sq_observed_diffusionSEXP, SEXP dxSEXP, SEXP hSEXP, SEXP drop_termsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type logdet_sq_observed_diffusion(logdet_sq_observed_diffusionSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sq_observed_diffusion(inv_sq_observed_diffusionSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type drop_terms(drop_termsSEXP);
    rcpp_result_gen = Rcpp::wrap(minusloglcpp_linear_state_space_theta1(logdet_sq_observed_diffusion, inv_sq_observed_diffusion, dx, h, drop_terms));
    return rcpp_result_gen;
END_RCPP
}
// minusloglcpp_linear_state_space_theta2
double minusloglcpp_linear_state_space_theta2(arma::mat observed_drift, arma::cube inv_sq_observed_diffusion, arma::mat dx, double h, int drop_terms);
RcppExport SEXP _yuima_minusloglcpp_linear_state_space_theta2(SEXP observed_driftSEXP, SEXP inv_sq_observed_diffusionSEXP, SEXP dxSEXP, SEXP hSEXP, SEXP drop_termsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type observed_drift(observed_driftSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type inv_sq_observed_diffusion(inv_sq_observed_diffusionSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type drop_terms(drop_termsSEXP);
    rcpp_result_gen = Rcpp::wrap(minusloglcpp_linear_state_space_theta2(observed_drift, inv_sq_observed_diffusion, dx, h, drop_terms));
    return rcpp_result_gen;
END_RCPP
}
// calc_inverce_square
arma::cube calc_inverce_square(arma::cube cube);
RcppExport SEXP _yuima_calc_inverce_square(SEXP cubeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type cube(cubeSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_inverce_square(cube));
    return rcpp_result_gen;
END_RCPP
}
// driftTermCpp
NumericMatrix driftTermCpp(ExpressionVector drift, CharacterVector modelstate, arma::mat data, Environment env);
RcppExport SEXP _yuima_driftTermCpp(SEXP driftSEXP, SEXP modelstateSEXP, SEXP dataSEXP, SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ExpressionVector >::type drift(driftSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type modelstate(modelstateSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Environment >::type env(envSEXP);
    rcpp_result_gen = Rcpp::wrap(driftTermCpp(drift, modelstate, data, env));
    return rcpp_result_gen;
END_RCPP
}
// linearDriftTermCpp
arma::mat linearDriftTermCpp(arma::mat slope, arma::vec intercept, arma::mat data);
RcppExport SEXP _yuima_linearDriftTermCpp(SEXP slopeSEXP, SEXP interceptSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type slope(slopeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(linearDriftTermCpp(slope, intercept, data));
    return rcpp_result_gen;
END_RCPP
}
// diffusionTermCpp
NumericVector diffusionTermCpp(List diffusion, CharacterVector modelstate, arma::mat data, Environment env);
RcppExport SEXP _yuima_diffusionTermCpp(SEXP diffusionSEXP, SEXP modelstateSEXP, SEXP dataSEXP, SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type diffusion(diffusionSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type modelstate(modelstateSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Environment >::type env(envSEXP);
    rcpp_result_gen = Rcpp::wrap(diffusionTermCpp(diffusion, modelstate, data, env));
    return rcpp_result_gen;
END_RCPP
}
// measureTermCpp
NumericVector measureTermCpp(List measure, CharacterVector modelstate, arma::mat data, Environment env);
RcppExport SEXP _yuima_measureTermCpp(SEXP measureSEXP, SEXP modelstateSEXP, SEXP dataSEXP, SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type measure(measureSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type modelstate(modelstateSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Environment >::type env(envSEXP);
    rcpp_result_gen = Rcpp::wrap(measureTermCpp(measure, modelstate, data, env));
    return rcpp_result_gen;
END_RCPP
}
// detcpp
double detcpp(NumericMatrix A);
RcppExport SEXP _yuima_detcpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(detcpp(A));
    return rcpp_result_gen;
END_RCPP
}
// Smake
NumericMatrix Smake(NumericVector b, int d);
RcppExport SEXP _yuima_Smake(SEXP bSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(Smake(b, d));
    return rcpp_result_gen;
END_RCPP
}
// solvecpp
NumericMatrix solvecpp(NumericMatrix A);
RcppExport SEXP _yuima_solvecpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(solvecpp(A));
    return rcpp_result_gen;
END_RCPP
}
// sub_f
double sub_f(NumericMatrix S, NumericVector b);
RcppExport SEXP _yuima_sub_f(SEXP SSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(sub_f(S, b));
    return rcpp_result_gen;
END_RCPP
}
// likndim
double likndim(NumericMatrix dx, NumericMatrix b, NumericMatrix A, double h);
RcppExport SEXP _yuima_likndim(SEXP dxSEXP, SEXP bSEXP, SEXP ASEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(likndim(dx, b, A, h));
    return rcpp_result_gen;
END_RCPP
}
// residualCpp
NumericVector residualCpp(NumericVector dx, NumericVector a, NumericVector b, double w, double h);
RcppExport SEXP _yuima_residualCpp(SEXP dxSEXP, SEXP aSEXP, SEXP bSEXP, SEXP wSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(residualCpp(dx, a, b, w, h));
    return rcpp_result_gen;
END_RCPP
}
