#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bibsynchro(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ctsubsampling(void *, void *, void *, void *, void *, void *, void *);
extern void HayashiYoshida(void *, void *, void *, void *, void *, void *, void *);
extern void hyavar(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hycrossavar(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void HYcrosscorr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void HYcrosscorr2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void HYcrosscov(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void HYcrosscov2(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void krprod(void *, void *, void *, void *);
extern void msrc(void *, void *, void *, void *, void *, void *, void *);
extern void pHayashiYoshida(void *, void *, void *, void *, void *, void *, void *, void *);
extern void refreshsampling(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void refreshsamplingphy(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rGIG(void *, void *, void *, void *, void *);
extern void rpts(void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _yuima_cpp_collapse(SEXP, SEXP);
extern SEXP _yuima_cpp_E(SEXP);
extern SEXP _yuima_cpp_ito(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_cpp_ito_outer(SEXP, SEXP);
extern SEXP _yuima_cpp_ito_product(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_cpp_label(SEXP);
extern SEXP _yuima_cpp_outer(SEXP, SEXP);
extern SEXP _yuima_cpp_paste(SEXP, SEXP, SEXP);
extern SEXP _yuima_cpp_split(SEXP, SEXP);
extern SEXP _yuima_detcpp(SEXP);
extern SEXP _yuima_evalKernelCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_evalKernelCpp2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_Irregular_PseudoLoglik_COG(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_is_zero(SEXP);
extern SEXP _yuima_likndim(SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_makeprop(SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_residualCpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_Smake(SEXP, SEXP);
extern SEXP _yuima_solvecpp(SEXP);
extern SEXP _yuima_sqnorm(SEXP);
extern SEXP _yuima_sub_f(SEXP, SEXP);
extern SEXP _yuima_W1(SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_W2(SEXP, SEXP, SEXP);
extern SEXP Cycle_Carma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP euler(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pseudoLoglik_COGARCH1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_driftTermCpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_diffusionTermCpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_measureTermCpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _yuima_minusloglcpp_linear_state_space_theta1(SEXP,SEXP,SEXP,SEXP);
extern SEXP _yuima_calc_inverse_square(SEXP);
extern SEXP _yuima_calc_kalman_bucy_filter_cpp(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);


static const R_CMethodDef CEntries[] = {
    {"bibsynchro",         (DL_FUNC) &bibsynchro,          9},
    {"ctsubsampling",      (DL_FUNC) &ctsubsampling,       7},
    {"HayashiYoshida",     (DL_FUNC) &HayashiYoshida,      7},
    {"hyavar",             (DL_FUNC) &hyavar,              9},
    {"hycrossavar",        (DL_FUNC) &hycrossavar,        19},
    {"HYcrosscorr",        (DL_FUNC) &HYcrosscorr,        12},
    {"HYcrosscorr2",       (DL_FUNC) &HYcrosscorr2,       11},
    {"HYcrosscov",         (DL_FUNC) &HYcrosscov,         10},
    {"HYcrosscov2",        (DL_FUNC) &HYcrosscov2,         9},
    {"krprod",             (DL_FUNC) &krprod,              4},
    {"msrc",               (DL_FUNC) &msrc,                7},
    {"pHayashiYoshida",    (DL_FUNC) &pHayashiYoshida,     8},
    {"refreshsampling",    (DL_FUNC) &refreshsampling,     9},
    {"refreshsamplingphy", (DL_FUNC) &refreshsamplingphy, 10},
    {"rGIG",               (DL_FUNC) &rGIG,                5},
    {"rpts",               (DL_FUNC) &rpts,                5},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_yuima_cpp_collapse",                    (DL_FUNC) &_yuima_cpp_collapse,                           2},
    {"_yuima_cpp_E",                           (DL_FUNC) &_yuima_cpp_E,                                  1},
    {"_yuima_cpp_ito",                         (DL_FUNC) &_yuima_cpp_ito,                                5},
    {"_yuima_cpp_ito_outer",                   (DL_FUNC) &_yuima_cpp_ito_outer,                          2},
    {"_yuima_cpp_ito_product",                 (DL_FUNC) &_yuima_cpp_ito_product,                        8},
    {"_yuima_cpp_label",                       (DL_FUNC) &_yuima_cpp_label,                              1},
    {"_yuima_cpp_outer",                       (DL_FUNC) &_yuima_cpp_outer,                              2},
    {"_yuima_cpp_paste",                       (DL_FUNC) &_yuima_cpp_paste,                              3},
    {"_yuima_cpp_split",                       (DL_FUNC) &_yuima_cpp_split,                              2},
    {"_yuima_detcpp",                          (DL_FUNC) &_yuima_detcpp,                                 1},
    {"_yuima_evalKernelCpp",                   (DL_FUNC) &_yuima_evalKernelCpp,                         10},
    {"_yuima_evalKernelCpp2",                  (DL_FUNC) &_yuima_evalKernelCpp2,                        13},
    {"_yuima_Irregular_PseudoLoglik_COG",      (DL_FUNC) &_yuima_Irregular_PseudoLoglik_COG,            15},
    {"_yuima_is_zero",                         (DL_FUNC) &_yuima_is_zero,                                1},
    {"_yuima_likndim",                         (DL_FUNC) &_yuima_likndim,                                4},
    {"_yuima_makeprop",                        (DL_FUNC) &_yuima_makeprop,                               4},
    {"_yuima_residualCpp",                     (DL_FUNC) &_yuima_residualCpp,                            5},
    {"_yuima_Smake",                           (DL_FUNC) &_yuima_Smake,                                  2},
    {"_yuima_solvecpp",                        (DL_FUNC) &_yuima_solvecpp,                               1},
    {"_yuima_sqnorm",                          (DL_FUNC) &_yuima_sqnorm,                                 1},
    {"_yuima_sub_f",                           (DL_FUNC) &_yuima_sub_f,                                  2},
    {"_yuima_W1",                              (DL_FUNC) &_yuima_W1,                                     4},
    {"_yuima_W2",                              (DL_FUNC) &_yuima_W2,                                     3},
    {"Cycle_Carma",                            (DL_FUNC) &Cycle_Carma,                                  12},
    {"euler",                                  (DL_FUNC) &euler,                                        11},
    {"pseudoLoglik_COGARCH1",                  (DL_FUNC) &pseudoLoglik_COGARCH1,                        14},
    {"_yuima_driftTermCpp",                           (DL_FUNC) &_yuima_driftTermCpp,                                  4},
    {"_yuima_diffusionTermCpp",                       (DL_FUNC) &_yuima_diffusionTermCpp,                              4},
    {"_yuima_measureTermCpp",                         (DL_FUNC) &_yuima_measureTermCpp,                                4},
    {"_yuima_minusloglcpp_linear_state_space_theta1", (DL_FUNC) &_yuima_minusloglcpp_linear_state_space_theta1,        4},
    {"_yuima_calc_inverse_square",                    (DL_FUNC) &_yuima_calc_inverse_square,                           1},
    {"_yuima_calc_kalman_bucy_filter_cpp",            (DL_FUNC) &_yuima_calc_kalman_bucy_filter_cpp,                  16},
    {NULL, NULL, 0}
};

void R_init_yuima(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
