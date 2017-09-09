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
extern SEXP Cycle_Carma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP euler(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pseudoLoglik_COGARCH1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP yuima_detcpp(SEXP);
extern SEXP yuima_Irregular_PseudoLoglik_COG(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP yuima_likndim(SEXP, SEXP, SEXP, SEXP);
extern SEXP yuima_makeprop(SEXP, SEXP, SEXP, SEXP);
extern SEXP yuima_Smake(SEXP, SEXP);
extern SEXP yuima_solvecpp(SEXP);
extern SEXP yuima_sqnorm(SEXP);
extern SEXP yuima_sub_f(SEXP, SEXP);
extern SEXP yuima_W1(SEXP, SEXP, SEXP, SEXP);
extern SEXP yuima_W2(SEXP, SEXP, SEXP);

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
    {"Cycle_Carma",                      (DL_FUNC) &Cycle_Carma,                      12},
    {"euler",                            (DL_FUNC) &euler,                            11},
    {"pseudoLoglik_COGARCH1",            (DL_FUNC) &pseudoLoglik_COGARCH1,            14},
    {"yuima_detcpp",                     (DL_FUNC) &yuima_detcpp,                      1},
    {"yuima_Irregular_PseudoLoglik_COG", (DL_FUNC) &yuima_Irregular_PseudoLoglik_COG, 15},
    {"yuima_likndim",                    (DL_FUNC) &yuima_likndim,                     4},
    {"yuima_makeprop",                   (DL_FUNC) &yuima_makeprop,                    4},
    {"yuima_Smake",                      (DL_FUNC) &yuima_Smake,                       2},
    {"yuima_solvecpp",                   (DL_FUNC) &yuima_solvecpp,                    1},
    {"yuima_sqnorm",                     (DL_FUNC) &yuima_sqnorm,                      1},
    {"yuima_sub_f",                      (DL_FUNC) &yuima_sub_f,                       2},
    {"yuima_W1",                         (DL_FUNC) &yuima_W1,                          4},
    {"yuima_W2",                         (DL_FUNC) &yuima_W2,                          3},
    {NULL, NULL, 0}
};

void R_init_yuima(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
