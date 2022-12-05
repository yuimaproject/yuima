//
//  euler.c
//  
//
//  Created by Yuta Koike on 2016/01/23.
//
//

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP euler(SEXP x0, SEXP t0, SEXP R, SEXP dt, SEXP dW, SEXP modeltime, SEXP modelstate,
           SEXP drift, SEXP diffusion, SEXP env, SEXP rho){
    
    int i, j, k, n, d, r;
    double *rdt, *rdW, *rX, *rx0, *b, *sigma;
    SEXP X, xpar, tpar, b0, sigma0;
    //SEXP X, xpar, tpar, xvar, tvar;
    
    PROTECT(x0 = AS_NUMERIC(x0));
    rx0 = REAL(x0);
    PROTECT(dt = AS_NUMERIC(dt));
    rdt = REAL(dt);
    PROTECT(dW = AS_NUMERIC(dW));
    rdW = REAL(dW);
    
    d = length(x0);
    n = length(dt);
    
    r = *INTEGER(R);
    
    /* PROTECT(X = allocVector(REALSXP, n + 1));
     rX = REAL(X);
     rX[0] = *REAL(x0);*/
    PROTECT(X = allocMatrix(REALSXP, d, n + 1));
    rX = REAL(X);
    
    for (j = 0; j < d; j++) {
        rX[j] = rx0[j];
    }
    
    PROTECT(tpar = allocVector(REALSXP, 1));
    //PROTECT(xpar = allocVector(REALSXP, 1));
    
    PROTECT(t0 = AS_NUMERIC(t0));
    REAL(tpar)[0] = REAL(t0)[0]; /* initial time */
    
    //PROTECT(b0 = allocVector(REALSXP, d));
    //PROTECT(sigma0 = allocVector(REALSXP, d*r));
    PROTECT(b0 = allocVector(VECSXP, 1));
    PROTECT(sigma0 = allocVector(VECSXP, 1));
    
    PROTECT(xpar = allocVector(REALSXP, 1));
     
    for (i = 0; i < n; i++) {
        
        /* assign the current variables */
        defineVar(installChar(STRING_ELT(modeltime, 0)), tpar, env); 
        // PROTECT(tvar = STRING_ELT(modeltime, 0));
        //defineVar(installChar(tvar), tpar, env); 
        
        for (j = 0; j < d; j++) {
            /* PROTECT(xpar = allocVector(REALSXP, 1));*/
            REAL(xpar)[0] = rX[j + i * d];
            /* defineVar(installChar(STRING_ELT(modelstate, j)), duplicate(xpar), env); */
            /* next two lines to fix possibile memory leaks proposed by Tomas Kalibera */
            SEXP ch = installChar(STRING_ELT(modelstate, j));
            defineVar(ch, duplicate(xpar), env);
            
            //PROTECT(xvar = STRING_ELT(modelstate, j));
            //defineVar(installChar(xvar), duplicate(xpar), env);
            //defineVar(installChar(STRING_ELT(modelstate, j)), xpar, env);
/*            UNPROTECT(1); */
        }
        
        /*defineVar(install("env"), env, rho);*/
        
        /* evaluate coefficients */
        /* PROTECT(b0 = eval(drift, rho));
        PROTECT(sigma0 = eval(diffusion, rho));
        AS_NUMERIC is added by YK (Dec 4, 2018) */
        /* PROTECT(b0 = AS_NUMERIC(eval(drift, rho)));
        PROTECT(sigma0 = AS_NUMERIC(eval(diffusion, rho)));
        */
        /* PROTECT(b0 = allocVector(REALSXP, d));  */
        //b0 = AS_NUMERIC(eval(drift, rho));
        //b0 = PROTECT(AS_NUMERIC(eval(drift, rho)));
        /* PROTECT(sigma0 = allocVector(REALSXP, d*r)); */
        //sigma0 = AS_NUMERIC(eval(diffusion, rho)); 
        //sigma0 = PROTECT(AS_NUMERIC(eval(diffusion, rho))); 
        SET_VECTOR_ELT(b0, 0, eval(drift, rho));
        SET_VECTOR_ELT(sigma0, 0, eval(diffusion, rho));
        
        //b = REAL(b0);
        //sigma = REAL(sigma0);
        b = REAL(VECTOR_ELT(b0, 0));
        sigma = REAL(VECTOR_ELT(sigma0, 0));
        
        for (j = 0; j < d; j++) {
            rX[j + (i + 1) * d] = rX[j + i * d] + b[j] * rdt[i];
            for (k = 0; k < r; k++) {
                rX[j + (i + 1) * d] += sigma[k + j * r] * rdW[k + i * r];
            }
        }
        
        /*defineVar(install("t"), tpar, rho);*/
        /*rX[i + 1] = rX[i] + *REAL(eval(drift, rho)) * rdt[i] + *REAL(eval(diffusion, rho)) * rdW[i];*/
        /*rX[i + 1] = rX[i] + *REAL(eval(drift, rho)) * REAL(dt)[i] + *REAL(eval(diffusion, rho)) * REAL(dW)[i];*/
        
        REAL(tpar)[0] += rdt[i];
       
       
    }
    UNPROTECT(1); /* xpar */
    UNPROTECT(2); /* b0, sigma0 */ 
    UNPROTECT(6);
    return(X);
}
