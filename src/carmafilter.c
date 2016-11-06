/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Code of this package: Copyright (C) 2006 S. M. Iacus
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * Exports
 *	sde_sim_xxx(...)
 *
 * to be called as  .C(.)  in ../R/sde.sim.xxx.R
 * where xxx is one among "euler", "milstein", "milstei2", "KPS"
 */



#include <R.h>
#include <Rmath.h>
#include <R_ext/Boolean.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Complex.h>

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)





SEXP carma_tmp(SEXP V, SEXP P, SEXP A);

SEXP Cycle_Carma(SEXP Y, SEXP StateVar, SEXP ExpA, SEXP Times_Obs, SEXP P,
                SEXP Qmatr, SEXP SigMatr, SEXP Zc, SEXP Result,SEXP Kgain,
                SEXP dum_zc, SEXP Mat22int);



SEXP carma_tmp(SEXP V, SEXP P, SEXP A){

    int p;
    int i, j, h;
    double *rV, *rA, *rB, *rC, *rSigma;
    SEXP  B, C, Sigma;

    if(!isInteger(P)) error("`P' must be integer");
    if(!isNumeric(V)) error("`V' must be numeric");
    if(!isNumeric(A)) error("`A' must be numeric");


    PROTECT(V = AS_NUMERIC(V));
    rV = REAL(V);
    PROTECT(A = AS_NUMERIC(A));
    rA = REAL(A);

    p = *INTEGER(P);


    PROTECT(B = allocMatrix(REALSXP, p, p));
    rB = REAL(B);

    PROTECT(C  = allocMatrix(REALSXP, p, p));
    rC = REAL(C);


    PROTECT(Sigma  = allocMatrix(REALSXP, p, p));
    rSigma = REAL(Sigma);

    /* B = A %*% V */
    for(i=0; i<p; i++){
        for(j=0; j<p; j++){
            rB[i+j*p] = 0;
            for(h=0; h<p; h++){
                rB[i+j*p] = rB[i+j*p] + rA[i+h*p] * rV[h+j*p];
            }
        }
   /*}
     C = B %*% A^T
    for(i=0; i<p; i++){*/
        for(j=0; j<p; j++){
            rC[i+j*p] = 0;
            for(h=0; h<p; h++){
                rC[i+j*p] = rC[i+j*p] + rB[i+h*p] * rA[j+h*p];
            }
          rSigma[i+j*p] = rV[i+j*p] - rC[i+j*p];
        }
    }
    /*for(i=0; i<p*p; i++){

        }*/
    UNPROTECT(5);
    return Sigma;
}

SEXP Cycle_Carma(SEXP Y, SEXP StateVar, SEXP ExpA, SEXP Times_Obs, SEXP P,
                SEXP Qmatr, SEXP SigMatr, SEXP Zc, SEXP Result, SEXP Kgain,
                SEXP dum_zc, SEXP Mat22int){

                /* Declare Integer Variable */

                int times_obs, p;
                int i, j, h,  t;


                /* Declare pointer */

                double *rY, *rStateVar, *rExpA, *rQmatr, *rSigMatr, *rZc;
                double *rKgain, *rdum_zc,  *rMat22int,  *rResult;
                double Uobs=0;
                double dummy=0;
                double dummy1=0;
                /*FILE *fd;*/
                /* Declare SEXP */



                /* Check the type of variables*/

              //  if(!isInteger(P)) error("`P' must be integer");
              //  if(!isInteger(Times_Obs)) error("`Times_Obs' must be integer");
              //  if(!isNumeric(Y)) error("`Y' must be numeric");
              //  if(!isNumeric(StateVar)) error("`StateVar' must be numeric");
              //  if(!isNumeric(ExpA)) error("`ExpA' must be numeric");
              //  if(!isNumeric(Qmatr)) error("`Qmatr' must be numeric");
              //  if(!isNumeric(SigMatr)) error("`SigMatr' must be numeric");
              //  if(!isNumeric(Zc)) error("`Zc' must be numeric");


                /* Protect Objects */

                PROTECT(Y = AS_NUMERIC(Y));
                rY = REAL(Y);

                PROTECT(StateVar = AS_NUMERIC(StateVar));
                rStateVar = REAL(StateVar);

                PROTECT(ExpA = AS_NUMERIC(ExpA));
                rExpA = REAL(ExpA);

                PROTECT(Qmatr = AS_NUMERIC(Qmatr));
                rQmatr = REAL(Qmatr);

                PROTECT(SigMatr = AS_NUMERIC(SigMatr));
                rSigMatr = REAL(SigMatr);

                PROTECT(Zc = AS_NUMERIC(Zc));
                rZc = REAL(Zc);

                PROTECT(Result = AS_NUMERIC(Result));
                rResult = REAL(Result);

                PROTECT(Kgain = AS_NUMERIC(Kgain));
                rKgain = REAL(Kgain);

                PROTECT(dum_zc = AS_NUMERIC(dum_zc));
                rdum_zc = REAL(dum_zc);

                PROTECT(Mat22int = AS_NUMERIC(Mat22int));
                rMat22int = REAL(Mat22int);


                /*PROTECT(Logstar = AS_NUMERIC(Logstar));
                rLoglstar = REAL(Logstar);*/


                /* Declare dimensions for the state variables and observations */
                times_obs = *INTEGER(Times_Obs);
                p = *INTEGER(P);


               /* Main Code
                 Dimension of Inputs:
                Y = Vector p dimension;
                StateVar = matrix p x 1;
                ExpA = matrix p x p;
                Times_Obs = integer;
                P = integer;
                Qmatr = matrix p x p;
                SigMatr = matrix p x p;
                Zc =  matrix 1 x p;
                Logstar = scalar.

                for(t in 1:rtime_obs){*/

                rResult[0]=0;

                 for(t=0; t<times_obs; t++){
              /*t=0;*/
                      /* prediction */
                      for(i=0; i<p; i++){
                        dummy = 0;
                        for(j=0; j<p; j++){
                          /*    statevar <- expA %*% statevar: px1=pxp px1 */
                          dummy += rExpA[i+j*p]*rStateVar[j];
                        }
                          rStateVar[i] = dummy;
                      /*SigMatr <- expA %*% SigMatr %*% expAT + Qmatr: pxp = pxp pxp pxp
                         First We compute rMat22int <- expA %*% SigMatr : pxp = pxp pxp
                      for(i=0; i<p; i++){*/
                        for(j=0; j<p; j++){
                          rMat22int[i+j*p]=0;
                          for(h=0; h<p; h++){
                              rMat22int[i+j*p]=rMat22int[i+j*p]+rExpA[i+h*p]*
                              rSigMatr[h+j*p];
                          }
                        }
                       }
                       /* Second We compute rMat22est <- rMat22int %*% t(expA)
                         + Qmatr: pxp = pxp pxp + pxp */
                      Uobs = 0;
                      for(i=0; i<p; i++){
                        for(j=0; j<p; j++){
                          rSigMatr[i+j*p]=0;
                          for(h=0; h<p; h++){
                            rSigMatr[i+j*p] += rMat22int[i+h*p]*rExpA[j+h*p];

                          }
                          rSigMatr[i+j*p] = rSigMatr[i+j*p]+rQmatr[i+j*p];
                        }
                      /*forecast
                      Uobs <- y[t] - zc %*% statevar # 1 - 1Xp px*/
                      Uobs += rZc[i]*rStateVar[i];
                      }
                      Uobs = rY[t]-Uobs;
                      /*   dum.zc <- zc %*% SigMatr  # 1xp pxp */
                      rResult[1] =0;
                      for(i=0; i<p; i++){
                        dummy = 0;
                        for(j=0; j<p; j++){
                          dummy += rSigMatr[i+j*p]*rZc[j];
                        }
                        rdum_zc[i]=dummy;
                      /*     Sd_2 <- dum.zc %*% zcT  # 1xp px1 */
                         rResult[1] += rdum_zc[i]*rZc[i];
                      }
                      /* correction */
                      /*   Kgain <- SigMatr %*% zcT %*% 1/Sd_2 # pxp px1*1 */
                      for(i=0; i<p; i++){
                        dummy1=0;
                        for(j=0; j<p; j++){
                          dummy1 += rSigMatr[i+j*p]*rZc[j];
                        }
                        rKgain[i]=dummy1/rResult[1];
                      /*     statevar <- statevar+Kgain %*% Uobs # px1+px1 */
                         rStateVar[i] +=  rKgain[i]*Uobs;
                      /*SigMatr <- SigMatr - Kgain %*% dum.zc # pxp-px1 1xp*/
                        for(j=0; j<p; j++){
                            rSigMatr[i+j*p] += -rKgain[i]*rdum_zc[j];
                        }
                       }
                       /*term_int = -0.5 * (log(Sd_2)+ Uobs * Uobs * 1/Sd_2)  every entries are scalars*/
                       /*fd=fopen("dueinteri.txt", "w+");*/
                      /*if(rResult[1]>0){*/
                      rResult[0] += -0.5 * (log(rResult[1])+ Uobs * Uobs /rResult[1]);
                      /*}else{
                        rResult[0] += -10000000000;
                      }*/
                      /*printf("\n res %.5f", rResult[0]);*/
                      /* manual debug */
                   }


                 /*fclose(fd);*/

                   UNPROTECT(10);
                   return Result;

}



/*
SEXP sde_sim_euler(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                   SEXP d,888 SEXP s, SEXP sx,
				   SEXP eta, SEXP alpha, SEXP corr, SEXP rho)
{
  SEXP X, Y1, Y2;
  double  *rY1, *rY2, T1, T2, *rX;
  double DELTA, ETA, ALPHA;
  double sdt, Z, tmp, d1, d2, *rx0;
  Rboolean CORR;
  int i, n, j, m;

  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");
  if(!isInteger(M)) error("`M' must be integer");

  if(!isFunction(d)) error("`d' must be a function");
  if(!isFunction(s)) error("`s' must be a function");
  if(!isFunction(sx)) error("`sx' must be a function");

  if(!isNumeric(eta)) error("`eta' must be numeric");
  if(!isNumeric(alpha)) error("`alpha' must be numeric");
  if(!isLogical(corr)) error("`corr' must be logical");

  if(!isEnvironment(rho)) error("`rho' must be an environment");

  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(eta = AS_NUMERIC(eta));
  PROTECT(alpha = AS_NUMERIC(alpha));
  PROTECT(corr = AS_LOGICAL(corr));

  n = *INTEGER(N);
  m = *INTEGER(M);

  PROTECT(Y1 = NEW_NUMERIC(m));
  PROTECT(Y2 = NEW_NUMERIC(m));
  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1));
  rX = REAL(X);
  rY1 = REAL(Y1);
  rY2 = REAL(Y2);
  rx0 = REAL(x0);
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j];
  T1 = *REAL(t0);
  DELTA = *REAL(delta);
  ETA = *REAL(eta);
  ALPHA = *REAL(alpha);
  CORR = *LOGICAL(corr);

  sdt = sqrt(DELTA);

  for(j=0; j<m; j++)
   rY1[j] = rX[j*(n+1)];

  GetRNGstate();
  if(CORR==TRUE){
   for(i=1; i< n+1; i++){
    T2 = T1 + DELTA;
    for(j=0; j<m; j++){
 	 Z = rnorm(0,sdt);
	 tmp = rX[i-1 +j*(n+1)];
	 rY2[j] = tmp + feval(T1,tmp,d, rho)*DELTA + feval(T1,tmp,s, rho)*Z;
	 d1 = feval(T2,rY2[j],d,rho) - ETA*feval(T2,rY2[j], s, rho)*feval(T2,rY2[j],sx,rho);
	 d2 = feval(T2,tmp,d,rho) - ETA*feval(T2,tmp, s, rho)*feval(T2,tmp,sx,rho);
	 rX[i +j*(n+1)] = tmp + (ALPHA*d1 + (1.0-ALPHA)*d2)*DELTA +
               (ETA * feval(T2,rY2[j],s,rho) + (1.0-ETA)*feval(T1,rY1[j],s,rho))*Z;
 	 rY1[j] = rY2[j];
	}
	T1 = T2;
   }
  } else {
    for(i=1; i< n+1; i++){
	 T1 = T1 + DELTA;
     for(j=0; j<m; j++){
	  Z = rnorm(0, sdt);
	  tmp = rX[i + j*(n+1) - 1];
      rX[i + (n+1)*j] = tmp + feval(T1,tmp,d, rho)*DELTA + feval(T1,tmp,s, rho)*Z;
	  }
	}
  }

  PutRNGstate();

  UNPROTECT(9);
  return(X);
}


*/
