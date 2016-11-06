/* The function computes the pseudo-loglikelihood of a COGARCH(p,q) model
 * See Iacus et al. 2015 for details
 */
#include <Rmath.h>
#include <R_ext/Boolean.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Complex.h>

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)

/*SEXP myfun1(SEXP a, SEXP b);*/

SEXP pseudoLoglik_COGARCH1(SEXP a0, SEXP bq, SEXP a1, SEXP stateMean, SEXP Q,
                           SEXP DeltaG2, SEXP Deltat, SEXP DeltaB1,
                           SEXP a, SEXP e,
                           SEXP V, SEXP nObs,
                           SEXP dummyMatr, SEXP dummyeB1);

/*SEXP myfun1(SEXP a, SEXP b){
     double *ra, *rb, *rvai;
     double dummy1=0;
     int i, j;
     SEXP vai;

  PROTECT(a = AS_NUMERIC(a));
  ra = REAL(a);

  PROTECT(b = AS_NUMERIC(b));
  rb = REAL(b);

  PROTECT(vai = NEW_NUMERIC(1));
  rvai = REAL(vai);

  printf("\n  Q  %f.5",  dummy1);
  for(i=0; i < 2; i++){
    dummy1 = 0;
    for(j=0; j < 2; j++){
      dummy1 += ra[i+j*2]*rb[j];
      printf("\n  Q c: %d, %d %f.5", i, j, dummy1);
    }
    rvai[i]= dummy1;

  }
 UNPROTECT(3);
 return vai;
 }
*/



SEXP pseudoLoglik_COGARCH1(SEXP a0, SEXP bq, SEXP a1, SEXP stateMean, SEXP Q,
                           SEXP DeltaG2, SEXP Deltat, SEXP DeltaB1,
                           SEXP a, SEXP e,
                           SEXP V, SEXP nObs,
                           SEXP dummyMatr, SEXP dummyeB1){

                          /* Declare Integer Variable */
                          int numb_Obs, q, t, i, j;
                          double *ra0, *rbq, *ra1, *rstateMean, *rstate;
                          double *rDeltaG2, *rDeltat, *rDeltaB1;
                          double *ra, *re, *rV, rVarDeltaG=0;
                          double *rdummyMatr;
                          /*rPseudologLik=0,*/
                          double *rdummyeB1;
                          double dummy1=0;
                          double dummy2=0;
                          double *res;
                          double *rstatedum;
                          SEXP PseudoLoglik;
                          SEXP state;
                          SEXP statedum;




                          /* Protect Objects */

                          PROTECT(a0 = AS_NUMERIC(a0));
                          ra0 = REAL(a0);

                          /*1*/

                          PROTECT(PseudoLoglik = NEW_NUMERIC(1));
                            res=REAL(PseudoLoglik);



                            /*2*/

                          PROTECT(bq = AS_NUMERIC(bq));
                          rbq = REAL(bq);

                          /*3*/

                          PROTECT(a1 = AS_NUMERIC(a1));
                          ra1 = REAL(a1);

                          /*4*/

                          PROTECT(stateMean = AS_NUMERIC(stateMean));
                          rstateMean = REAL(stateMean);

                          /*5*/



                          PROTECT(DeltaG2 = AS_NUMERIC(DeltaG2));
                          rDeltaG2 = REAL(DeltaG2);

                          /*6*/

                          PROTECT(Deltat = AS_NUMERIC(Deltat));
                          rDeltat = REAL(Deltat);

                          /*7*/

                          PROTECT(DeltaB1 = AS_NUMERIC(DeltaB1));
                          rDeltaB1 = REAL(DeltaB1);

                          /*8*/

                          PROTECT(a = AS_NUMERIC(a));
                          ra = REAL(a);

                          /*9*/


                          PROTECT(e = AS_NUMERIC(e));
                          re = REAL(e);

                          /*10*/


                          PROTECT(V = AS_NUMERIC(V));
                          rV = REAL(V);




                          /*11*/


                          PROTECT(dummyMatr = AS_NUMERIC(dummyMatr));
                          rdummyMatr = REAL(dummyMatr);

                          /*12*/

                          PROTECT(dummyeB1 = AS_NUMERIC(dummyeB1));
                          rdummyeB1 = REAL(dummyeB1);

                          /*13*/


                          /* Declare dimensions for the state variables and observations */
                          numb_Obs = *INTEGER(nObs);
                          q = *INTEGER(Q);

                          PROTECT(state = NEW_NUMERIC(q));
                          rstate = REAL(state);

                          PROTECT(statedum = NEW_NUMERIC(q));
                          rstatedum = REAL(statedum);

                          for(i=0; i<q; i++){
                            rstate[i]=rstateMean[i];
                          }

                          /*printf("\n  Q c: %d, %d ", q, numb_Obs);
                           for(i=0; i<q; i++){
                            printf("\n  RSTATE: %.5f %.5f %.5f", rstate[i], rstateMean[i], rdummyMatr[i]);
                            printf("\n  RMAtr: %.5f %.5f" , rDeltaB1[i], rdummyeB1[i]);
                          }*/

                          *res=0;
                          /*fd = fopen("dueinteri.txt", "r");
                           printf("\n %p %p", &rstate, &rstateMean);*/
                          for(t=0; t<numb_Obs; t++){
                            /* VarDeltaG <- as.numeric(a0*Deltat*bq/(bq-a1)+dummyMatr%*%(state-stateMean)) */
                            dummy1 = 0;
                            for(i=0; i<q; i++){
                              dummy1 += rdummyMatr[i]*(rstate[i]-rstateMean[i]);
                              /*printf("\n %d  %.10f %.10f %.10f %.10f", i, dummy1, rdummyMatr[i], rstate[i], rstateMean[i]);*/
                            }

                            /*printf("\n  dummy1 c: %.10f", dummy1);*/

                            rVarDeltaG = ra0[0]*rDeltat[0]*rbq[0]/(rbq[0]-ra1[0])+dummy1;

                            /* state <- DeltatB1%*%state+DeltaG2[i]/V*dummyeB1%*%state+a0*DeltaG2[i]/V*e */

                            for(i=0; i<q; i++){
                              dummy1 = 0;
                              dummy2 = 0;
                              for(j=0; j<q; j++){
                                dummy1 += rDeltaB1[i+j*q]*rstate[j];
                                dummy2 += rdummyeB1[i+j*q]*rstate[j];

                              }
                              rstatedum[i]= dummy1+rDeltaG2[t]/rV[0]*dummy2+ra0[0]*rDeltaG2[t]/rV[0]*re[i];
                              /*rstate[i]= dummy1+dummy2;
                              printf("\n d1 %.10f d2 %.10f", dummy1, dummy2);*/
                            }

                            for(i=0; i<q; i++){
                              rstate[i]=rstatedum[i];
                            }
                            /* V <- as.numeric(a0+ta%*% state) */
                            rV[0] = 0;
                            for(i=0; i<q; i++){
                              rV[0] += ra[i]*rstate[i];
                            }
                            rV[0] = rV[0]+ra0[0];

                            /* PseudologLik<- -1/2*(DeltaG2[i]/VarDeltaG+log(VarDeltaG)+log(2*pi)) */

                            *res += 0.5*(-rDeltaG2[t]/rVarDeltaG-log(rVarDeltaG)-log(2.*3.14159265));
                          /*fscanf(fd, "%lf", &res);
                          printf("\n c: %.10f -  %.5f %.5f  -  %.5f", rVarDeltaG, rstate[0], rstate[1], rV[0]);*/

                          }

                          /*fclose(fd);*/


                          UNPROTECT(15);
                          return PseudoLoglik;

                          }


/*SEXP pseudoLoglik_COGARCH(SEXP a0, SEXP bq, SEXP a1, SEXP stateMean, SEXP Q,
                          SEXP state, SEXP DeltaG2, SEXP Deltat, SEXP DeltatB,
                          SEXP B, SEXP a,SEXP e, SEXP Btilde,
                          SEXP InvBtilde, SEXP V, SEXP I, SEXP VarDeltaG,
                          SEXP PseudologLik, SEXP nObs, SEXP expr, SEXP rho){

                           Declare Integer Variable



                         int numb_Obs, q, t;
                         SEXP ans, A;
                        double rA, rans;

                        if(!isNewList(DeltatB)) error("'DeltatB' must be a list");

                         if(!isEnvironment(rho)) error("`rho' must be an environment");
                          Protect Objects

                         numb_Obs = *INTEGER(nObs);

                         q = *INTEGER(Q);

                         PROTECT(B = allocMatrix(REALSXP, q, q));
                         PROTECT(DeltatB);

                         ans = PROTECT(allocVector(VECSXP, numb_Obs));

                         PROTECT(A = allocMatrix(REALSXP, q, q));
                         rA = *REAL(A);
                         rPseudologLik[0]=0;
                         for(t=0; t<numb_Obs; t++){

                            VarDeltaG <- as.numeric(a0*Deltat[i]*bq/(bq-a1)+t(a)%*
                            %expm(Btilde*Deltat[i]) %*%InvBtilde%*% (I-expm(-Btilde*Deltat[i]))
                            %*%(state-stateMean))



                           state <- (I+DeltaG2[i]/V*e%*%t(a))%*%expm(B*Deltat[i])
                                    %*%state+a0*DeltaG2[i]/V*e


                           V <- as.numeric(a0+t(a)%*% state)


                           PseudologLik<--1/2*(DeltaG2[i]/VarDeltaG+
                             log(VarDeltaG)+log(2*pi))

                              defineVar(install("x"), VECTOR_ELT(DeltatB,t), rho);
                              SET_VECTOR_ELT(ans, t, eval(expr, rho));

                         }


                         UNPROTECT(4);
                         return PseudologLik;
                         return ans;

 }
 */
