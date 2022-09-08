#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "R.h"



void rpts(int *x, double *alpha, double *a, double *b, double *rn);
void rGIG(int *x, double *lambda, double *chi, double *psi, double *rn);
double qdens(double m, double *lambda, double beta);
double rps(double *alpha, double *a, double *b);



void rGIG(int *x, double *lambda, double *chi, double *psi, double *rn)
{
    GetRNGstate();
  
    double beta = sqrt((*chi)*(*psi));
    int i = 0;
    /*Rejection method for non T-1/2-concave part*/
  // if(0.0<=*lambda<1.0 && 0.0<beta<=2.0/3.0*sqrt(1.0-*lambda)){  //warning: comparisons like ‘X<=Y<=Z’ do not have their mathematical meaning [-Wparentheses]
    if((*lambda>=0.0) && (*lambda<1.0) && (beta> 0.0) && (beta<=2.0/3.0*sqrt(1.0-*lambda)) ){
        double m = beta/((1.0-*lambda)+sqrt(pow(1.0-*lambda,2.0)+pow(beta,2.0)));
        double x0 = beta/(1.0-*lambda);
        double xstar = fmax(x0,2.0/beta);
        double k1 = qdens(m,lambda,beta);
        double A1 = k1*x0;
        double k2 = 0.0;
        double A2 = 0.0;
        if(x0<2.0/beta && *lambda>0.0){
            k2 = exp(-beta);
            A2 = k2*(pow(2.0/beta,*lambda)-pow(x0,*lambda))/(*lambda);
        }
        else if(x0<2.0/beta && *lambda == 0.0){
            k2 = exp(-beta);
            A2 = k2*log(2.0/pow(beta,2.0));
        }
        else{
            k2 = 0.0;
            A2 = 0.0;
        }
        double k3 = pow(xstar,*lambda-1.0);
        double A3 = 2.0*k3*exp(-xstar*beta/2.0)/beta;
        double A = A1+A2+A3;
        while(i<*x){
            double U = 1;
            double V = 0;
            double h = 1.0;
            double X = 10.0;
            while(U*h > qdens(X,lambda,beta)){
                U = unif_rand();
                V = A*unif_rand();
                if(V <= A1){
                    X = x0*V/A1;
                    h = k1;
                }
                else if(V <= A1+A2 && *lambda > 0.0){
                    V = V-A1;
                    X = pow((pow(x0,*lambda)+V*(*lambda)/k2),1.0/(*lambda));
                    h = k2*pow(X,*lambda-1.0);
                }
                else if(V <= A1+A2 && *lambda == 0.0){
                    V = V-A1;
                    X = beta*exp(V*exp(beta));
                    h = k2*pow(X,*lambda-1.0);
                }
                else{
                    V = V - (A1+A2);
                    X = -2.0/beta*log(exp(-xstar*beta/2.0)-V*beta/(2.0*k3));
                    h = k3*exp(-X*beta/2.0);
                }
            }
            rn[i] = X*sqrt(*chi/(*psi));
            i++;
        } /*Ratio-of-Uniforms without mode shift*/
//    }else if(0.0<=*lambda<=1.0 && fmin(1.0/2.0,2.0/3.0*sqrt(1.0-*lambda))<=beta<=1.0){ //warning: comparisons like ‘X<=Y<=Z’ do not have their mathematical meaning [-Wparentheses]
    } else if((*lambda>=0.0) && (*lambda<=1.0) && (beta>=fmin(1.0/2.0,2.0/3.0*sqrt(1.0-*lambda))) && (beta<=1.0)){
        double m = beta/((1.0-*lambda)+sqrt(pow(1.0-*lambda,2.0)+pow(beta,2.0)));
        double xplus = ((1.0+*lambda)+sqrt(pow(1.0+*lambda,2.0)+pow(beta,2.0)))/beta;
        double vplus = sqrt(qdens(m,lambda,beta));
        double uplus = xplus*sqrt(qdens(xplus,lambda,beta));
        while(i<*x){
            double V = 1.0;
            double X = 10.0;
            double U = 0.0;
            while(pow(V,2.0) > qdens(X,lambda,beta)){
                U = uplus*unif_rand();
                V = vplus*unif_rand();
                X = U/V;
            }
            rn[i] = X*sqrt(*chi/(*psi));
            i++;
        }
    }/*Ratio-of-Uniforms with mode shift (Dagpunar-Lehner)*/
    else{
        double m = (sqrt(pow(*lambda-1.0,2.0)+pow(beta,2.0))+(*lambda-1.0))/beta;
        double a = -2.0*(*lambda+1.0)/beta-m;
        double b = 2.0*(*lambda-1.0)/beta*m-1.0;
        double c = m;
        double p = b-pow(a,2.0)/3.0;
        double q = 2.0*pow(a,3.0)/27.0-a*b/3.0+c;
        double phi = acos(-q/2.0*sqrt(-27.0/pow(p,3.0)));
        double xminus = sqrt(-4.0/3.0*p)*cos(phi/3.0+4.0/3.0*M_PI)-a/3.0;
        double xplus = sqrt(-4.0/3.0*p)*cos(phi/3.0)-a/3.0;
        double vplus = sqrt(qdens(m,lambda,beta));
        double uminus = (xminus-m)*sqrt(qdens(xminus,lambda,beta));
        double uplus = (xplus-m)*sqrt(qdens(xplus,lambda,beta));
        while(i<*x){
            double V = 1.0;
            double X = 0.0;
            while(pow(V,2.0) > qdens(X,lambda,beta)){
                double U = uminus+(uplus-uminus)*unif_rand();
                V = vplus*unif_rand();
                X = U/V+m;
            }
            rn[i] = X*sqrt(*chi/(*psi));
            i++;
        }
    }
    PutRNGstate();
}



void rpts(int *x, double *alpha, double *a, double *b, double *rn)
{   
    GetRNGstate();
  
    int i=0;
    double y;
    double x1,y1,z1,uni;
    while(i<*x){
        uni=-M_PI/2.0+M_PI*unif_rand();
        x1=pow((*a)*tgamma(1.0-*alpha)*cos(M_PI*(*alpha)/(2.0))/(*alpha),1.0/(*alpha));
        y1=sin((*alpha)*uni+M_PI*(*alpha)/2.0)/pow(cos(uni)*cos(M_PI*(*alpha)/2.0),1.0/(*alpha));
        z1=pow(cos((1.0-*alpha)*uni-M_PI*(*alpha)/2.0)/exp_rand(),(1.0-*alpha)/(*alpha));
        y=x1*y1*z1; /* here input variables are pointa type*/
        if(unif_rand()<=exp(-(*b)*y))
        {
            rn[i]=y;
            i++;
        }
        
            }
    PutRNGstate();
}

double qdens(double m, double *lambda, double beta)
{
    
    return pow(m,*lambda-1.0)*exp(-beta/2.0*(m+1.0/m));
}




