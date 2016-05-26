#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "R.h"
#include "Mt.h"


void rpts(int *x, double *alpha, double *a, double *b, double *rn);
double rps(double *alpha, double *a, double *b);
double rexp();

void rpts(int *x, double *alpha, double *a, double *b, double *rn)
{
    int i=0;
    double y;
    init_genrand((unsigned)time(NULL));
    while(i<*x){
        y=rps(alpha,a,b); /* here input variables are pointa type*/
        if(genrand_real3()<=exp(-(*b)*y))
        {
            rn[i]=y;
            i++;
        }
        
            }
}



double rps(double *alpha, double *a, double *b)
{
    double x1,y1,z1,uni;
    uni=-M_PI/2.0+M_PI*genrand_real3();
    x1=pow((*a)*tgamma(1.0-*alpha)*cos(M_PI*(*alpha)/(2.0))/(*alpha),1.0/(*alpha));
    y1=sin((*alpha)*uni+M_PI*(*alpha)/2.0)/pow(cos(uni)*cos(M_PI*(*alpha)/2.0),1.0/(*alpha));
    z1=pow(cos((1.0-*alpha)*uni-M_PI*(*alpha)/2.0)/rexp(),(1.0-*alpha)/(*alpha));
    
    return x1*y1*z1;
}

double rexp()
{
    return -log(genrand_real3());
}


