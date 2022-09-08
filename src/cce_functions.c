#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <float.h>

void ctsubsampling(double *znum, double *ztime, int *frequency, int *nsparse, int *n, double *grid, double *result)
{
  int t, i, k;
   
  for(t = 0; t < *frequency; t++){
    for(i = 0; i < *nsparse; i++){
      /*k = 1;*/
      /*while((ztime[k]<=grid[i])&&(k<*n)){*/
      for(k = 1; k < *n; k++){
        /*k++;*/
        if(ztime[k] > grid[i]){
          break;
        }
      }
      result[(*nsparse)*t+i] += znum[k-1];
      grid[i]++;
    }
    
  }
}


void refreshsampling(int *Dim, int* I, double *Times, double *rtimes,
                     int *Lengths, int *D, double *MinTime, int *MinL,
                     int *Samplings)
{
  int d, i, J;
  double Tmp;
  
  for(d = 0; d < *Dim; d++) {
    /*while(Times[Lengths[d] * d + (I[d]+1)] <= rtimes[0]){*/
    while(Times[D[d] + (I[d]+1)] <= rtimes[0]){
      I[d]++;
      /*if((I[d]+1) >= Lengths[d + 1]){*/
      if((I[d]+1) >= Lengths[d]){
        break;
      }
    }
    Samplings[(*MinL) * d] = I[d] + 1;
  }
  
  i = 0;
  
  while(rtimes[i] < *MinTime) {
    
    rtimes[i + 1] = rtimes[i];
    
    for(d = 0; d < *Dim; d++) {
      Tmp = rtimes[i];
      J = I[d];
      /*while((J < (Lengths[d + 1]-1)) && (Tmp <= rtimes[i])) {*/
      while((J < (Lengths[d]-1)) && (Tmp <= rtimes[i])) {
        J++;
        /*Tmp = Times[Lengths[d] * d + J];*/
        Tmp = Times[D[d] + J];
      }
      if(Tmp > rtimes[i + 1])
        rtimes[i + 1] = Tmp;
    }
    
    i++;
      
    for(d = 0; d < *Dim; d++) {
      /*while(Times[Lengths[d] * d + (I[d]+1)] <= rtimes[i]){*/
      while(Times[D[d] + (I[d]+1)] <= rtimes[i]){
        I[d]++;
        /*if((I[d]+1) >= Lengths[d + 1]){*/
        if((I[d]+1) >= Lengths[d]){
          break;
        }
      }
      Samplings[(*MinL) * d + i] = I[d] + 1;
    }
  }
        
}


void refreshsamplingphy(int *Dim, int* I, double *Times, double *rtimes,
                        int *Lengths, int *D, double *MinTime, int *MinL,
                        int *Samplings, int *rNum)
{
  int d, i;
  double Tmp;
  
  for(d = 0; d < *Dim; d++){
    Samplings[(*MinL + 1) * d] = 1;
  }
  
  for(i = 0; rtimes[i] < *MinTime; i++) {
    rtimes[i + 1] = rtimes[i];
    for(d = 0; d < *Dim; d++) {
      /*while(I[d] < (Lengths[d + 1] - 1)){*/
      while(I[d] < (Lengths[d] - 1)){
        I[d]++;
        /*if(Times[Lengths[d] * d + I[d]] > rtimes[i]){*/
        if(Times[D[d] + I[d]] > rtimes[i]){
          Samplings[(*MinL + 1) * d + (i + 1)] = I[d] + 1;
          /*Tmp = Times[Lengths[d] * d + I[d]];*/
          Tmp = Times[D[d] + I[d]];
          if(rtimes[i + 1] < Tmp){
            rtimes[i + 1] = Tmp; 
          }
          break;
        }
      }
    }
  }
  
  *rNum = i + 1;
  
  for(d = 0; d < *Dim; d++) {
    /*while(I[d] < (Lengths[d + 1] -1)){*/
    while(I[d] < (Lengths[d] -1)){
      I[d]++;
      /*if(Times[Lengths[d] * d + I[d]] > rtimes[i]){*/
      if(Times[D[d] + I[d]] > rtimes[i]){
        Samplings[(*MinL + 1) * d + (i + 1)] = I[d] + 1;
        break;
      }
    }
  }
  
}


void bibsynchro(double *xtime, double *ytime, int *xlength, int *ylength,
                int *mu, int *w, int *q, int *r, int *Num)
{
  int I, mu0, w0, i;
  
  if(xtime[0] < ytime[0]){
    I = 0;
    while(xtime[I] < ytime[0]){
      I++;
      if(I >= (*xlength-1)){
        break;
      }
    }
    mu0 = I;
    if(ytime[0] < xtime[mu0]){
      q[0] = mu0;
    }else{
      q[0] = mu0 + 1;
    }
    r[0] = 1;
  }else if(xtime[0] > ytime[0]){
    I = 0;
    while(xtime[0] > ytime[I]){
      I++;
      if(I >= (*ylength-1)){
        break;
      }
    }
    w0 = I;
    q[0] = 1;
    if(xtime[0] < ytime[w0]){
      r[0] = w0;
    }else{
      r[0] = w0 + 1;
    }
  }else{
    q[0] = 1;
    r[0] = 1;
  }
  
  for(i = 0; (q[i] < (*xlength-1)) && (r[i] < (*ylength-1)); i++) {
    if(xtime[q[i]] < ytime[r[i]]){
      /*if(xtime[*xlength - 1] < ytime[r[i]]){*/
      if(xtime[*xlength - 2] < ytime[r[i]]){
        break;
      }
      I = q[i];
      while(xtime[I] < ytime[r[i]]){
        I++;
      }
      mu[i] = I;
      w[i] = r[i];
      if(ytime[r[i]] < xtime[mu[i]]){
        q[i+1] = mu[i];
      }else{
        q[i+1] = mu[i] + 1;
      }
      r[i+1] = r[i] + 1;
    }else if(xtime[q[i]] > ytime[r[i]]){
      /*if(xtime[q[i]] > ytime[*ylength - 1]){*/
      if(xtime[q[i]] > ytime[*ylength - 2]){
        break;
      }
      mu[i] = q[i];
      I = r[i];
      while(xtime[q[i]] > ytime[I]){
        I++;
      }
      w[i] = I;
      q[i+1] = q[i] + 1;
      if(xtime[q[i]] < ytime[w[i]]){
        r[i+1] = w[i];
      }else{
        r[i+1] = w[i] + 1;
      }
    }else{
      mu[i] = q[i];
      w[i] = r[i];
      q[i+1] = q[i] + 1;
      r[i+1] = r[i] + 1;
    }
    
  }
  
  mu[i] = q[i];
  w[i] = r[i];
  *Num = i + 1;
}


void HayashiYoshida(int *xlength, int *ylength, double *xtime, double *ytime,
                    double *dX, double *dY, double *value)
{
    int I, J;
    I = 0;
    J = 0;
    
    /* Checking Starting Point */
    while((I < (*xlength-1)) && (J < (*ylength-1))){
        if(xtime[I] >= ytime[J + 1]){
            J++;
        }else if(xtime[I + 1] <= ytime[J]){
            I++;
        }else{
            break;
        }
    }
    
    /* Main Component */
    while((I < (*xlength-1)) && (J < (*ylength-1))) {
        *value += dX[I] * dY[J];
        if(xtime[I + 1] > ytime[J + 1]){
            J++;
        }else if(xtime[I + 1] < ytime[J + 1]){
            I++;
        }else{
            I++;
            J++;
        }
    }
}


void pHayashiYoshida(int *kn, int *xlength, int *ylength, 
                     double *xtime, double *ytime, double *barX,
                     double *barY, double *value)
{
  int k, l, start, end;
  start = *kn;
  end = 0;
  
  for(k = 0; k < *xlength; k++){
    while((xtime[k] >= ytime[start]) && ((start - *kn) < (*ylength-1))) {
      start++;
    }
    while((xtime[k + *kn] > ytime[end + 1]) && (end < (*ylength-1))) {
      end++;
    }
    for(l = (start - *kn); l <= end; l++){
      *value += barX[k] * barY[l];
    }
  }
}


void msrc(int *M, int *N, double *xg, double *xl, double *ygamma, double *ylambda,
          double *result)
{
  int m, i;
  
  for(m = 0; m < *M; m++) {
    for(i = m; i < *N; i++){
      result[m] += (xg[i] - xl[i-m]) * (ygamma[i] - ylambda[i-m]);
    }
  }
  
}


void HYcrosscov(int *gridL, int *xL, int *yL, double *grid, double *xtime,
                double *ytime, double *tmptime, double *dX, double *dY, double *value)
{
    int i, j, I, J;
    
    for(i = 0; i < *gridL; i++){
        
        for(j = 0; j < *yL; j++){
            tmptime[j] = round(ytime[j] + grid[i]);
        }
        
        I = 0;
        J = 0;
        
        /* Checking Starting Point */
        while((I < (*xL-1)) && (J < (*yL-1))){
            if(round(xtime[I]) >= tmptime[J + 1]){
                J++;
            }else if(round(xtime[I + 1]) <= tmptime[J]){
                I++;
            }else{
                break;
            }
        }
        
        /* Main Component */
        while((I < (*xL-1)) && (J < (*yL-1))) {
            value[i] += dX[I] * dY[J];
            if(round(xtime[I + 1]) > tmptime[J + 1]){
                J++;
            }else if(round(xtime[I + 1]) < tmptime[J + 1]){
                I++;
            }else{
                I++;
                J++;
            }
        }
    }
}


void HYcrosscorr(int *gridL, int *xL, int *yL, double *grid, double *xtime,
                 double *ytime, double *tmptime, double *dX, double *dY,
                 double *xvol, double *yvol, double *value)
{
    int i, j, I, J;
    double A, B, C, s;
    
    for(i = 0; i < *gridL; i++){
        
        for(j = 0; j < *yL; j++){
            tmptime[j] = round(ytime[j] + grid[i]);
        }
        
        I = 0;
        J = 0;
        
        /* Checking Starting Point */
        while((I < (*xL-1)) && (J < (*yL-1))){
            if(round(xtime[I]) >= tmptime[J + 1]){
                J++;
            }else if(round(xtime[I + 1]) <= tmptime[J]){
                I++;
            }else{
                break;
            }
        }
        
        /* Main Component */
        while((I < (*xL-1)) && (J < (*yL-1))) {
            value[i] += dX[I] * dY[J];
            if(round(xtime[I + 1]) > tmptime[J + 1]){
                J++;
            }else if(round(xtime[I + 1]) < tmptime[J + 1]){
                I++;
            }else{
                I++;
                J++;
            }
        }
        
        /* Positive semi-definite correction */
        A = (*xvol)*(*xvol) + value[i]*value[i];
        B = (*xvol + *yvol)*value[i];
        C = (*yvol)*(*yvol) + value[i]*value[i];
        
        /* correct a bug (2016-06-30) */
        /* s = sqrt(A*C-B*B); */
        s = A*C-B*B;
        if(s > 0){
            s = sqrt(s);
        }else{
            s = 0;
        }
        
        
        value[i] = B/sqrt((A + s)*(C + s));
        
    }
}


void HYcrosscov2(int *gridL, int *xL, int *yL, double *grid, double *xtime,
                 double *ytime, double *dX, double *dY, double *value)
{
    int i, j, I, J;
    double *ytime2;
    
    for (j = 0; j < *xL; j++) {
        xtime[j] = round(xtime[j]);
    }
    
    for(i = 0; i < *gridL; i++){
        
        ytime2 = (double *)malloc(sizeof(double) * (*yL));
        
        for(j = 0; j < *yL; j++){
            ytime2[j] = round(ytime[j] + grid[i]);
        }
        
        I = 0;
        J = 0;
        
        /* Checking Starting Point */
        while((I < (*xL-1)) && (J < (*yL-1))){
            if(xtime[I] >= ytime2[J + 1]){
                J++;
            }else if(xtime[I + 1] <= ytime2[J]){
                I++;
            }else{
                break;
            }
        }
        
        /* Main Component */
        while((I < (*xL-1)) && (J < (*yL-1))) {
            value[i] += dX[I] * dY[J];
            if(xtime[I + 1] > ytime2[J + 1]){
                J++;
            }else if(xtime[I + 1] < ytime2[J + 1]){
                I++;
            }else{
                I++;
                J++;
            }
        }
        
        free(ytime2);
    }
}


void HYcrosscorr2(int *gridL, int *xL, int *yL, double *grid, double *xtime,
                  double *ytime, double *dX, double *dY,
                  double *xvol, double *yvol, double *value)
{
    int i, j, I, J;
    double A, B, C, s;
    double *ytime2;
    
    for (j = 0; j < *xL; j++) {
        xtime[j] = round(xtime[j]);
    }
    
    for(i = 0; i < *gridL; i++){
        
        ytime2 = (double *)malloc(sizeof(double) * (*yL));
        
        for(j = 0; j < *yL; j++){
            ytime2[j] = round(ytime[j] + grid[i]);
        }
        
        I = 0;
        J = 0;
        
        /* Checking Starting Point */
        while((I < (*xL-1)) && (J < (*yL-1))){
            if(xtime[I] >= ytime2[J + 1]){
                J++;
            }else if(xtime[I + 1] <= ytime2[J]){
                I++;
            }else{
                break;
            }
        }
        
        /* Main Component */
        while((I < (*xL-1)) && (J < (*yL-1))) {
            value[i] += dX[I] * dY[J];
            if(xtime[I + 1] > ytime2[J + 1]){
                J++;
            }else if(xtime[I + 1] < ytime2[J + 1]){
                I++;
            }else{
                I++;
                J++;
            }
        }
        
        /* Positive semi-definite correction */
        A = (*xvol)*(*xvol) + value[i]*value[i];
        B = (*xvol + *yvol)*value[i];
        C = (*yvol)*(*yvol) + value[i]*value[i];
        
        /* correct a bug (2016-06-30) */
        /* s = sqrt(A*C-B*B); */
        s = A*C-B*B;
        if(s > 0){
            s = sqrt(s);
        }else{
            s = 0;
        }
        
        
        value[i] = B/sqrt((A + s)*(C + s));
        
        free(ytime2);
    }
}



void hyavar(double *xtime, double *ytime, int *xlength, int *ylength,
            double *x, double *y, int *N, double *h, double *result)
{
    int I, mu0, w0, i, L, m;
    double dSigma11, dSigma12, dSigma22;
    int mu[*N], w[*N], q[*N], r[*N];
    double rtimes[*N], Sigma11[*N], Sigma12[*N], Sigma22[*N], H1[*N], H2[*N], H12[*N], H3[*N];
//    double rtimes[*N], Sigma11[*N], Sigma12[*N], Sigma22[*N], H1[*N], H2[*N], H12[*N], H3[*N], dS2[*N], dxdy[*N];
    
    Sigma11[0] = 0; Sigma12[0] = 0; Sigma22[0] = 0;
    
    if(xtime[0] < ytime[0]){
        rtimes[0] = ytime[0];
        I = 0;
        while(xtime[I] < ytime[0]){
            I++;
            if(I >= (*xlength-1)){
                break;
            }
        }
        mu0 = I;
        if(ytime[0] < xtime[mu0]){
            q[0] = mu0;
        }else{
            q[0] = mu0 + 1;
        }
        r[0] = 1;
    }else if(xtime[0] > ytime[0]){
        rtimes[0] = xtime[0];
        I = 0;
        while(xtime[0] > ytime[I]){
            I++;
            if(I >= (*ylength-1)){
                break;
            }
        }
        w0 = I;
        q[0] = 1;
        if(xtime[0] < ytime[w0]){
            r[0] = w0;
        }else{
            r[0] = w0 + 1;
        }
    }else{
        rtimes[0] = xtime[0];
        q[0] = 1;
        r[0] = 1;
    }
    
    for(i = 0; (q[i] < (*xlength - 1)) && (r[i] < (*ylength - 1)); i++) {
        H1[i] = 0; H2[i] = 0; H12[i] = 0; H3[i] = 0;
        if(xtime[q[i]] < ytime[r[i]]){
            /*if(xtime[*xlength - 1] < ytime[r[i]]){*/
            if(xtime[*xlength - 2] < ytime[r[i]]){
                break;
            }
            Sigma22[i] += pow(y[r[i]] - y[r[i] - 1], 2);
            H2[i] = pow(ytime[r[i]] - ytime[r[i] - 1], 2);
            I = q[i];
            while(xtime[I] < ytime[r[i]]){
                Sigma11[i] += pow(x[I] - x[I - 1], 2);
                H1[i] += pow(xtime[I] - xtime[I - 1], 2);
                I++;
            }
            mu[i] = I;
            w[i] = r[i];
            r[i+1] = r[i] + 1;
            if(ytime[r[i]] < xtime[mu[i]]){
                q[i+1] = mu[i];
            }else{
                q[i+1] = mu[i] + 1;
                Sigma11[i] += pow(x[mu[i]] - x[mu[i] - 1], 2);
                H1[i] += pow(xtime[mu[i]] - xtime[mu[i] - 1], 2);
            }
            H12[i] = H1[i] - pow(xtime[q[i]] - xtime[q[i] - 1], 2) + pow(xtime[q[i]] - rtimes[i], 2);
        }else if(xtime[q[i]] > ytime[r[i]]){
            /*if(xtime[q[i]] > ytime[*ylength - 1]){*/
            if(xtime[q[i]] > ytime[*ylength - 2]){
                break;
            }
            Sigma11[i] += pow(x[q[i]] - x[q[i] - 1], 2);
            H1[i] = pow(xtime[q[i]] - xtime[q[i] - 1], 2);
            mu[i] = q[i];
            I = r[i];
            while(xtime[q[i]] > ytime[I]){
                Sigma22[i] += pow(y[I] - y[I - 1], 2);
                H2[i] += pow(ytime[I] - ytime[I - 1], 2);
                I++;
            }
            w[i] = I;
            q[i+1] = q[i] + 1;
            if(xtime[q[i]] < ytime[w[i]]){
                r[i+1] = w[i];
            }else{
                r[i+1] = w[i] + 1;
                Sigma22[i] += pow(y[w[i]] - y[w[i] - 1], 2);
                H2[i] += pow(ytime[w[i]] - ytime[w[i] - 1], 2);
            }
            H12[i] = H2[i] - pow(ytime[r[i]] - ytime[r[i] - 1], 2) + pow(ytime[r[i]] - rtimes[i], 2);
        }else{
            mu[i] = q[i];
            w[i] = r[i];
            q[i+1] = q[i] + 1;
            r[i+1] = r[i] + 1;
            Sigma11[i] += pow(x[q[i]] - x[q[i] - 1], 2);
            Sigma22[i] += pow(y[r[i]] - y[r[i] - 1], 2);
            H1[i] = pow(xtime[q[i]] - xtime[q[i] - 1], 2);
            H2[i] = pow(ytime[r[i]] - ytime[r[i] - 1], 2);
            H12[i] = pow(xtime[q[i]] - rtimes[i], 2);
        }
        
        //dxdy[i] = (x[mu[i]] - x[q[i] - 1]) * (y[w[i]] - y[r[i] - 1]);
        Sigma12[i] += (x[mu[i]] - x[q[i] - 1]) * (y[w[i]] - y[r[i] - 1]);
        H3[i] = (xtime[mu[i]] - xtime[q[i] - 1]) * (ytime[w[i]] - ytime[r[i] - 1]);
        rtimes[i + 1] = fmin(xtime[mu[i]], ytime[w[i]]);
        //dS2[i] = pow(rtimes[i + 1] - rtimes[i], 2);
        
        Sigma11[i + 1] = Sigma11[i];
        Sigma22[i + 1] = Sigma22[i];
        Sigma12[i + 1] = Sigma12[i];
        
    }
    
    mu[i] = q[i];
    w[i] = r[i];
    //dxdy[i] = (x[mu[i]] - x[q[i] - 1]) * (y[w[i]] - y[r[i] - 1]);
    
    L = 0;
    
    for(m = 0; m < (i - 1); m++){
        
        while(rtimes[L + 1] <= (rtimes[m + 1] - *h)){
            L++;
        }
        
        dSigma11 = (Sigma11[m] - Sigma11[L])/(*h);
        dSigma12 = (Sigma12[m] - Sigma12[L])/(*h);
        dSigma22 = (Sigma22[m] - Sigma22[L])/(*h);
        
        /*result[0] += dxdy[m] * (dxdy[m] + 2 * dxdy[m + 1]) - 2 * dSigma12 * dSigma12 * dS2[m];*/
        result[0] += dSigma11 * dSigma22 * H3[m] + dSigma12 * dSigma12 * (H1[m] + H2[m] - H12[m]);
        result[1] += 2 * dSigma11 * dSigma12 * H1[m];
        result[2] += 2 * dSigma12 * dSigma22 * H2[m];
        result[3] += 2 * dSigma12 * dSigma12 * H12[m];
        
    }
    
}


void hycrossavar(double *grid, double *xtime, double *ytime, int *gridL, int *xlength, int *ylength, double *x, double *y,
                 int *N, double *h, double *v1, double *v2, double *d11, double *d22, double *d12,
                 double * d23, double *d31, double *covar, double *corr)
{
    int I, mu0, w0, i, L, m, j, k;
    double dSigma11, dSigma12, dSigma22, avar1, avar2, avar3, avar4;
    int mu[*N], w[*N], q[*N], r[*N];
    double rtimes[*N], Sigma11[*N], Sigma12[*N], Sigma22[*N], H1[*N], H2[*N], H12[*N], H3[*N], tmptime[*ylength];
//    double rtimes[*N], Sigma11[*N], Sigma12[*N], Sigma22[*N], H1[*N], H2[*N], H12[*N], H3[*N], dS2[*N], dxdy[*N], tmptime[*ylength];
    
    for(i = 0; i < *xlength; i++){
      xtime[i] = round(xtime[i]);
    }
    
    for (k = 0; k < *gridL; k++) {
        
        for(j = 0; j < *ylength; j++){
            tmptime[j] = round(ytime[j] + grid[k]);
        }
        
        Sigma11[0] = 0; Sigma12[0] = 0; Sigma22[0] = 0;
        
        if(xtime[0] < tmptime[0]){
            rtimes[0] = tmptime[0];
            I = 0;
            while(xtime[I] < tmptime[0]){
                I++;
                if(I >= (*xlength-1)){
                    break;
                }
            }
            mu0 = I;
            if(tmptime[0] < xtime[mu0]){
                q[0] = mu0;
            }else{
                q[0] = mu0 + 1;
            }
            r[0] = 1;
        }else if(xtime[0] > tmptime[0]){
            rtimes[0] = xtime[0];
            I = 0;
            while(xtime[0] > tmptime[I]){
                I++;
                if(I >= (*ylength-1)){
                    break;
                }
            }
            w0 = I;
            q[0] = 1;
            if(xtime[0] < tmptime[w0]){
                r[0] = w0;
            }else{
                r[0] = w0 + 1;
            }
        }else{
            rtimes[0] = xtime[0];
            q[0] = 1;
            r[0] = 1;
        }
        
        for(i = 0; (q[i] < (*xlength - 1)) && (r[i] < (*ylength - 1)); i++) {
            H1[i] = 0; H2[i] = 0; H12[i] = 0; H3[i] = 0;
            if(xtime[q[i]] < tmptime[r[i]]){
                /*if(xtime[*xlength - 1] < tmptime[r[i]]){*/
                if(xtime[*xlength - 2] < tmptime[r[i]]){
                    break;
                }
                Sigma22[i] += pow(y[r[i]] - y[r[i] - 1], 2);
                H2[i] = pow(tmptime[r[i]] - tmptime[r[i] - 1], 2);
                I = q[i];
                while(xtime[I] < tmptime[r[i]]){
                    Sigma11[i] += pow(x[I] - x[I - 1], 2);
                    H1[i] += pow(xtime[I] - xtime[I - 1], 2);
                    I++;
                }
                mu[i] = I;
                w[i] = r[i];
                r[i+1] = r[i] + 1;
                if(tmptime[r[i]] < xtime[mu[i]]){
                    q[i+1] = mu[i];
                }else{
                    q[i+1] = mu[i] + 1;
                    Sigma11[i] += pow(x[mu[i]] - x[mu[i] - 1], 2);
                    H1[i] += pow(xtime[mu[i]] - xtime[mu[i] - 1], 2);
                }
                H12[i] = H1[i] - pow(xtime[q[i]] - xtime[q[i] - 1], 2) + pow(xtime[q[i]] - rtimes[i], 2);
            }else if(xtime[q[i]] > tmptime[r[i]]){
                /*if(xtime[q[i]] > tmptime[*ylength - 1]){*/
                if(xtime[q[i]] > tmptime[*ylength - 2]){
                    break;
                }
                Sigma11[i] += pow(x[q[i]] - x[q[i] - 1], 2);
                H1[i] = pow(xtime[q[i]] - xtime[q[i] - 1], 2);
                mu[i] = q[i];
                I = r[i];
                while(xtime[q[i]] > tmptime[I]){
                    Sigma22[i] += pow(y[I] - y[I - 1], 2);
                    H2[i] += pow(tmptime[I] - tmptime[I - 1], 2);
                    I++;
                }
                w[i] = I;
                q[i+1] = q[i] + 1;
                if(xtime[q[i]] < tmptime[w[i]]){
                    r[i+1] = w[i];
                }else{
                    r[i+1] = w[i] + 1;
                    Sigma22[i] += pow(y[w[i]] - y[w[i] - 1], 2);
                    H2[i] += pow(tmptime[w[i]] - tmptime[w[i] - 1], 2);
                }
                H12[i] = H2[i] - pow(tmptime[r[i]] - tmptime[r[i] - 1], 2) + pow(tmptime[r[i]] - rtimes[i], 2);
            }else{
                mu[i] = q[i];
                w[i] = r[i];
                q[i+1] = q[i] + 1;
                r[i+1] = r[i] + 1;
                Sigma11[i] += pow(x[q[i]] - x[q[i] - 1], 2);
                Sigma22[i] += pow(y[r[i]] - y[r[i] - 1], 2);
                H1[i] = pow(xtime[q[i]] - xtime[q[i] - 1], 2);
                H2[i] = pow(tmptime[r[i]] - tmptime[r[i] - 1], 2);
                H12[i] = pow(xtime[q[i]] - rtimes[i], 2);
            }
            
            //dxdy[i] = (x[mu[i]] - x[q[i] - 1]) * (y[w[i]] - y[r[i] - 1]);
            Sigma12[i] += (x[mu[i]] - x[q[i] - 1]) * (y[w[i]] - y[r[i] - 1]);
            H3[i] = (xtime[mu[i]] - xtime[q[i] - 1]) * (tmptime[w[i]] - tmptime[r[i] - 1]);
            rtimes[i + 1] = fmin(xtime[mu[i]], tmptime[w[i]]);
            //dS2[i] = pow(rtimes[i + 1] - rtimes[i], 2);
            
            Sigma11[i + 1] = Sigma11[i];
            Sigma22[i + 1] = Sigma22[i];
            Sigma12[i + 1] = Sigma12[i];
            
        }
        
        mu[i] = q[i];
        w[i] = r[i];
        //dxdy[i] = (x[mu[i]] - x[q[i] - 1]) * (y[w[i]] - y[r[i] - 1]);
        
        L = 0;
        avar1 = 0; avar2 = 0; avar3 = 0; avar4 = 0;
        
        for(m = 0; m < (i - 1); m++){
            
            while(rtimes[L + 1] <= (rtimes[m + 1] - *h)){
                L++;
            }
            
            dSigma11 = (Sigma11[m] - Sigma11[L])/(*h);
            dSigma12 = (Sigma12[m] - Sigma12[L])/(*h);
            dSigma22 = (Sigma22[m] - Sigma22[L])/(*h);
            
            /*result[0] += dxdy[m] * (dxdy[m] + 2 * dxdy[m + 1]) - 2 * dSigma12 * dSigma12 * dS2[m];*/
            avar1 += dSigma11 * dSigma22 * H3[m] + dSigma12 * dSigma12 * (H1[m] + H2[m] - H12[m]);
            avar2 += 2 * dSigma11 * dSigma12 * H1[m];
            avar3 += 2 * dSigma12 * dSigma22 * H2[m];
            avar4 += 2 * dSigma12 * dSigma12 * H12[m];
            
        }
        
        covar[k] = avar1;
        
        corr[k] = (*v1) * d11[k] + avar1 + (*v2) * d22[k] + avar2 * d12[k] + avar3 * d23[k] + avar4 * d31[k];
        
    }
    
}


void krprod(double *X, int *nrow, int *ncol, double *result)
{
    int i, j, k, tmpcol;
    
    for (k = 0; k < *ncol; k++) {
        
        tmpcol = k * (*nrow);
        
        for (i = 0; i < *nrow; i++) {
            
            for (j = i; j < *nrow; j++) {
                
                result[tmpcol * (*nrow) + i * (*nrow) + j] = X[tmpcol + i] * X[tmpcol + j];
                result[tmpcol * (*nrow) + j * (*nrow) + i] = result[tmpcol * (*nrow) + i * (*nrow) + j];
                
            }
            
        }
        
    }
    
}
