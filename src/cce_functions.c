#include <Rinternals.h>
#include <math.h>

void ctsubsampling(double *znum, double *ztime, int *frequency, int *nsparse,
                   int *n, double *grid, double *result)
{
  int t, i, k;
  
  for(t = 0; t < *frequency; t++){
    k = 1;
    for(i = 0; i < *nsparse; i++){
      while((ztime[k]<=grid[i])&&(k < *n)){
        k++;
      }
      result[(*nsparse)*t+i] += znum[k-1];
      grid[i]++;
    }
    
  }
}


void refreshsampling(int *Dim, int* I, double *Times, double *rtimes,
                     int *Lengths, double *MinTime, int *MinL,
                     int *Samplings)
{
  int d, i, J;
  double Tmp;
  
  for(d = 0; d < *Dim; d++) {
    while(Times[Lengths[d] * d + (I[d]+1)] <= rtimes[0]){
      I[d]++;
      if((I[d]+1) >= Lengths[d + 1]){
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
      while((J < (Lengths[d + 1]-1)) && (Tmp <= rtimes[i])) {
        J++;
        Tmp = Times[Lengths[d] * d + J];
      }
      if(Tmp > rtimes[i + 1])
        rtimes[i + 1] = Tmp;
    }
    
    i++;
      
    for(d = 0; d < *Dim; d++) {
      while(Times[Lengths[d] * d + (I[d]+1)] <= rtimes[i]){
        I[d]++;
        if((I[d]+1) >= Lengths[d + 1]){
          break;
        }
      }
      Samplings[(*MinL) * d + i] = I[d] + 1;
    }
  }
        
}


void refreshsamplingphy(int *Dim, int* I, double *Times, double *rtimes,
                        int *Lengths, double *MinTime, int *MinL,
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
      while(I[d] < (Lengths[d + 1] - 1)){
        I[d]++;
        if(Times[Lengths[d] * d + I[d]] > rtimes[i]){
          Samplings[(*MinL + 1) * d + (i + 1)] = I[d] + 1;
          Tmp = Times[Lengths[d] * d + I[d]];
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
    while(I[d] < (Lengths[d + 1] -1)){
      I[d]++;
      if(Times[Lengths[d] * d + I[d]] > rtimes[i]){
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
      if(xtime[*xlength - 1] < ytime[r[i]]){
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
      if(xtime[q[i]] > ytime[*ylength - 1]){
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
      tmptime[j] = ytime[j] + grid[i];
    }
    
    I = 0;
    J = 0;
    
    /* Checking Starting Point */
    while((I < (*xL-1)) && (J < (*yL-1))){
        if(xtime[I] >= tmptime[J + 1]){
            J++;
        }else if(xtime[I + 1] <= tmptime[J]){
            I++;
        }else{
            break;
        }
    }
    
    /* Main Component */
    while((I < (*xL-1)) && (J < (*yL-1))) {
        value[i] += dX[I] * dY[J];
        if(xtime[I + 1] > tmptime[J + 1]){
            J++;
        }else if(xtime[I + 1] < tmptime[J + 1]){
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
      tmptime[j] = ytime[j] + grid[i];
    }
    
    I = 0;
    J = 0;
    
    /* Checking Starting Point */
    while((I < (*xL-1)) && (J < (*yL-1))){
        if(xtime[I] >= tmptime[J + 1]){
            J++;
        }else if(xtime[I + 1] <= tmptime[J]){
            I++;
        }else{
            break;
        }
    }
    
    /* Main Component */
    while((I < (*xL-1)) && (J < (*yL-1))) {
        value[i] += dX[I] * dY[J];
        if(xtime[I + 1] > tmptime[J + 1]){
            J++;
        }else if(xtime[I + 1] < tmptime[J + 1]){
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
    
    s = sqrt(A*C-B*B);
    
    value[i] = B/sqrt((A + s)*(C + s));
    
  }
}
