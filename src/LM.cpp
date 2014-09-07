#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/LM.hpp"

//
// Class creator
//
template <typename real>
LM<real>::LM(const int m) {
  zavg = new real[m+1];
  alpha0 = new real[m+1];
  a = new real[m+1];
  aold = new real[m+1];
  H = new real[(m+1)*(m+1)];
  dmlnL = new real[m+1];
  x = new real[m+1];
  dmlnLold = new real[m+1];
  Hdx_help = new real[m+1];
  y_help = new real[m+1];
  r_help = new real[m+1];
  S = new real[(m+1)*(m+1)];
  K = new real[(m+1)];
  K_help = new real[(m+1)*(m+1)];
  KK_help = new real[(m+1)*(m+1)];
  x_help = new real[3*(m+1)];
  p_help = new real[(m+1)*(m+1)];
}

//
// Class destructor
//
template <typename real>
LM<real>::~LM() {
  delete [] zavg;
  delete [] alpha0;
  delete [] a;
  delete [] aold;
  delete [] H;
  delete [] dmlnL;
  delete [] x;
  delete [] dmlnLold;
  delete [] Hdx_help;
  delete [] y_help;
  delete [] r_help;
  delete [] S;
  delete [] K;
  delete [] K_help;
  delete [] KK_help;
  delete [] x_help;
  delete [] p_help;
}

template <typename real>
real LM<real>::randomf(real a, real b){
  real aa;
  aa = ((real)( rand()%10001))/(real)(10000.0);
  return a+(b-a)*aa;
}

//
// a[m+1], z[M]
//
template <typename real>
inline real LM<real>::rxn_coor(const int m, const int *cv, const real *a, const real *z){
  real rc = a[0];
  for(int i=0;i < m;i++) {
    rc += z[cv[i]]*a[i+1];
  }
  return rc;
}

//
// a[m+1], zA[nalist][M], zB[nblist][M]
//
template <typename real>
real LM<real>::mlnL(const int m, const int *cv, const int nalist, const int nblist, const int M,
	       const real *a, const real *zA, const real *zB){
  real lnl = ((real)(nalist + nblist))*log(0.5);

  for(int i=0;i < nalist;i++){
    real rk = rxn_coor(m, cv, a, &zA[i*M]);
    lnl += log((real)1 - tanh(rk));
  }

  for(int i=0;i < nblist;i++){
    real rk = rxn_coor(m, cv, a, &zB[i*M]);
    lnl += log((real)1 + tanh(rk));
  }

  return (-lnl);
}

//
// a[m+1], zA[nalist][M], zB[nblist][M], dmlnL[m+1]
//
template <typename real>
real LM<real>::grad(const int m, const int *cv, const int nalist, const int nblist, const int M,
		const real *a, const real *zA, const real *zB, real *dmlnL){

  //real lnl = -mlnL(m, cv, nalist, nblist, M, a, zA, zB);
  real lnl = ((real)(nalist + nblist))*log((real)0.5);

  for(int i=0;i <= m;i++) {
    dmlnL[i] = (real)0;
  }

#pragma simd
  for(int i=0;i < nalist;i++) {
    real rk = rxn_coor(m, cv, a, &zA[i*M]);
    real tanh_rk = tanh(rk);
    lnl += log((real)1 - tanh_rk);
    real C = (real)1 + tanh_rk;
    dmlnL[0] += C;
    for(int j=0;j < m;j++){
      dmlnL[j+1] += C*zA[i*M + cv[j]];
    }
  }

#pragma simd
  for(int i=0;i < nblist;i++){
    real rk = rxn_coor(m, cv, a, &zB[i*M]);
    real tanh_rk = tanh(rk);
    lnl += log((real)1 + tanh_rk);
    real C = (real)1 - tanh_rk;
    dmlnL[0] -= C;
    for(int j=0;j < m;j++) {
      dmlnL[j+1] -= C*zB[i*M + cv[j]];
    }
  }

  return (-lnl);
}

//
// v[m+1], u[m+1]
//
template <typename real>
real LM<real>::dot(const int m, const real *v, const real *u){
  int k;
  real dp;
  dp = (real)0;
  for(k=0;k<=m;k=k+1){
    dp = dp+v[k]*u[k];
  }
  return dp;
}

//
// dx[m+1], g_new[m+1], g_old[m+1], h[m+1][m+1]
//
// Help arrays: Hdx[1+m], y[1+m], r[1+m]
//
//
template <typename real>
void LM<real>::BFGS_update(const bool debug, const int m, const real *dx, const real *g_new, const real *g_old,
		     real *h, real *Hdx, real *y, real *r) {

  if (debug) printf("\n  DAMPED BFGS HESSIAN UPDATE p.541 NOCEDAL & WRIGHT");
  if (debug) printf("\n    SECANT");
  for(int j=0;j <= m;j++){
    y[j] = g_new[j] - g_old[j];
  }
  if (debug) printf("  H.dx:");
  for(int i=0;i <= m;i++){
    Hdx[i] = dot(m, &h[i*(m+1)], dx);
    if (debug) printf(" %.2f",Hdx[i]);
  }
  real dxHdx = dot(m, dx, Hdx);
  real dxy = dot(m, dx, y);
  real theta;
  if( dxy > ((real)0.2)*dxHdx ){
    theta = (real)1;
  }else{
    theta = ((real)0.8)*dxHdx/(dxHdx-dxy);
  }
  if (debug) printf("  theta: %.3f", theta);
  for(int i=0;i <= m;i++){
    r[i] = theta*y[i]+((real)1-theta)*Hdx[i];
  }
  real dxr = dot(m, r, dx);
  for(int i=0;i <= m;i++){
    for(int j=0;j <= m;j++){
      h[i*(m+1)+j] += -Hdx[i]*Hdx[j]/dxHdx + r[i]*r[j]/dxr;
    }
  }
}

//
// x[1+m]
//
template <typename real>
real LM<real>::normalize(const int m, real *x){
  real numer;
  int j;
  numer = (real)0;
  for(j=0;j<=m;j=j+1){
    numer = numer + x[j]*x[j];
  }
  numer = sqrt(numer);
  for(j=0;j<=m;j=j+1){
    x[j] = x[j]/numer;
  }
  return numer;
}

//
// vector[1+m], hmwc[1+m][1+m]
// Helpers: p[1+m][1+m], K[1+m][1+m];
//
template <typename real>
void LM<real>::projectfrommatrix(const int m, const real *vector, real *hmwc,
			   real *p, real *K){

  for(int i=0;i<=m;i=i+1){
    for(int j=0;j<=m;j=j+1){
      if(i == j){
        p[i*(m+1)+j] = (real)1 - vector[i]*vector[j];
      }
      else{
        p[i*(m+1)+j] = - vector[i]*vector[j];
      }
    }
  }
  for(int i=0;i<=m;i=i+1){
    for(int j=0;j<=m;j=j+1){
      K[i*(m+1)+j] = (real)0;
      for(int k=0;k<=m;k=k+1){
        K[i*(m+1)+j] += hmwc[i*(m+1)+k]*p[k*(m+1)+j];
      }
    }
  }
  for(int i=0;i<=m;i=i+1){
    for(int j=0;j<=m;j=j+1){
      hmwc[i*(m+1)+j] = (real)0;
      for(int k=0;k<=m;k=k+1){
        hmwc[i*(m+1)+j] += p[i*(m+1)+k]*K[k*(m+1)+j];
      }
    }
  }
}

int mod(int a, int b)
{
  int r = a % b;
  return r < 0 ? r + b : r;
}

//
// hmwc[1+m][1+m], smwc[1+m][1+m], w2[1+m]
// Helpers: K[1+m][1+m], x[3][1+m], p[1+m][1+m], KK[1+m][1+m]
//
template <typename real>
void LM<real>::diagonalize(const bool debug, const int m, const real *hmwc, real *smwc, real *w2,
		     real *K, real *x, real *KK, real *p) {

  int ndiag = 1000;

  for(int j=0;j <= m;j++){
    for(int k=0;k <= m;k++){
      K[j*(m+1)+k] = hmwc[j*(m+1)+k];
    }
  }
  real tolerance = ((real)0.000001);
  if (debug) printf("\n  DIAGONALIZING:");
  fflush(stdout);
  for(int i=0;i <= m;i++){
    for(int j=0;j <= m;j++){
      x[0*(m+1)+j]=randomf(((real)(-1.0)), ((real)1.0));
      x[1*(m+1)+j] = (real)0;
      x[2*(m+1)+j] = (real)0;
    }
    normalize(m, &x[0*(m+1)]);
    int test = 0;
    for(int k=1;test == 0;k++){
      for(int j=0;j <= m;j++){
        x[mod(k,3)*(m+1)+j] = dot(m, &K[j*(m+1)], &x[mod(k-1,3)*(m+1)]);
      }
      normalize(m, &x[mod(k,3)*(m+1)]);
      real relerror = (real)0;
      for(int j=0;j <= m;j++){
        real diff = x[mod(k,3)*(m+1)+j];
	diff -= x[mod(k-2,3)*(m+1)+j];
        relerror = relerror + fabs(diff);
      }
      if((relerror < tolerance && k > 100) || (k == ndiag)){
        for(int j=0;j <= m;j++){
          x[((k+1)%3)*(m+1)+j] = dot(m, &K[j*(m+1)], &x[(k%3)*(m+1)]);
        }
        w2[i] = dot(m, &x[(k%3)*(m+1)], &x[((k+1)%3)*(m+1)]);
	if (debug) printf(" (%.2f) ",w2[i]);
        for(int j=0;j <= m;j++){
          smwc[i*(m+1)+j] = x[(k%3)*(m+1)+j];
        }
        projectfrommatrix(m, &smwc[i*(m+1)], K, p, KK);
        test = 1;
      }
    }
  }
  if (debug) printf("\n");
}

//
// Performs Likelihood Maximization calculation for m collective variables cv[0...m-1] = {0...M-1}
// zA[nalist][M]
// zB[nblist][M]
//
// alnLmax[m+2]
//
// crit_move = 0.001; crit_grad = 0.01; crit_dlnL = 0.0001; maxsize = 0.1;
//
template <typename real>
void LM<real>::calc_lm(const bool debug, const int m, const int *cv, const int nalist,
		       const int nblist, const int M, const real *zA, const real *zB,
		       const real crit_move, const real crit_grad, const real crit_dlnL, const real maxsize,
		       real *alnLmax) {

  for(int j=0;j <= m;j++) zavg[j] = (real)0;

  for(int i=0;i < nalist;i++){
    for(int j=0;j < m;j++){
      int cvj = cv[j];
      real zAval = zA[i*M + cvj];
      zavg[j] += zAval; //zA[i*M + cv[j-1]]; //zA[i][j];
    }
  }

  for(int j=0;j < m;j++){
    zavg[j] /= ((real)nalist);
  }

  /*****  SCREEN FOR STARTING VALUES  *****/


  real mlnl, mlnlold;

  if (debug) printf("\n  GET INITIAL VALUES BY RANDOM NUMBERS");
  for (int j=0;j <= m;j++) a[j] = (real)0;
  for (int j=0;j <= m;j++) alpha0[j] = (real)0;
  for(int i=1;i <= 16;i++){
    a[1] = randomf(((real)(-2.0)), ((real)2.0));
    a[0] = -a[1]*zavg[m];
    mlnl = mlnL(m, cv, nalist, nblist, M, a, zA, zB);
    if( mlnl < mlnL(m, cv, nalist, nblist, M, alpha0, zA, zB) ){
      for(int j=0;j <= 1;j++){
	alpha0[j] = a[j];
	if (debug) printf("\n improvement %lf a: %.3f %.3f",mlnl, a[0], a[1]);
      }
    }
  }

  //alpha0[0] = 0.0;
  //alpha0[1] = 0.0;

  if (debug) printf("\n  initial a: ");
  for(int j=0;j <= m;j++){
    a[j] = alpha0[j];
    if (debug) printf(" %.3f ", a[j]);
    for(int k=0;k <= m;k++){
      H[j*(m+1)+k] = (real)0;
    }
    H[j*(m+1)+j] = (real)1 + randomf((real)0, ((real)0.1));
  }

  int test = 0;
  int i;
  for(i=1;test==0;i++){
    if(i != 1){
      mlnlold = mlnl;
      for(int j=0;j<=m;j++){
	dmlnLold[j] = dmlnL[j];
      }
    }
    real mlnl = grad(m, cv, nalist, nblist, M, a, zA, zB, dmlnL);
    if (isnan(mlnl)) {
      printf("mlnl NaN\n");
      //exit(1);
      for(int j=0;j<=m;j++){
	alnLmax[j] = ((real)(-1.0e10));
      }
      return;
    }
    if (debug) printf("\n dmlnL = %f %f", dmlnL[0], dmlnL[1]);
    if(i != 1){
      BFGS_update(debug, m, x, dmlnL, dmlnLold, H, Hdx_help, y_help, r_help);
    }
    diagonalize(debug, m, H, S, K, K_help, x_help, KK_help, p_help);
    
    // What the heck is this supposed to do?
    //double dlnLsize = dmlnL - dmlnLold;
    real dlnLsize = (real)20.0;

    real gradsize = sqrt(dot(m, dmlnL, dmlnL));
    if (debug) printf("\n\n  gradients: %d  mlnL: %f  gradsize: %f", i, mlnl, gradsize);
    fflush(stdout);
    for(int j=0;j <= m;j++){
      x[j] = (real)0;
      for(int k=0;k <= m;k++){
	x[j] = x[j] - (dot(m, &S[k*(m+1)], dmlnL)/K[k])*S[k*(m+1)+j];
      }
    }
    real movesize = sqrt(dot(m,x,x));
    if (debug) printf("\n  RAW MOVESIZE/(B): %f", movesize);
    if(movesize > maxsize){
      for(int j=0;j <= m;j++){
	x[j] = x[j]*maxsize/movesize;
      }
    }      
    if (debug) printf("\n  NEW MOVESIZE/(B): %f  ALPHA: ", sqrt(dot(m,x,x)));
    for(int j=0;j <= m;j++){
      aold[j] = a[j];
      a[j] = a[j] + x[j];
      if (debug) printf(" %.3f", a[j]);
    }
    if (debug) printf("\n  CONVERGENCE TESTING");
    if (debug) printf("\n   movesize:tolerance  %f::%f", movesize, crit_move);
    if (debug) printf("\n   gradsize:tolerance  %f::%f", gradsize, crit_grad);
    if (debug) printf("\n   deltaE:|tolerance|  %f::%f", dlnLsize, crit_dlnL);
    if(movesize < crit_move && gradsize < crit_grad){
      test = 1;
    }
    if(movesize < crit_move && fabs(dlnLsize) < crit_dlnL){
      test = 1;
    }
    if(fabs(dlnLsize) < crit_dlnL && gradsize < crit_grad){
      test = 1;
    }
    if(test == 1){
      if (debug) printf("\n  CONVERGED ");
      for(int j=0;j<=m;j++){
	alnLmax[j] = a[j];
      }
      alnLmax[m+1] = -mlnl;
    }

  }

  //std::cout << "" << std::endl;

}

template class LM<float>;
template class LM<double>;
