
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LM.hpp"

// Class creator
LM::LM(const int m) {
  zavg = new double[m+1];
  alpha0 = new double[m+1];
  a = new double[m+1];
  aold = new double[m+1];
  H = new double[(m+1)*(m+1)];
  dmlnL = new double[m+1];
  x = new double[m+1];
  dmlnLold = new double[m+1];
  Hdx_help = new double[m+1];
  y_help = new double[m+1];
  r_help = new double[m+1];
  S = new double[(m+1)*(m+1)];
  K = new double[(m+1)];
  K_help = new double[(m+1)*(m+1)];
  KK_help = new double[(m+1)*(m+1)];
  x_help = new double[3*(m+1)];
  p_help = new double[(m+1)*(m+1)];
}

// Class destructor
LM::~LM() {
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

double LM::randomf(double a, double b){
  double aa;
  aa = ((double)( rand()%10001))/10000.0;
  return a+(b-a)*aa;
}

//
// a[m+1], z[M]
//
double LM::rxn_coor(const int m, const int *cv, const double *a, const double *z){
  double rc = a[0];
  for(int i=0;i < m;i++) {
    rc += z[cv[i]]*a[i+1];
  }
  return rc;
}

//
// a[m+1], zA[nalist][M], zB[nblist][M]
//
double LM::mlnL(const int m, const int *cv, const int nalist, const int nblist, const int M,
	       const double *a, const double *zA, const double *zB){
  double lnl = ((double)(nalist + nblist))*log(0.5);

  for(int i=0;i < nalist;i++){
    double rk = rxn_coor(m, cv, a, &zA[i*M]);
    lnl += log(1.0 - tanh(rk));
  }

  for(int i=0;i < nblist;i++){
    double rk = rxn_coor(m, cv, a, &zB[i*M]);
    lnl += log(1.0 + tanh(rk));
  }

  return (-lnl);
}

//
// a[m+1], zA[nalist][M], zB[nblist][M], dmlnL[m+1]
//
double LM::grad(const int m, const int *cv, const int nalist, const int nblist, const int M,
		const double *a, const double *zA, const double *zB, double *dmlnL){

  double lnl = -mlnL(m, cv, nalist, nblist, M, a, zA, zB);

  for(int i=0;i <= m;i++) {
    dmlnL[i] = 0.0;
  }

  for(int i=0;i < nalist;i++) {
    double rk = rxn_coor(m, cv, a, &zA[i*M]);
    double C = 1.0 + tanh(rk);
    dmlnL[0] += C;
    for(int j=0;j < m;j++){
      dmlnL[j+1] += C*zA[i*M + cv[j]];
    }
  }

  for(int i=0;i < nblist;i++){
    double rk = rxn_coor(m, cv, a, &zB[i*M]);
    double C = 1.0 - tanh(rk);
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
double LM::dot(const int m, const double *v, const double *u){
  int k;
  double dp;
  dp = 0.0;
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
void LM::BFGS_update(const bool debug, const int m, const double *dx, const double *g_new, const double *g_old,
		     double *h, double *Hdx, double *y, double *r) {

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
  double dxHdx = dot(m, dx, Hdx);
  double dxy = dot(m, dx, y);
  double theta;
  if( dxy > 0.2*dxHdx ){
    theta = 1.0;
  }else{
    theta = 0.8*dxHdx/(dxHdx-dxy);
  }
  if (debug) printf("  theta: %.3f", theta);
  for(int i=0;i <= m;i++){
    r[i] = theta*y[i]+(1.0-theta)*Hdx[i];
  }
  double dxr = dot(m, r, dx);
  for(int i=0;i <= m;i++){
    for(int j=0;j <= m;j++){
      h[i*(m+1)+j] += -Hdx[i]*Hdx[j]/dxHdx + r[i]*r[j]/dxr;
    }
  }
}

//
// x[1+m]
//
double LM::normalize(const int m, double *x){
  double numer;
  int j;
  numer = 0.0;
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
void LM::projectfrommatrix(const int m, const double *vector, double *hmwc,
			   double *p, double *K){

  for(int i=0;i<=m;i=i+1){
    for(int j=0;j<=m;j=j+1){
      if(i == j){
        p[i*(m+1)+j] = 1.0 - vector[i]*vector[j];
      }
      else{
        p[i*(m+1)+j] = - vector[i]*vector[j];
      }
    }
  }
  for(int i=0;i<=m;i=i+1){
    for(int j=0;j<=m;j=j+1){
      K[i*(m+1)+j] = 0.0;
      for(int k=0;k<=m;k=k+1){
        K[i*(m+1)+j] += hmwc[i*(m+1)+k]*p[k*(m+1)+j];
      }
    }
  }
  for(int i=0;i<=m;i=i+1){
    for(int j=0;j<=m;j=j+1){
      hmwc[i*(m+1)+j] = 0.0;
      for(int k=0;k<=m;k=k+1){
        hmwc[i*(m+1)+j] += p[i*(m+1)+k]*K[k*(m+1)+j];
      }
    }
  }
}

//
// hmwc[1+m][1+m], smwc[1+m][1+m], w2[1+m]
// Helpers: K[1+m][1+m], x[3][1+m], p[1+m][1+m], KK[1+m][1+m]
//
void LM::diagonalize(const bool debug, const int m, const double *hmwc, double *smwc, double *w2,
		     double *K, double *x, double *KK, double *p) {

  int ndiag = 1000;

  for(int j=0;j <= m;j++){
    for(int k=0;k <= m;k++){
      K[j*(m+1)+k] = hmwc[j*(m+1)+k];
    }
  }
  double tolerance = .000001;
  if (debug) printf("\n  DIAGONALIZING:");
  fflush(stdout);
  for(int i=0;i <= m;i++){
    for(int j=0;j <= m;j++){
      x[0*(m+1)+j]=randomf(-1.0,1.0);
      x[1*(m+1)+j] = 0.0;
      x[2*(m+1)+j] = 0.0;
    }
    normalize(m, &x[0*(m+1)]);
    int test = 0;
    for(int k=1;test == 0;k++){
      for(int j=0;j <= m;j++){
        x[(k%3)*(m+1)+j] = dot(m, &K[j*(m+1)], &x[((k-1)%3)*(m+1)]);
      }
      normalize(m, &x[(k%3)*(m+1)]);
      double relerror = 0.0;
      for(int j=0;j <= m;j++){
        double diff = x[(k%3)*(m+1)+j] - x[((k-2)%3)*(m+1)+j];
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
void LM::calc_lm(const bool debug, const int m, const int *cv, const int nalist, const int nblist, const int M,
		 const double *zA, const double *zB,
		 const double crit_move, const double crit_grad, const double crit_dlnL, const double maxsize,
		 double *alnLmax) {

  for(int j=0;j <= m;j++) zavg[j] = 0.0;

  for(int i=1;i <= nalist;i++){
    for(int j=1;j <= m;j++){
      zavg[j] = zavg[j] + zA[i*M + cv[j-1]]; //zA[i][j];
    }
  }

  for(int j=1;j <= m;j++){
    zavg[j] = zavg[j]/((double)nalist);
    alpha0[j] = 0.0;
  }

  /*****  SCREEN FOR STARTING VALUES  *****/


  double mlnl, mlnlold;

  if (debug) printf("\n  GET INITIAL VALUES BY RANDOM NUMBERS");
  for(int i=1;i <= 16;i++){
    a[1] = randomf(-2.0, 2.0);
    a[0] = -a[1]*zavg[m];
    mlnl = mlnL(m, cv, nalist, nblist, M, a, zA, zB);
    if( mlnl < mlnL(m, cv, nalist, nblist, M, alpha0, zA, zB) ){
      for(int j=0;j <= 1;j++){
	alpha0[j] = a[j];
	if (debug) printf("\n improvement %lf a: %.3f %.3f",mlnl, a[0], a[1]);
      }
    }
  }

  alpha0[0] = 0.0;
  alpha0[1] = 0.0;

  if (debug) printf("\n  initial a: ");
  for(int j=0;j <= m;j++){
    a[j] = alpha0[j];
    if (debug) printf(" %.3f ", a[j]);
    for(int k=0;k <= m;k++){
      H[j*(m+1)+k] = 0.0;
    }
    H[j*(m+1)+j] = 1.0 + randomf(0.0, 0.1);
  }

  int test = 0;
  for(int i=1;test==0;i++){
    if(i != 1){
      mlnlold = mlnl;
      for(int j=0;j<=m;j++){
	dmlnLold[j] = dmlnL[j];
      }
    }
    double mlnl = grad(m, cv, nalist, nblist, M, a, zA, zB, dmlnL);
    if (debug) printf("\n dmlnL = %f %f", dmlnL[0], dmlnL[1]);
    if(i != 1){
      BFGS_update(debug, m, x, dmlnL, dmlnLold, H, Hdx_help, y_help, r_help);
    }
    diagonalize(debug, m, H, S, K, K_help, x_help, KK_help, p_help);
    
    // What the heck is this supposed to do?
    //double dlnLsize = dmlnL - dmlnLold;
    double dlnLsize = 20.0;

    double gradsize = sqrt(dot(m, dmlnL, dmlnL));
    if (debug) printf("\n\n  gradients: %d  mlnL: %f  gradsize: %f", i, mlnl, gradsize);
    fflush(stdout);
    for(int j=0;j <= m;j++){
      x[j] = 0.0;
      for(int k=0;k <= m;k++){
	x[j] = x[j] - (dot(m, &S[k*(m+1)], dmlnL)/K[k])*S[k*(m+1)+j];
      }
    }
    double movesize = sqrt(dot(m,x,x));
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

}
