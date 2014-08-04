
#include <stdio.h>
#include "cv_util.h"
#include "LM.hpp"

int main() {

  // Number of shooting points
  const int N = 10000;

  // Number of CVs
  const int M = 149;
  
  // q[M][N]
  double* q = new double[M*N];
  int* hA = new int[N];
  int* hB = new int[N];

  read_data("data.txt", N, M, q, hA, hB);

  printf("\n  SCANNED INPUT FROM data.txt");

  int* alist = new int[N];
  int* blist = new int[N];
  int nalist, nblist;
  AsBs(N, hA, hB, nalist, alist, nblist, blist);

  printf("\n  SORTED DATA BY A AND B OUTCOMES");

  // z[M][N]
  double* z = new double[M*N];
  double* qmax = new double[M];
  double* qmin = new double[M];
  double* qspan = new double[M];
  for(int i=0;i < M;i++){
    qmax[i] = getmax(N, &q[i*N]);
    qmin[i] = getmin(N, &q[i*N]);
    qspan[i] = qmax[i] - qmin[i];
    printf("\n %d %.3f %.3f %.3f ", i+1, qmin[i], qmax[i], qspan[i]); 
  }

  printf("\n  REDUCED VARIABLES: Z = (Q-Qmin)/(Qmax-Qmin) \n");
  printf("\n  SAMPLED Z's IN [0, 1], Q = Z(Qmax-Qmin)+Qmin\n"); 
  for(int i=0;i < M;i++){
    for(int j=0;j < N;j++){
      z[i*N + j] = (q[i*N + j]-qmin[i])/(qmax[i]-qmin[i]);
    }
  }

  delete [] q;
  delete [] qmax;
  delete [] qmin;
  delete [] qspan;
  delete [] hA;
  delete [] hB;

  // z[i][alist[...]] belong to alist

  // zA[nalist][M]
  // zB[nblist][M]
  double *zA = new double[nalist*M];
  double *zB = new double[nblist*M];

  for(int i=0;i < M;i++){
    for(int j=0;j < nalist;j++){
      zA[j*M + i] = z[i*N + alist[j]];
    }
    for(int j=0;j < nblist;j++){
      zB[j*M + i] = z[i*N + blist[j]];
    }
  }

  delete [] z;
  delete [] alist;
  delete [] blist;

  //
#define m_max 3
  int m;
  int cv[m_max];
  double crit_move = 0.001;
  double crit_grad = 0.01;
  double crit_dlnL = 0.0001;
  double maxsize = 0.1;
  double alnLmax[m_max+2];
  double max_val[m_max+2];
  int max_ind[m_max];
  LM *lm;

  m = 1;
  lm = new LM(m);
  max_val[m+1] = -1.0e100;
  max_ind[0] = -1;
  for (int i=0;i < M;i++) {
    cv[0] = i;
    lm->calc_lm(false, m, cv, nalist, nblist, M, zA, zB, crit_move, crit_grad, crit_dlnL, maxsize, alnLmax);
    printf("%d rxncoor",i+1);
    for (int jj=0;jj < m+2;jj++) printf(" %f",alnLmax[jj]);
    printf("\n");
    if (alnLmax[m+1] > max_val[m+1]) {
      for (int jj=0;jj < m+2;jj++) max_val[jj] = alnLmax[jj];
      max_ind[0] = i;
    }
  }
  if (max_ind[0] == -1) {
    printf("ERROR\n");
    return 0;
  }
  printf("best rxncoor %d",max_ind[0]+1);
  for (int jj=0;jj < m+2;jj++) printf(" %f",max_val[jj]);
  printf("\n");
  delete lm;

  m = 2;
  lm = new LM(m);
  max_val[m+1] = -1.0e100;
  max_ind[0] = -1;
  max_ind[1] = -1;
  for (int i=0;i < M;i++) {
    cv[0] = i;
    for (int j=0;j < M;j++) {
      cv[1] = j;
      lm->calc_lm(false, m, cv, nalist, nblist, M, zA, zB, crit_move, crit_grad, crit_dlnL, maxsize, alnLmax);
      printf("%d %d rxncoor",i+1,j+1);
      for (int jj=0;jj < m+2;jj++) printf(" %f",alnLmax[jj]);
      printf("\n");
      if (alnLmax[m+1] > max_val[m+1]) {
	for (int jj=0;jj < m+2;jj++) max_val[jj] = alnLmax[jj];
	max_ind[0] = i;
	max_ind[1] = j;
      }
    }
  }
  if (max_ind[0] == -1) {
    printf("ERROR\n");
    return 0;
  }
  printf("best rxncoor %d %d",max_ind[0]+1,max_ind[1]+1);
  for (int jj=0;jj < m+2;jj++) printf(" %f",max_val[jj]);
  printf("\n");
  delete lm;

  m = 3;
  lm = new LM(m);
  max_val[m+1] = -1.0e100;
  max_ind[0] = -1;
  max_ind[1] = -1;
  max_ind[2] = -1;
  for (int i=0;i < M;i++) {
    cv[0] = i;
    for (int j=0;j < M;j++) {
      cv[1] = j;
      for (int k=0;k < M;k++) {
	cv[2] = k;
	lm->calc_lm(false, m, cv, nalist, nblist, M, zA, zB, crit_move, crit_grad, crit_dlnL, maxsize, alnLmax);
	printf("%d %d %d rxncoor",i+1,j+1,k+1);
	for (int jj=0;jj < m+2;jj++) printf(" %f",alnLmax[jj]);
	printf("\n");
	if (alnLmax[m+1] > max_val[m+1]) {
	  for (int jj=0;jj < m+2;jj++) max_val[jj] = alnLmax[jj];
	  max_ind[0] = i;
	  max_ind[1] = j;
	  max_ind[2] = k;
	}
      }
    }
  }
  if (max_ind[0] == -1) {
    printf("ERROR\n");
    return 0;
  }
  printf("best rxncoor %d %d %d",max_ind[0]+1,max_ind[1]+1,max_ind[1]+1);
  for (int jj=0;jj < m+2;jj++) printf(" %f",max_val[jj]);
  printf("\n");
  delete lm;

  delete [] zA;
  delete [] zB;
  
  return 1;
}
