
#include <omp.h>
#include <stdio.h>
#include "cv_util.h"
#include "LM.hpp"

int LM_from_data();
int LM_from_coord(char *filename);

int main(int argc, char **argv) {

  if (argc == 1) {
    LM_from_data();
  } else {
    LM_from_coord(argv[1]);
  }

  return 1;
}

//
//
//
int LM_from_coord(char *filename) {
  printf("LM_from_coord = %s\n",filename);

  return 1;
}

//
// Does likelihood maximization from data.txt
//
int LM_from_data() {
  // Number of shooting points
  const int N = 20000;

  // Number of CVs
  const int M = 149;

  // Number of CV combinations
  const int m = 1;

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

#define MAX_NTHREAD 24
  int nthread;
#pragma omp parallel
  {
    if (omp_get_thread_num() == 0) nthread = omp_get_num_threads();
  }
  printf("Using %d threads\n",nthread);

  //
#define m_max 3
  double crit_move = 0.001;
  double crit_grad = 0.01;
  double crit_dlnL = 0.0001;
  double maxsize = 0.1;
  double max_val[MAX_NTHREAD][m_max+2];
  int max_ind[MAX_NTHREAD][m_max];

  if (m == 1) {
#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      double alnLmax[m_max+2];
      //LM *lm = new LM(m);
      LM lm(m);
      max_val[tid][m+1] = -1.0e100;
      max_ind[tid][0] = -1;
      int i;
#pragma omp for private(i) schedule(dynamic)
      for (i=0;i < M;i++) {
	int cv[m_max];
	cv[0] = i;
	lm.calc_lm(false, m, cv, nalist, nblist, M, zA, zB, crit_move, crit_grad, crit_dlnL, maxsize, alnLmax);
	printf("%d rxncoor",i+1);
	for (int jj=0;jj < m+2;jj++) printf(" %f",alnLmax[jj]);
	printf("\n");
	if (alnLmax[m+1] > max_val[tid][m+1]) {
	  for (int jj=0;jj < m+2;jj++) max_val[tid][jj] = alnLmax[jj];
	  max_ind[tid][0] = i;
	}
      }
      fprintf(stderr,"Hello!\n");
      //delete lm;
    }

    fprintf(stderr,"Hello\n");

    for (int i=1;i < nthread;i++) {
      if (max_val[i][m+1] > max_val[0][m+1]) {
	for (int jj=0;jj < m+2;jj++) max_val[0][jj] = max_val[i][jj];
	max_ind[0][0] = max_ind[i][0];
      }
    }
    
    if (max_ind[0][0] == -1) {
      printf("ERROR\n");
      return 0;
    }
    printf("best rxncoor %d",max_ind[0][0]+1);
    for (int jj=0;jj < m+2;jj++) printf(" %f",max_val[0][jj]);
    printf("\n");
  } else if (m == 2) {
#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      double alnLmax[m_max+2];
      LM *lm = new LM(m);
      max_val[tid][m+1] = -1.0e100;
      max_ind[tid][0] = -1;
      max_ind[tid][1] = -1;
      int i;
#pragma omp for private(i) schedule(dynamic)
      for (i=0;i < M;i++) {
	int cv[m_max];
	cv[0] = i;
	for (int j=0;j < M;j++) {
	  cv[1] = j;
	  lm->calc_lm(false, m, cv, nalist, nblist, M, zA, zB, crit_move, crit_grad, crit_dlnL, maxsize, alnLmax);
	  printf("%d %d rxncoor",i+1,j+1);
	  for (int jj=0;jj < m+2;jj++) printf(" %f",alnLmax[jj]);
	  printf("\n");
	  if (alnLmax[m+1] > max_val[tid][m+1]) {
	    for (int jj=0;jj < m+2;jj++) max_val[tid][jj] = alnLmax[jj];
	    max_ind[tid][0] = i;
	    max_ind[tid][1] = j;
	  }
	}
      }
      delete lm;
    }
    
    for (int i=1;i < nthread;i++) {
      if (max_val[i][m+1] > max_val[0][m+1]) {
	for (int jj=0;jj < m+2;jj++) max_val[0][jj] = max_val[i][jj];
	max_ind[0][0] = max_ind[i][0];
	max_ind[0][1] = max_ind[i][1];
      }
    }
    
    if (max_ind[0][0] == -1) {
      printf("ERROR\n");
      return 0;
    }
    printf("best rxncoor %d %d",max_ind[0][0]+1,max_ind[0][1]+1);
    for (int jj=0;jj < m+2;jj++) printf(" %f",max_val[0][jj]);
    printf("\n");
  } else if (m == 3) {
#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      double alnLmax[m_max+2];
      LM *lm = new LM(m);
      max_val[tid][m+1] = -1.0e100;
      max_ind[tid][0] = -1;
      max_ind[tid][1] = -1;
      max_ind[tid][2] = -1;
      int i;
#pragma omp for private(i) schedule(dynamic)
      for (i=0;i < M;i++) {
	int cv[m_max];
	cv[0] = i;
	for (int j=0;j < M;j++) {
	  cv[1] = j;
	  for (int k=0;k < M;k++) {
	    cv[2] = k;
	    lm->calc_lm(false, m, cv, nalist, nblist, M, zA, zB, crit_move, crit_grad, crit_dlnL, maxsize, alnLmax);
	    printf("%d %d %d rxncoor",i+1,j+1,k+1);
	    for (int jj=0;jj < m+2;jj++) printf(" %f",alnLmax[jj]);
	    printf("\n");
	    if (alnLmax[m+1] > max_val[tid][m+1]) {
	      for (int jj=0;jj < m+2;jj++) max_val[tid][jj] = alnLmax[jj];
	      max_ind[tid][0] = i;
	      max_ind[tid][1] = j;
	      max_ind[tid][1] = k;
	    }
	  }
	}
      }
      delete lm;
    }
    
    for (int i=1;i < nthread;i++) {
      if (max_val[i][m+1] > max_val[0][m+1]) {
	for (int jj=0;jj < m+2;jj++) max_val[0][jj] = max_val[i][jj];
	max_ind[0][0] = max_ind[i][0];
	max_ind[0][1] = max_ind[i][1];
	max_ind[0][1] = max_ind[i][2];
      }
    }
    
    if (max_ind[0][0] == -1) {
      printf("ERROR\n");
      return 0;
    }
    printf("best rxncoor %d %d %d",max_ind[0][0]+1,max_ind[0][1]+1,max_ind[0][2]+1);
    for (int jj=0;jj < m+2;jj++) printf(" %f",max_val[0][jj]);
    printf("\n");
  } else {
    printf("Incorrect m value\n");
  }

  /*
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
  */

  delete [] zA;
  delete [] zB;
 
  return 1;
}
