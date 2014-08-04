
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// hA[N], hB[N]
// alist[N], blist[N]
//
void AsBs(const int N, const int *hA, const int *hB, int &nalist, int *alist, int &nblist, int *blist) {
  int i, j, check, postime, negtime, flagshort, flagoverlap;
  nalist = 0;
  nblist = 0;
  for(i=0;i < N;i++){
    /* diagnostic: overlapping basins */
    flagoverlap = hA[i]*hB[i];
    if( flagoverlap == 0 ){
      if( hA[i] == 1 ){
	alist[nalist++] = i;
      }
      if( hB[i] == 1 ){
	blist[nblist++] = i;
      }
    }
  }
  printf("\n  n(A)=%d  n(B)=%d ", nalist, nblist);
}

double getmin(const int N, const double *q) {
  int i;
  double temp;
  temp = q[0];
  for(i=1;i < N;i++){
    if(q[i] < temp){
      temp = q[i];
    }
  }
  return temp;
}

double getmax(const int N, const double *q) {
  int i;
  double temp;
  temp = q[0];
  for(i=1;i < N;i++){
    if(q[i] > temp){
      temp = q[i];
    }
  }
  return temp;
}

//
//
// q[M][N]
// hA[N]
// hB[N]
//
void read_data(const char *filename, const int N, const int M, double *q, int *hA, int *hB) {

  FILE *data;
  int i, j, k, bad_apples, check, hbA, hbB;
  double temp;
  
  data = fopen(filename,"r");
  bad_apples = 0;
  for(i=0;i < N;i++){
    /* shooting point number */
    fscanf(data, "%d ", &k);
    if(k != i+1){
      printf("\n  ERROR READING DATA FROM %s\n",filename);
      printf("i=%d\n",i);
      exit(1);
    }
    /* backward to A ? */
    fscanf(data, "%d ", &hbA);
    /* backward to B ? */
    fscanf(data, "%d ", &hbB);
    /* forward to A ? */
    fscanf(data, "%d ", &hA[i]);
    /* forward to B ? */
    fscanf(data, "%d ", &hB[i]);
    check = hbA*hB[i] + hA[i]*hbB;
    if(check == 2){
      printf("\n\n  your h-values are wrong or your basins overlap! \n");
      printf("i=%d\n",i);
      exit(1);
    }
    bad_apples = bad_apples + check;
    for(j=0;j < M;j++){
      fscanf(data, "%lf ", &temp);
      q[j*N + i] = temp;
    }
  }
  printf("\n\n  accepted %d of %d\n", bad_apples, N);
}
