
void AsBs(const int N, const int *hA, const int *hB, int &nalist, int *alist, int &nblist, int *blist);
double getmin(const int N, const double *q);
double getmax(const int N, const double *q);
void read_data(const char *filename, const int N, const int M, double *q, int *hA, int *hB);
void calc_lm(const bool debug, const int m, const int *cv, const int nalist, const int nblist, const int M,
	     const double *zA, const double *zB,
	     const double crit_move, const double crit_grad, const double crit_dlnL, const double maxsize,
	     double *alnLmax);
