//
// Likelyhood Maximization class
//
// Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
// Based on code by Baron Peters
//

class LM {

private:
  // Private variables
  double *zavg;
  double *alpha0;
  double *a;
  double *aold;
  double *H;
  double *dmlnL;
  double *x;
  double *dmlnLold;
  double *Hdx_help;
  double *y_help;
  double *r_help;
  double *S;
  double *K;
  double *K_help;
  double *KK_help;
  double *x_help;
  double *p_help;

  // Private subroutines
  double randomf(double a, double b);
  double rxn_coor(const int m, const int *cv, const double *a, const double *z);
  double mlnL(const int m, const int *cv, const int nalist, const int nblist, const int M,
	      const double *a, const double *zA, const double *zB);
  double grad(const int m, const int *cv, const int nalist, const int nblist, const int M,
	      const double *a, const double *zA, const double *zB, double *dmlnL);
  double dot(const int m, const double *v, const double *u);
  void BFGS_update(const bool debug, const int m, const double *dx, const double *g_new, const double *g_old,
		   double *h, double *Hdx, double *y, double *r);
  double normalize(const int m, double *x);
  void projectfrommatrix(const int m, const double *vector, double *hmwc,
			 double *p, double *K);
  void diagonalize(const bool debug, const int m, const double *hmwc, double *smwc, double *w2,
		   double *K, double *x, double *KK, double *p);

public:
  LM(const int m);
  ~LM();

  void calc_lm(const bool debug, const int m, const int *cv, const int nalist, const int nblist, const int M,
	       const double *zA, const double *zB,
	       const double crit_move, const double crit_grad, const double crit_dlnL, const double maxsize,
	       double *alnLmax);

};
