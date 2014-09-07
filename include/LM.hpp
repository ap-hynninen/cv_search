#ifndef LM_HPP
#define LM_HPP
//
// Likelyhood Maximization class
//
// Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
// Based on code by Baron Peters
//

template <typename real>
class LM {

private:
  // Private variables
  real *zavg;
  real *alpha0;
  real *a;
  real *aold;
  real *H;
  real *dmlnL;
  real *x;
  real *dmlnLold;
  real *Hdx_help;
  real *y_help;
  real *r_help;
  real *S;
  real *K;
  real *K_help;
  real *KK_help;
  real *x_help;
  real *p_help;

  // Private subroutines
  real randomf(real a, real b);
  real rxn_coor(const int m, const int *cv, const real *a, const real *z);
  real mlnL(const int m, const int *cv, const int nalist, const int nblist, const int M,
	      const real *a, const real *zA, const real *zB);
  real grad(const int m, const int *cv, const int nalist, const int nblist, const int M,
	      const real *a, const real *zA, const real *zB, real *dmlnL);
  real dot(const int m, const real *v, const real *u);
  void BFGS_update(const bool debug, const int m, const real *dx, const real *g_new, const real *g_old,
		   real *h, real *Hdx, real *y, real *r);
  real normalize(const int m, real *x);
  void projectfrommatrix(const int m, const real *vector, real *hmwc,
			 real *p, real *K);
  void diagonalize(const bool debug, const int m, const real *hmwc, real *smwc, real *w2,
		   real *K, real *x, real *KK, real *p);

public:
  LM(const int m);
  ~LM();

  void calc_lm(const bool debug, const int m, const int *cv, const int nalist, const int nblist, const int M,
	       const real *zA, const real *zB,
	       const real crit_move, const real crit_grad, const real crit_dlnL, const real maxsize,
	       real *alnLmax);

};

#endif // LM_HPP
