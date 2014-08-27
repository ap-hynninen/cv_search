#ifndef GA_HPP
#define GA_HPP
//
// Genetic Algorithm for collective variable (CV) finding
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include <random>
#include "Coord.hpp"
#include "Pair.hpp"

struct lnval_t {
  int ind;
  double data[3];
};

class GA {

private:

  // Population sizes
  int npair;

  // Coordinates
  Coord *coord;

  int nalist;
  int* alist;

  int nblist;
  int* blist;

  double *zA;
  double *zB;

  // Pairs
  Pair *pair;
  float *pair_cv;

  // Random number engine
  std::mt19937 rand_eng;

  void eval_cv_values();
  //  void build_next_generation_old(const int ntop_pair, const int *top_pair_ind, Pair *top_pair);
  void build_next_generation(std::vector<lnval_t> &lnval);

public:

  GA(Coord *coord, const int npair, const int *hA, const int *hB);
  ~GA();

  void init_population();
  void run(const int niter);

};

#endif // GA_HPP
