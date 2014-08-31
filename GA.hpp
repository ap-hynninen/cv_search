#ifndef GA_HPP
#define GA_HPP
//
// Genetic Algorithm for collective variable (CV) finding
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#ifdef USE_RANDOM
#include <random>
#endif
#include "Coord.hpp"
#include "Genome.hpp"

#define max_num_cv 4

struct lnval_t {
  int ind;
  double data[max_num_cv+2];
  double key;
};

class GA {

private:

  // Number of genomes = population size
  int ngenome;

  // Number of collective variables (CV) in each genome (typically 1...4 or so)
  int num_cv;

  // Genome
  Genome *genome;

  // Genome CV values
  float *genome_cv;

  // Coordinates
  Coord *coord;

  int nalist;
  int* alist;

  int nblist;
  int* blist;

  double *zA;
  double *zB;

  // Random number engine
#ifdef USE_RANDOM
  std::mt19937 rand_eng;
#endif

  void eval_cv_values();
  //  void build_next_generation_old(const int ntop_pair, const int *top_pair_ind, Pair *top_pair);
  void build_next_generation(std::vector<lnval_t> &lnval);

public:

  GA(Coord *coord, const int num_cv, const int ngenome, const int *hA, const int *hB);
  ~GA();

  void init_population();
  void run(const int niter);

};

#endif // GA_HPP
