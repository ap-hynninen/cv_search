#ifndef RW_HPP
#define RW_HPP
//
// Roulette-wheel selection.
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
// Uses stochastic acceptance, see:
// Roulette-wheel selection via stochastic acceptance by Adam Lipowski, Dorota Lipowska
// Physica A 391 (2012) pp. 2193-2196
// http://arxiv.org/abs/1109.3627
//
#ifdef USE_RANDOM
#include <random>
#endif

class RW {

private:

  // Number of individuals
  int n;

  // Weights
  double *w;

  // Inverse of the maximum weight
  double inv_w_max;

public:

  RW(const int n_in, const double *w_in);
  ~RW();

#ifdef USE_RANDOM
  int pick(std::mt19937 &rand_eng);
#else
  int pick();
#endif

};

#endif // RW_HPP
