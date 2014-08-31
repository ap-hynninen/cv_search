//
// Roulette-wheel selection.
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
// Uses stochastic acceptance, see:
// Roulette-wheel selection via stochastic acceptance by Adam Lipowski, Dorota Lipowska
// Physica A 391 (2012) pp. 2193-2196
// http://arxiv.org/abs/1109.3627
//
#include <cassert>
#include <iostream>
#ifdef USE_RANDOM
#include <random>
#else
#include <cstdlib>
#endif
#include "RW.hpp"

//
// Class creator
//
RW::RW(const int n_in, const double *w_in) {
  this->n = n_in;
  w = new double[n];
  double w_max = w_in[0];
  for (int i=0;i < n;i++) {
    w[i] = w_in[i];
    w_max = (w[i] > w_max) ? w[i] : w_max;
    assert(w[i] > 0.0);
  }
  inv_w_max = 1.0/w_max;
}

//
// Class destructor
//
RW::~RW() {
  delete [] w;
}

//
// Pick a random individual
//
#ifdef USE_RANDOM
int RW::pick(std::mt19937 &rand_eng) {
  std::uniform_int_distribution<> pick_indv(0, n-1);
  std::uniform_real_distribution<> pick_prob(0.0, 1.0);
  int i;
  while (true) {
    i = pick_indv(rand_eng);
    if (pick_prob(rand_eng) < w[i]*inv_w_max) break;
  }
  return i;
}
#else
int RW::pick() {
  int i;
  while (true) {
    i = rand() % n;
    if ( ((double)rand())/((double)RAND_MAX) < w[i]*inv_w_max) break;
  }
  return i;
}
#endif
