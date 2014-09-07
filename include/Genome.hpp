#ifndef GENOME_HPP
#define GENOME_HPP
//
// Genome class
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include "Pair.hpp"

class Genome {

private:

  // Number of pairs in this genome (typically num_cv <= 3 or so)
  int num_cv;

  // Pairs in this genome (pair[0..num_cv-1])
  Pair *pair;

public:

  Genome(const int num_cv);
  Genome();
  ~Genome();

  void set_num_cv(const int num_cv);
  Pair* get_pair(const int ipair) {return &pair[ipair];}

};

#endif // GENOME_HPP
