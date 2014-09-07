//
// Genome class
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include <iostream>
#include "../include/Genome.hpp"

//
// Class creators
//
Genome::Genome(const int num_cv) {
  set_num_cv(num_cv);
}

Genome::Genome() : num_cv(0), pair(NULL) {
}

//
// Class destructor
//
Genome::~Genome() {
  if (pair != NULL) delete [] pair;
}

//
// Sets the number of pairs in this genome
//
void Genome::set_num_cv(const int num_cv) {
  this->num_cv = num_cv;
  pair = new Pair[num_cv];
}
