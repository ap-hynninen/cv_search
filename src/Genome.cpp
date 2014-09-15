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

Genome::Genome() : num_cv(0), gene(NULL) {
}

//
// Class destructor
//
Genome::~Genome() {
  if (gene != NULL) delete [] gene;
}

//
// Sets the number of genes in this genome
//
void Genome::set_num_cv(const int num_cv) {
  this->num_cv = num_cv;
  gene = new Gene[num_cv];
}
