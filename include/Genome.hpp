#ifndef GENOME_HPP
#define GENOME_HPP
//
// Genome class
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include "Gene.hpp"

class Genome {

private:

  // Number of genes in this genome (typically num_cv <= 3 or so)
  int num_cv;

  // Genes in this genome (gene[0..num_cv-1])
  Gene *gene;

public:

  Genome(const int num_cv);
  Genome();
  ~Genome();

  void set_num_cv(const int num_cv);
  Gene* get_gene(const int igene) {return &gene[igene];}

};

#endif // GENOME_HPP
