#ifndef GENE_HPP
#define GENE_HPP
//
// Cluster gene
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include <iostream>
#include <vector>
#include "struct.h"

class Gene {

private:

  // Atom indices
  std::vector<int> atom;

public:
  Gene() {}
  ~Gene() {}

  void set(const std::vector<int> *cluster) {atom = *cluster;}
  void add(const int i) {atom.push_back(i);}
  void replace(const int j, const int i) {
    if (j >= atom.size()) {
      std::cerr << "Gene::replace, j out of bounds" << std::endl;
      exit(1);
    }
    atom[j] = i;
  }
  std::vector<int>* get_atom() {return &atom;}
  int get_natom() {return atom.size();}

  bool has_duplicate();

  float eval(const float3 *coord);

  void print();
};

#endif // GENE_HPP
