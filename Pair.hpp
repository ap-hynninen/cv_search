#ifndef PAIR_HPP
#define PAIR_HPP
//
// Cluster pair
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include <vector>
#include "struct.h"

class Pair {

private:
  int npair;

  // Atom indices for each pair atom1[0...npair-1]
  std::vector<int> *atom1;
  std::vector<int> *atom2;

  void calc_com(const float3 *coord, const float *mass, const int natom, const int* atom,
		float &x, float &y, float &z);

public:
  Pair(const int npair);
  ~Pair();

  void set(const int ipair,
	   const int atom1_start, const int atom1_stop,
	   const int atom2_start, const int atom2_stop);

  void set(const int ipair,
	   const std::vector<int> *cluster1,
	   const std::vector<int> *cluster2);

  void add1(const int ipair, const int i);
  void add2(const int ipair, const int i);

  void flip_atom1(const int ipair, const int iatom);
  void flip_atom2(const int ipair, const int iatom);

  bool contains1(const int ipair, const int iatom);
  bool contains2(const int ipair, const int iatom);

  std::vector<int>* get_atom1(const int ipair) {return &atom1[ipair];}
  std::vector<int>* get_atom2(const int ipair) {return &atom2[ipair];}

  int get_natom1(const int ipair) {return atom1[ipair].size();}
  int get_natom2(const int ipair) {return atom2[ipair].size();}

  void eval(const float3 *coord, const float *mass, const int N, float *r);
};

#endif // PAIR_HPP
