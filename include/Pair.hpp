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

  // Atom indices for each pair atom1[0...npair-1]
  std::vector<int> atom1;
  std::vector<int> atom2;

  void calc_com(const float3 *coord, const float *mass, const int natom, const int* atom,
		float &x, float &y, float &z);

public:
  Pair();
  ~Pair();

  void set(const int atom1_start, const int atom1_stop,
	   const int atom2_start, const int atom2_stop);

  void set(const std::vector<int> *cluster1,
	   const std::vector<int> *cluster2);

  void add1(const int i);
  void add2(const int i);

  void flip_atom1(const int iatom);
  void flip_atom2(const int iatom);

  bool contains1(const int iatom);
  bool contains2(const int iatom);

  std::vector<int>* get_atom1() {return &atom1;}
  std::vector<int>* get_atom2() {return &atom2;}

  int get_natom1() {return atom1.size();}
  int get_natom2() {return atom2.size();}

  float eval(const float3 *coord, const float *mass);

  void print();
};

#endif // PAIR_HPP
