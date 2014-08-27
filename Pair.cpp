#include <cassert>
#include <cmath>
#include "Pair.hpp"

//
// Class constructor
//
Pair::Pair(const int npair) {
  this->npair = npair;

  atom1 = new std::vector<int>[npair];
  atom2 = new std::vector<int>[npair];

}

//
// Class destructor
//
Pair::~Pair() {
  delete [] atom1;
  delete [] atom2;
}

//
// Sets cluster atoms from index range
//
void Pair::set(const int ipair,
	       const int atom1_start, const int atom1_stop,
	       const int atom2_start, const int atom2_stop) {
  assert(ipair >= 0);
  assert(ipair < npair);
  assert(atom1_stop >= atom1_start);
  assert(atom2_stop >= atom2_start);

  atom1[ipair].resize(atom1_stop - atom1_start + 1);
  atom2[ipair].resize(atom2_stop - atom2_start + 1);

  for (int i=0;i < (atom1_stop - atom1_start + 1);i++) {
    atom1[ipair][i] = i + atom1_start;
  }

  for (int i=0;i < (atom2_stop - atom2_start + 1);i++) {
    atom2[ipair][i] = i + atom2_start;
  }

}

//
// Sets cluster atoms from two other clusters
//
void Pair::set(const int ipair,
	       const std::vector<int> *cluster1,
	       const std::vector<int> *cluster2) {
  assert(ipair >= 0);
  assert(ipair < npair);

  atom1[ipair] = *cluster1;
  atom2[ipair] = *cluster2;
}

//
// Adds atom i
//
void Pair::add1(const int ipair, const int i) {
  atom1[ipair].push_back(i);
}

//
// Adds atom i
//
void Pair::add2(const int ipair, const int i) {
  atom2[ipair].push_back(i);
}

//
// Flip atom1
//
void Pair::flip_atom1(const int ipair, const int iatom) {
  
  for (int i=0;i < atom1[ipair].size();i++) {
    if (atom1[ipair][i] == iatom) {
      atom1[ipair].erase(atom1[ipair].begin()+i);
      return;
    }
  }

  atom1[ipair].push_back(iatom);

}

//
// Flip atom2
//
void Pair::flip_atom2(const int ipair, const int iatom) {
  
  for (int i=0;i < atom2[ipair].size();i++) {
    if (atom2[ipair][i] == iatom) {
      atom2[ipair].erase(atom2[ipair].begin()+i);
      return;
    }
  }

  atom2[ipair].push_back(iatom);

}

//
// Returns true if atom1 contains iatom
//
bool Pair::contains1(const int ipair, const int iatom) {

  for (int i=0;i < atom1[ipair].size();i++) {
    if (atom1[ipair][i] == iatom) {
      return true;
    }
  }

  return false;
}

//
// Returns true if atom2 contains iatom
//
bool Pair::contains2(const int ipair, const int iatom) {

  for (int i=0;i < atom2[ipair].size();i++) {
    if (atom2[ipair][i] == iatom) {
      return true;
    }
  }

  return false;
}

//
// Evaluate
//
void Pair::eval(const float3 *coord, const float *mass, const int N, float *r) {

  for (int ipair=0;ipair < npair;ipair++) {

    // Calculate center-of-mass
    float x1, y1, z1;
    calc_com(coord, mass, atom1[ipair].size(), atom1[ipair].data(), x1, y1, z1);

    float x2, y2, z2;
    calc_com(coord, mass, atom2[ipair].size(), atom2[ipair].data(), x2, y2, z2);

    float dx = x1-x2;
    float dy = y1-y2;
    float dz = z1-z2;

    r[N*ipair] = sqrt(dx*dx + dy*dy + dz*dz);
  }

}

//
// Calculates center-of-mass
//
void Pair::calc_com(const float3 *coord, const float *mass, const int natom, const int* atom,
		    float &x, float &y, float &z) {
  x = 0.0f;
  y = 0.0f;
  z = 0.0f;
  float masstot = 0.0f;
  for (int i=0;i < natom;i++) {
    int ii = atom[i];
    float massi = mass[i];
    x += coord[ii].x*massi;
    y += coord[ii].y*massi;
    z += coord[ii].z*massi;
    masstot += massi;
  }
  float inv_masstot = 1.0f/masstot;
  x *= inv_masstot;
  y *= inv_masstot;
  z *= inv_masstot;
}
