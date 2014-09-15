#include <cassert>
#include <iostream>
#include <cmath>
#include "../include/Gene.hpp"

//
// Returns true if atom list has duplicate atoms
//
bool Gene::has_duplicate() {
  if (atom.size() == 2) {
    return (atom[0] == atom[1]);
  } else if (atom.size() == 3) {
    return (atom[0] == atom[1]) || (atom[0] == atom[2]) || (atom[1] == atom[2]);
  } else {
    std::cerr << "Gene::has_duplicate, size > 3 not supported " << std::endl;
    exit(1);
  }
}

//
// Evaluate
//
float Gene::eval(const float3 *coord) {

  if (atom.size() == 2) {
    int i = atom[0];
    int j = atom[1];
    float dx = coord[i].x - coord[j].x;
    float dy = coord[i].y - coord[j].y;
    float dz = coord[i].z - coord[j].z;
    return sqrt(dx*dx + dy*dy + dz*dz);
  } else if (atom.size() == 3) {
    int i = atom[0];
    int j = atom[1];
    int k = atom[2];
    float x0 = coord[i].x - coord[j].x;
    float y0 = coord[i].y - coord[j].y;
    float z0 = coord[i].z - coord[j].z;
    float x2 = coord[k].x - coord[j].x;
    float y2 = coord[k].y - coord[j].y;
    float z2 = coord[k].z - coord[j].z;
    float x0len = sqrtf(x0*x0 + y0*y0 + z0*z0);
    float x2len = sqrtf(x2*x2 + y2*y2 + z2*z2);
    if (x0len < 1.0e-14 || x2len < 1.0e-14) return 0.0f;
    return acosf((x0*x2 + y0*y2 + z0*z2)/(x0len*x2len));
  } else {
    std::cerr << "Gene::eval, size > 3 not supported " << std::endl;
    exit(1);
  }
}

//
// Prints out the atoms involved in this gene
//
void Gene::print() {
  for (int i=0;i < atom.size();i++) {
    std::cout << "[ " << atom[i] + 1 << " ] ";
  }
  std::cout << std::endl;
}
