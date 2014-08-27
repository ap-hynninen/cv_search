#ifndef COORD_HPP
#define COORD_HPP
//
// Coordinate class
// Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include "struct.h"

class Coord {

private:
  
  // Number of coordinates
  int ncoord;

  // Number of shooting points
  int nshoot;

  // Coordinates: coord[0...nshoot-1][0...ncoord-1]
  float3 **coord;

  // Masses: mass[0...ncoord-1]
  float *mass;

  // Residue ids
  int *resid;

  // Number of residues
  int nresidue;

  // Maximum residue size
  int max_residue_size;

  // Residue starts
  int *residue_start;

  void load_coord(const char *filename);
  void setup_residues();

public:

  Coord(const char *coord_filename, const int ncoord, const int nshoot);
  ~Coord();

  int get_ncoord() {return ncoord;}
  int get_nresidue() {return nresidue;}
  int get_residue_start(int i) {return residue_start[i];}
  int get_residue_stop(int i) {return residue_start[i+1]-1;}
  int get_max_residue_size() {return max_residue_size;}
  int get_nshoot() {return nshoot;}
  float3 *get_coord(int ishoot) {return coord[ishoot];}
  float *get_mass() {return mass;}

};

#endif // COORD_HPP
