#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "Coord.hpp"

//
// Class constructor
//
Coord::Coord(const char *coord_filename) {

  coord = NULL;
  resid = NULL;

  load_coord(coord_filename);

}

//
// Class destructor
//
Coord::~Coord() {

  if (coord != NULL) {
    for (int i=0;i < nshoot;i++)
      delete [] coord[i];
    delete [] coord;
  }

  if (resid != NULL) delete [] resid;
}

//
// Load coordinates
//
void Coord::load_coord(const char *filename) {

  std::ifstream file;
  file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  try {
    // Open file
    file.open(filename);

    // Read number of coordinates and number of shooting points
    file >> ncoord >> nshoot;

    // Allocate coord
    coord = new float3*[nshoot];
    for (int i=0;i < nshoot;i++)
      coord[i] = new float3[ncoord];

    // Allocate resid
    resid = new int[ncoord];

    // Read coordinates
    for (int k=0;k < nshoot;k++) {
      for (int i=0;i < ncoord;i++) {
	char atomstr[6], atomname[6], resname[8];
	int atomid, resid_tmp;
	file >> atomstr >> atomid >> atomname >> resname >> resid_tmp
	     >> coord[k][i].x >> coord[k][i].y >> coord[k][i].z;
	if (k == 0) resid[i] = resid_tmp;
      }
    }

    // Close file
    file.close();
  }
  catch(std::ifstream::failure e) {
    std::cerr << "Error opening/reading/closing file " << filename << std::endl;
    exit(1);
  }

}

