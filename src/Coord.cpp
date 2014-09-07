#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "../include/Coord.hpp"

//
// Class constructor
//
Coord::Coord(const char *coord_filename, const int ncoord, const int nshoot) {

  this->ncoord = ncoord;
  this->nshoot = nshoot;

  coord = NULL;
  resid = NULL;
  mass = NULL;
  residue_start = NULL;

  load_coord(coord_filename);

  setup_residues();
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
  if (mass != NULL) delete [] mass;
  if (residue_start != NULL) delete [] residue_start;
}

//
// Setups residues
//
void Coord::setup_residues() {

  // Count the number of residues
  nresidue = 0;
  int prev_resid = -1;
  for (int i=0;i < ncoord;i++) {
    if (resid[i] != prev_resid) {
      prev_resid = resid[i];
      nresidue++;
    }
  }

  residue_start = new int[nresidue+1];

  // Mark residue starts
  prev_resid = -1;
  int j = 0;
  for (int i=0;i < ncoord;i++) {
    if (resid[i] != prev_resid) {
      prev_resid = resid[i];
      residue_start[j] = i;
      j++;
    }
  }
  residue_start[nresidue] = ncoord;

  // Calculate the maximum residue size
  max_residue_size = 0;
  for (int i=0;i < nresidue;i++) {
    //std::cout << i << " " << get_residue_start(i) << " " << get_residue_stop(i) << std::endl;
    int residue_size = get_residue_stop(i) - get_residue_start(i) + 1;
    max_residue_size = (max_residue_size < residue_size) ? residue_size : max_residue_size;
  }

  std::cout << "Coord::setup_residues, nresidue = " << nresidue
	    << " max_residue_size = " << max_residue_size << std::endl;

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

    std::cout << "Coord::load_coord, ncoord = " << ncoord << " nshoot = " << nshoot << std::endl;

    // Allocate coord
    coord = new float3*[nshoot];
    for (int i=0;i < nshoot;i++)
      coord[i] = new float3[ncoord];

    // Allocate resid
    resid = new int[ncoord];

    // Allocate mass
    mass = new float[ncoord];

    // Read coordinates
    for (int k=0;k < nshoot;k++) {
      for (int i=0;i < ncoord;i++) {
	char atomstr[6], atomname[6], resname[8];
	int atomid, resid_tmp;
	double dummy1, dummy2;
	file >> atomstr >> atomid >> atomname >> resname >> resid_tmp
	     >> coord[k][i].x >> coord[k][i].y >> coord[k][i].z >> dummy1 >> dummy2;
	if (k == 0) {
	  // Set mass
	  if (atomname[0] == 'H') {
	    mass[i] = 1.008;
	  } else if (atomname[0] == 'N') {
	    mass[i] = 14.007;
	  } else if (atomname[0] == 'C') {
	    mass[i] = 12.011;
	  } else if (atomname[0] == 'O') {
	    mass[i] = 15.999;
	  } else {
	    std::cerr << "Unidentified atom type: i=" << i << " atomname="<< atomname << std::endl;
	    exit(1);
	  }
	  // Set residue ID
	  resid[i] = resid_tmp;
	}
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

