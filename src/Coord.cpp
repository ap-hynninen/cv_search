#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "../include/Coord.hpp"

//
// Class constructor
//
Coord::Coord(const char *coord_filename, const int ncoord, const int nshoot, const float rnn) {

  this->ncoord = ncoord;
  this->nshoot = nshoot;

  coord = NULL;
  resid = NULL;
  mass = NULL;
  residue_start = NULL;
  nnlist_pos = NULL;
  nnlist = NULL;

  load_coord(coord_filename);
  
  build_nnlist(rnn);

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
  if (nnlist_pos != NULL) delete [] nnlist_pos;
  if (nnlist != NULL) delete [] nnlist;
}

//
// Builds neighborlist using the first shooting point
//
void Coord::build_nnlist(const float rnn) {
  // Estimate number of neighbors each atom has
  float xmin = 1.0e10;
  float ymin = 1.0e10;
  float zmin = 1.0e10;
  float xmax = -1.0e10;
  float ymax = -1.0e10;
  float zmax = -1.0e10;
  for (int i=0;i < ncoord;i++) {
    xmin = (xmin < coord[0][i].x) ? xmin : coord[0][i].x;
    ymin = (ymin < coord[0][i].y) ? ymin : coord[0][i].y;
    zmin = (zmin < coord[0][i].z) ? zmin : coord[0][i].z;
    xmax = (xmax > coord[0][i].x) ? xmax : coord[0][i].x;
    ymax = (ymax > coord[0][i].y) ? ymax : coord[0][i].y;
    zmax = (zmax > coord[0][i].z) ? zmax : coord[0][i].z;
  }
  double vol = (xmax-xmin)*(ymax-ymin)*(zmax-zmin);
  nnlist_len = (int)((4.0/3.0*3.14159265358979323846*rnn*rnn*rnn)*((double)ncoord)*1.2/vol) + 1;
  nnlist_pos = new int[ncoord+1];
  nnlist = new int[nnlist_len];
  // Build neighborlist
  int pos = 0;
  float rnn2 = rnn*rnn;
  for (int i=0;i < ncoord;i++) {
    float xi = coord[0][i].x;
    float yi = coord[0][i].y;
    float zi = coord[0][i].z;
    nnlist_pos[i] = pos;
    for (int j=0;j < ncoord;j++) {
      if (i != j) {
	float xj = coord[0][j].x;
	float yj = coord[0][j].y;
	float zj = coord[0][j].z;
	float dx = xi-xj;
	float dy = yi-yj;
	float dz = zi-zj;
	float r2 = dx*dx + dy*dy + dz*dz;
	if (r2 < rnn2) {
	  // Reallocate if needed
	  if (pos == nnlist_len) {
	    int nnlist_len_new = (int)((double)nnlist_len*1.5);
	    int *nnlist_new = new int[nnlist_len_new];
	    for (int t=0;t < nnlist_len;t++) {
	      nnlist_new[t] = nnlist[t];
	    }
	    delete [] nnlist;
	    nnlist = nnlist_new;
	    nnlist_len = nnlist_len_new;
	  }
	  // Add to list
	  nnlist[pos++] = j;
	}
      }
    }
  }
  nnlist_pos[ncoord] = pos;
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

