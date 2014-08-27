//
// Genetic Algorithm for collective variable (CV) finding
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include <iostream>
#include <algorithm>
#include <vector>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include "cv_util.h"
#include "LM.hpp"
#include "GA.hpp"

//
// Class constructor
//
GA::GA(Coord *coord, const int npair, const int *hA, const int *hB) {
  this->coord = coord;
  this->npair = npair;

  const int N = coord->get_nshoot();

  pair = new Pair(npair);
  pair_cv = new float[N*npair];

  // Seed random number generator
  rand_eng.seed(123456);

  alist = new int[N];
  blist = new int[N];
  AsBs(N, hA, hB, nalist, alist, nblist, blist);

  const int M = npair;
  zA = new double[nalist*M];
  zB = new double[nblist*M];
}

//
// Class destructor
//
GA::~GA() {
  delete pair;
  delete [] pair_cv;
  delete [] alist;
  delete [] blist;
  delete [] zA;
  delete [] zB;
}

//
// Initializes population
//
void GA::init_population() {

  /*
  std::uniform_int_distribution<> rand_dist(0,coord->get_nresidue()-1);

  for (int i=0;i < npair;i++) {
    // Pick randomly two residues
    int a = 0;
    int b = 0;
    //while (a == b) {
      a = rand_dist(rand_eng);
      b = rand_dist(rand_eng);
      //}

    if (a == b) {
      int start = coord->get_residue_start(a);
      int stop = coord->get_residue_stop(a);
      int n1 = (stop - start + 1)/2 - 1;
      pair->set(i, a, start, start+n1, a, start+n1+1, stop);
    } else {
      pair->set(i, a, coord->get_residue_start(a), coord->get_residue_stop(a),
		b, coord->get_residue_start(b), coord->get_residue_stop(b));
    }

  }
  */

  std::uniform_real_distribution<> pick_coord(0.0, 1.0);

  int ncoord = coord->get_ncoord();
  double coord_pick_rate = 5.0/(double)ncoord;

  // Glu 217 OE2 - Glu 212 HE2
  pair->add1(0, 27);
  pair->add2(0, 28);

  for (int i=1;i < npair;i++) {
    for (int j=0;j < ncoord;j++) {
      if (pick_coord(rand_eng) < coord_pick_rate) {
	pair->add1(i, j);
      }
      if (pick_coord(rand_eng) < coord_pick_rate) {
	pair->add2(i, j);
      }
    }
  }

  for (int i=0;i < 10;i++) {
    std::vector<int> *p1 = pair->get_atom1(i);
    std::vector<int> *p2 = pair->get_atom2(i);
    std::cout << "[ ";
    for (int i=0;i < p1->size();i++) {
      std::cout << (*p1)[i] << " ";
    }
    std::cout << "] [ ";
    for (int i=0;i < p2->size();i++) {
      std::cout << (*p2)[i] << " ";
    }
    std::cout << "]" << std::endl;
  }

}

//
// Evaluates Collective Variable (CV) values
//
void GA::eval_cv_values() {

  const int N = coord->get_nshoot();
  for (int ishoot=0;ishoot < N;ishoot++) {
    pair->eval(coord->get_coord(ishoot), coord->get_mass(), N, &pair_cv[ishoot]);
  }

  /*
  std::cout << "pair_cv[0]    =" << pair_cv[0] << std::endl;
  std::cout << "pair_cv[9999] =" << pair_cv[9999] << std::endl;
  std::cout << "pair_cv[10000]=" << pair_cv[10000] << std::endl;
  std::cout << "pair_cv[19999]=" << pair_cv[19999] << std::endl;
  */

}

/*
//
// Builds next generation
//
void GA::build_next_generation_old(const int ntop_pair, const int *top_pair_ind, Pair *top_pair) {

  std::uniform_int_distribution<> pick_pair(0,ntop_pair-1);
  std::uniform_int_distribution<> pick_action(0.0, 1.0);
  
  // Top contenders go through without changes
  for (int i=0;i < ntop_pair;i++) {
    pair->set(i, top_pair->get_resid1(i), top_pair->get_atom1(i),
	      top_pair->get_resid2(i), top_pair->get_atom2(i));
  }

  // Rest are created by mutating and mixing top contenders
  for (int i=ntop_pair;i < npair;i++) {
    if (pick_action(rand_eng) < 0.8) {
      //-----------------------------------------------------------------
      // Pick a random top contender and mutate it to produce offspring
      //-----------------------------------------------------------------
      int a = pick_pair(rand_eng);
      // Mutate by randomly flipping the cluster atoms on/off
      for (int d=0;d < 2;d++) {
	int resid = (d == 0) ? top_pair->get_resid1(a) : top_pair->get_resid2(a);
	int start = coord->get_residue_start(resid);
	int stop = coord->get_residue_stop(resid);
	for (int j=start;j <= stop;j++) {
	  if (pick_action(rand_eng) < 0.05) {
	    if (d == 0)
	      top_pair->flip_atom1(a, j);
	    else
	      top_pair->flip_atom2(a, j);
	  }
	}
      }
      pair->set(i, top_pair->get_resid1(a), top_pair->get_atom1(a),
		top_pair->get_resid2(a), top_pair->get_atom2(a));
    } else {
      //-------------------------------------------------------------------
      // Pick two random top contenders and pair them to produce offspring
      //-------------------------------------------------------------------
      int a = 0;
      int b = 0;
      while (a == b) {
	a = pick_pair(rand_eng);
	b = pick_pair(rand_eng);
      }
      // Mix a and b
      if (pick_action(rand_eng) < 0.5) {
	// Pick cluster 1 from a and cluster 2 from b
	pair->set(i, top_pair->get_resid1(a), top_pair->get_atom1(a),
		  top_pair->get_resid2(b), top_pair->get_atom2(b));
      } else {
	// Pick cluster 1 from b and cluster 2 from a
	pair->set(i, top_pair->get_resid1(b), top_pair->get_atom1(b),
		  top_pair->get_resid2(1), top_pair->get_atom2(1));
      }
    }
    
  }

}
*/

//
//
//
void GA::build_next_generation(std::vector<lnval_t> &lnval) {

  int ntop_pair = npair/10;
  int ncoord = coord->get_ncoord();

  std::uniform_int_distribution<> pick_pair(0,ntop_pair-1);
  std::uniform_real_distribution<> pick_action(0.0, 1.0);

  Pair new_pair(npair);

  std::cout << "-------------- New Top of the Crop -------------" << std::endl;
  std::vector<int> *p1 = pair->get_atom1(lnval[0].ind);
  std::vector<int> *p2 = pair->get_atom2(lnval[0].ind);
  std::cout << "[ ";
  for (int i=0;i < p1->size();i++) {
    std::cout << (*p1)[i] << " ";
  }
  std::cout << "] [ ";
  for (int i=0;i < p2->size();i++) {
    std::cout << (*p2)[i] << " ";
  }
  std::cout << "]" << std::endl;
  for (int i=0;i < 10;i++) {
    std::cout << lnval[i].ind << " " << lnval[i].data[2] << std::endl;
  }

  // Elitism, first 2 go through without changes
  for (int i=0;i < 2;i++) {
    int a = lnval[i].ind;
    new_pair.set(i, pair->get_atom1(a), pair->get_atom2(a));
  }

  for (int i=2;i < npair;i++) {
    if (pick_action(rand_eng) < 0.1) {
      // Mutate
      int a = lnval[pick_pair(rand_eng)].ind;

      new_pair.set(i, pair->get_atom1(a), pair->get_atom2(a));

      for (int j=0;j < ncoord;j++) {
	if (pick_action(rand_eng) < 0.01) {
	  new_pair.flip_atom1(i, j);
	}
	if (pick_action(rand_eng) < 0.01) {
	  new_pair.flip_atom2(i, j);
	}
      }
    } else {
      // Crossover
      int a = lnval[pick_pair(rand_eng)].ind;
      int b = lnval[pick_pair(rand_eng)].ind;
      for (int j=0;j < ncoord;j++) {
	if (pick_action(rand_eng) < 0.5) {
	  if (pair->contains1(a, j)) new_pair.add1(i, j);
	} else {
	  if (pair->contains1(b, j)) new_pair.add1(i, j);
	}
	if (pick_action(rand_eng) < 0.5) {
	  if (pair->contains2(a, j)) new_pair.add2(i, j);
	} else {
	  if (pair->contains2(b, j)) new_pair.add2(i, j);
	}
      }
    }
  }

  // Copy new_pair -> pair
  for (int i=0;i < npair;i++) {
    pair->set(i, new_pair.get_atom1(i), new_pair.get_atom2(i));
  }

}

bool compare_val(const lnval_t& a, const lnval_t& b) {
  return (a.data[2] > b.data[2]);
}

//
// Run genetic algorithm
//
void GA::run(const int niter) {

  //  int ntop_pair = npair/10;
  //int *top_pair_ind = new int[ntop_pair];
  //Pair top_pair(ntop_pair, coord->get_max_residue_size());

  for (int iter=0;iter < niter;iter++) {
    // Evaluate CV values
    eval_cv_values();

    // ----------------------------
    // Evaluate population fitness
    // ----------------------------

    // Prepare CVs
    const int M = npair;
    const int N = coord->get_nshoot();
    for(int i=0;i < M;i++){
      // Get maximum and minimum among all shooting points
      float qmin, qmax;
      getminmax(N, &pair_cv[i*N], qmin, qmax);
      // Reduce variables to [0, 1]
      for(int j=0;j < N;j++){
	pair_cv[i*N + j] = (pair_cv[i*N + j]-qmin)/(qmax-qmin);
      }
    }
    
    for(int i=0;i < M;i++){
      for(int j=0;j < nalist;j++){
	zA[j*M + i] = pair_cv[i*N + alist[j]];
      }
      for(int j=0;j < nblist;j++){
	zB[j*M + i] = pair_cv[i*N + blist[j]];
      }
    }

    double crit_move = 0.001;
    double crit_grad = 0.01;
    double crit_dlnL = 0.0001;
    double maxsize = 0.1;

    std::vector<lnval_t> lnval(npair);

#pragma omp parallel
    {
      LM<double> lm1(1);
      double alnLmax[1+2];
      // Loop through pairs
      int i;
#pragma omp for private(i) schedule(dynamic)
      for (i=0;i < npair;i++) {
	int cv[1];
	cv[0] = i;
	bool debug = false;
	if (pair->get_natom1(i) == 0 || pair->get_natom2(i) == 0) {
	  for (int jj=0;jj < 1+2;jj++) alnLmax[jj] = -1.0e10;
	} else {
	  lm1.calc_lm(debug, 1, cv, nalist, nblist, npair, zA, zB,
		      crit_move, crit_grad, crit_dlnL, maxsize, alnLmax);
	}
	lnval[i].ind = i;
	for (int jj=0;jj < 1+2;jj++) lnval[i].data[jj] = alnLmax[jj];
      }
    }

    // Sort values
    std::sort(lnval.begin(), lnval.end(), compare_val);
    
    /*
    // Select top 10% (with highest likelihood value)
    std::cout << "-------------- New Top of the Crop -------------" << std::endl;
    int nprint = (ntop_pair <= 10) ? ntop_pair : 10;
    for (int i=0;i < ntop_pair;i++) {
      if (i < nprint) {
	std::cout << lnval[i].ind << " " << lnval[i].data[2] << std::endl;
      }
      top_pair_ind[i] = lnval[i].ind;
      top_pair.set(i, pair->get_resid1(top_pair_ind[i]), pair->get_atom1(top_pair_ind[i]),
		   pair->get_resid2(top_pair_ind[i]), pair->get_atom2(top_pair_ind[i]));
    }

    // Build new generation
    build_next_generation_old(ntop_pair, top_pair_ind, &top_pair);
    */

    build_next_generation(lnval);

  }

  //delete [] top_pair_ind;

}
