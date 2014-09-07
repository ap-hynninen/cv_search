//
// Genetic Algorithm for collective variable (CV) finding
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include <iostream>
#include <cassert>
#include <algorithm>
#include <vector>
#ifdef USE_RANDOM
#include <random>
#else
#include <cstdlib>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include "../include/cv_util.h"
#include "../include/LM.hpp"
#include "../include/RW.hpp"
#include "../include/GA.hpp"

//
// Class constructor
//
GA::GA(Coord *coord, const int num_cv, const int ngenome, const int *hA, const int *hB) : 
  ngenome(ngenome), coord(coord), num_cv(num_cv) {
  assert(num_cv <= max_num_cv);

  const int nshoot = coord->get_nshoot();

  genome = new Genome[ngenome];
  genome_cv = new float[nshoot*ngenome*num_cv];
  for (int i=0;i < ngenome;i++) {
    genome[i].set_num_cv(num_cv);
  }

  // Seed random number generator
#ifdef USE_RANDOM
  rand_eng.seed(123456);
#else
  srand(123456);
#endif

  alist = new int[nshoot];
  blist = new int[nshoot];
  AsBs(nshoot, hA, hB, nalist, alist, nblist, blist);

  zA = new double[nalist*ngenome*num_cv];
  zB = new double[nblist*ngenome*num_cv];
}

//
// Class destructor
//
GA::~GA() {
  delete [] genome;
  delete [] genome_cv;
  delete [] alist;
  delete [] blist;
  delete [] zA;
  delete [] zB;
}

//
// Initializes population
//
void GA::init_population() {

#ifdef USE_RANDOM
  std::uniform_real_distribution<> pick_coord(0.0, 1.0);
#endif

  int ncoord = coord->get_ncoord();
  double coord_pick_rate = 5.0/(double)ncoord;

  int istart = 0;
  /*
  if (num_cv == 1) {
    // Glu 217 OE2 - Glu 217 HE2
    genome[0].get_pair(0)->add1(27);
    genome[0].get_pair(0)->add2(28);
    istart = 1;
  } else if (num_cv == 2) {
    // -1 anomeric carbon - Glu 212 OE2
    genome[0].get_pair(0)->add1(53-1);
    genome[0].get_pair(0)->add2(13-1);
    // -1 anomeric carbon to glycosidic oxygen
    genome[0].get_pair(1)->add1(53-1);
    genome[0].get_pair(1)->add2(55-1);
    istart = 1;
  } else if (num_cv == 3) {
    // -1 anomeric carbon - Glu 212 OE2
    genome[0].get_pair(0)->add1(53-1);
    genome[0].get_pair(0)->add2(13-1);
    // -1 anomeric carbon to glycosidic oxygen
    genome[0].get_pair(1)->add1(53-1);
    genome[0].get_pair(1)->add2(55-1);
    // Glu 217 OE2 - Glu 217 HE2
    genome[0].get_pair(2)->add1(28-1);
    genome[0].get_pair(2)->add2(29-1);
    istart = 1;
  }
  */

  for (int i=istart;i < ngenome;i++) {
    for (int icv=0;icv < num_cv;icv++) {
      for (int j=0;j < ncoord;j++) {
	double r;
#ifdef USE_RANDOM
	r = pick_coord(rand_eng);
#else
	r = (double)rand()/(double)RAND_MAX;
#endif
	if (r < coord_pick_rate) {
	  genome[i].get_pair(icv)->add1(j);
	}
#ifdef USE_RANDOM
	r = pick_coord(rand_eng);
#else
	r = (double)rand()/(double)RAND_MAX;
#endif
	if (r < coord_pick_rate) {
	  genome[i].get_pair(icv)->add2(j);
	}
      }
    }

    /*
    if (genome[i].get_pair(0)->get_natom1() <= 5 &&
	genome[i].get_pair(0)->get_natom2() <= 5) {
      i++;
    }
    */

  }

  int nprint = (ngenome >= 10) ? 10 : 1;
  for (int i=0;i < nprint;i++) {
    for (int icv=0;icv < num_cv;icv++) {
      genome[i].get_pair(icv)->print();
    }
  }

}

//
// Advance GA
//
void GA::build_next_generation(std::vector<lnval_t> &lnval) {

  int ntop_genome = (ngenome >= 10) ? ngenome/10 : 1;
  int ncoord = coord->get_ncoord();

#ifdef USE_RANDOM
  std::uniform_int_distribution<> pick_genome(0,ntop_genome-1);
  std::uniform_real_distribution<> pick_action(0.0, 1.0);
#endif

  std::cout << "-------------- New Top of the Crop -------------" << std::endl;
  for (int icv=0;icv < num_cv;icv++) {
    genome[lnval[0].ind].get_pair(icv)->print();
  }
  int nprint = (ngenome >= 10) ? 10 : 1;
  for (int i=0;i < nprint;i++) {
    std::cout << lnval[i].ind << " " << lnval[i].key << std::endl;
  }

  //#define USE_RW
#ifdef USE_RW
  // Setup roulette wheel selection
  double* w = new double[ngenome];
  for (int i=0;i < ngenome;i++) w[i] = 1.0/((-lnval[i].key)*(-lnval[i].key));
  RW rw(ngenome, w);
  delete [] w;
#endif

  Genome new_genome[ngenome];
  for (int i=0;i < ngenome;i++) {
    new_genome[i].set_num_cv(num_cv);
  }

  // Elitism, first 2 go through without changes
  for (int i=0;i < 2;i++) {
    int a = lnval[i].ind;
    for (int icv=0;icv < num_cv;icv++) {
      new_genome[i].get_pair(icv)->set(genome[a].get_pair(icv)->get_atom1(),
				       genome[a].get_pair(icv)->get_atom2());
    }
  }

  double mutate_rate = 0.4/(double)ncoord;

  int i=2;
  //  for (int i=2;i < ngenome;i++) {
  while (i < ngenome) {
    double r;
#ifdef USE_RANDOM
    r = pick_action(rand_eng);
#else
    r = (double)rand()/(double)RAND_MAX;
#endif
    //if (r < 0.2 || i == ngenome-1) {
    if (r < 0.2) {
      // Mutate
      int rpair;
#ifdef USE_RW
#ifdef USE_RANDOM
      rpair = rw.pick(rand_eng);
#else
      rpair = rw.pick();
#endif
#else
#ifdef USE_RANDOM
      rpair = pick_genome(rand_eng);
#else
      rpair = rand() % ntop_pair;
#endif
#endif
      int a = lnval[rpair].ind;

      for (int icv=0;icv < num_cv;icv++) {
	new_genome[i].get_pair(icv)->set(genome[a].get_pair(icv)->get_atom1(),
					 genome[a].get_pair(icv)->get_atom2());

	for (int j=0;j < ncoord;j++) {
#ifdef USE_RANDOM
	  r = pick_action(rand_eng);
#else
	  r = (double)rand()/(double)RAND_MAX;
#endif      
	  if (r < mutate_rate) {
	    new_genome[i].get_pair(icv)->flip_atom1(j);
	  }
#ifdef USE_RANDOM
	  r = pick_action(rand_eng);
#else
	  r = (double)rand()/(double)RAND_MAX;
#endif      
	  if (r < mutate_rate) {
	    new_genome[i].get_pair(icv)->flip_atom2(j);
	  }
	}
      }
      i++;
    } else {
      // Crossover
      int rpair;
#ifdef USE_RW
#ifdef USE_RANDOM
      rpair = rw.pick(rand_eng);
#else
      rpair = rw.pick();
#endif
#else
#ifdef USE_RANDOM
      rpair = pick_genome(rand_eng);
#else
      rpair = rand() % ntop_pair;
#endif
#endif
      int a = lnval[rpair].ind;

#ifdef USE_RW
#ifdef USE_RANDOM
      rpair = rw.pick(rand_eng);
#else
      rpair = rw.pick();
#endif
#else
#ifdef USE_RANDOM
      rpair = pick_genome(rand_eng);
#else
      rpair = rand() % ntop_pair;
#endif
#endif
      int b = lnval[rpair].ind;

      for (int icv=0;icv < num_cv;icv++) {
	/*
#ifdef USE_RANDOM
	r = pick_action(rand_eng);
#else
	r = (double)rand()/(double)RAND_MAX;
#endif      
	if (r < 0.25) {
	  new_genome[i].get_pair(icv)->set(genome[a].get_pair(icv)->get_atom1(),
					   genome[a].get_pair(icv)->get_atom2());
	  new_genome[i+1].get_pair(icv)->set(genome[b].get_pair(icv)->get_atom1(),
					     genome[b].get_pair(icv)->get_atom2());
	} else if (r < 0.5) {
	  new_genome[i].get_pair(icv)->set(genome[a].get_pair(icv)->get_atom1(),
					   genome[b].get_pair(icv)->get_atom2());
	  new_genome[i+1].get_pair(icv)->set(genome[b].get_pair(icv)->get_atom1(),
					     genome[a].get_pair(icv)->get_atom2());
	} else if (r < 0.75) {
	  new_genome[i].get_pair(icv)->set(genome[b].get_pair(icv)->get_atom1(),
					   genome[a].get_pair(icv)->get_atom2());
	  new_genome[i+1].get_pair(icv)->set(genome[a].get_pair(icv)->get_atom1(),
					     genome[b].get_pair(icv)->get_atom2());
	} else {
	  new_genome[i].get_pair(icv)->set(genome[b].get_pair(icv)->get_atom1(),
					   genome[b].get_pair(icv)->get_atom2());
	  new_genome[i+1].get_pair(icv)->set(genome[a].get_pair(icv)->get_atom1(),
					     genome[a].get_pair(icv)->get_atom2());
	}
	*/


	for (int j=0;j < ncoord;j++) {
#ifdef USE_RANDOM
	  r = pick_action(rand_eng);
#else
	  r = (double)rand()/(double)RAND_MAX;
#endif      
	  if (r < 0.5) {
	    if (genome[a].get_pair(icv)->contains1(j)) new_genome[i].get_pair(icv)->add1(j);
	  } else {
	    if (genome[b].get_pair(icv)->contains1(j)) new_genome[i].get_pair(icv)->add1(j);
	  }
#ifdef USE_RANDOM
	  r = pick_action(rand_eng);
#else
	  r = (double)rand()/(double)RAND_MAX;
#endif      
	  if (r < 0.5) {
	    if (genome[a].get_pair(icv)->contains2(j)) new_genome[i].get_pair(icv)->add2(j);
	  } else {
	    if (genome[b].get_pair(icv)->contains2(j)) new_genome[i].get_pair(icv)->add2(j);
	  }
	}

      }
      //i += 2;
      i++;
    }

    /*
    if (new_genome[i].get_pair(0)->get_natom1() <= 5 &&
	new_genome[i].get_pair(0)->get_natom2() <= 5) {
      i++;
    }
    */

  }

  for (int i=0;i < ngenome;i++) {
    for (int icv=0;icv < num_cv;icv++) {
      genome[i].get_pair(icv)->set(new_genome[i].get_pair(icv)->get_atom1(),
				   new_genome[i].get_pair(icv)->get_atom2());
    }
  }

}

//
// Evaluates Collective Variable (CV) values
//
void GA::eval_cv_values() {

  const int nshoot = coord->get_nshoot();
  for (int ishoot=0;ishoot < nshoot;ishoot++) {
    for (int i=0;i < ngenome;i++) {
      for (int icv=0;icv < num_cv;icv++) {
	genome_cv[ishoot + (i*num_cv + icv)*nshoot] = 
	  genome[i].get_pair(icv)->eval(coord->get_coord(ishoot), coord->get_mass());
      }
    }
  }

}


bool compare_val(const lnval_t& a, const lnval_t& b) {
  //return (a.data[2] > b.data[2]);
  return (a.key > b.key);
}

//
// Run genetic algorithm
//
void GA::run(const int niter) {

  for (int iter=0;iter < niter;iter++) {
    // Evaluate CV values
    eval_cv_values();

    // ----------------------------
    // Evaluate population fitness
    // ----------------------------

    // Prepare CVs
    const int nshoot = coord->get_nshoot();
    for(int i=0;i < ngenome*num_cv;i++) {
      // Get maximum and minimum among all shooting points
      float qmin, qmax;
      getminmax(nshoot, &genome_cv[i*nshoot], qmin, qmax);
      // Reduce variables to [0, 1]
      for(int j=0;j < nshoot;j++) {
	genome_cv[i*nshoot + j] =  (genome_cv[i*nshoot + j]-qmin)/(qmax-qmin);
      }
    }

    // alist[0...nalist-1] = list of shooting points for A
    // blist[0...nblist-1] = list of shooting points for B
    for(int i=0;i < ngenome*num_cv;i++) {
      for(int j=0;j < nalist;j++) {
	zA[j*ngenome*num_cv + i] = genome_cv[i*nshoot + alist[j]];
      }
      for(int j=0;j < nblist;j++){
	zB[j*ngenome*num_cv + i] = genome_cv[i*nshoot + blist[j]];
      }
    }

    double crit_move = 0.001;
    double crit_grad = 0.01;
    double crit_dlnL = 0.0001;
    double maxsize = 0.1;

    std::vector<lnval_t> lnval(ngenome);

#pragma omp parallel
    {
      LM<double> lm(num_cv);
      double alnLmax[num_cv+2];
      // Loop through genomes
      int i;
#pragma omp for private(i) schedule(dynamic)
      for (i=0;i < ngenome;i++) {
	int cv[num_cv];
	for(int icv=0;icv < num_cv;icv++) cv[icv] = i*num_cv + icv;
	bool zero_natom = false;
	for(int icv=0;icv < num_cv;icv++) {
	  if (genome[i].get_pair(icv)->get_natom1() == 0 ||
	      genome[i].get_pair(icv)->get_natom2() == 0) zero_natom = true;
	  if (genome[i].get_pair(icv)->get_natom1() >= 5 ||
	      genome[i].get_pair(icv)->get_natom2() >= 5) zero_natom = true;
	}
	if (zero_natom) {
	  for (int jj=0;jj < num_cv+2;jj++) alnLmax[jj] = -1.0e10;
	} else {
	  lm.calc_lm(false, num_cv, cv, nalist, nblist, ngenome*num_cv, zA, zB,
		     crit_move, crit_grad, crit_dlnL, maxsize, alnLmax);
	}
	lnval[i].ind = i;
	for (int jj=0;jj < num_cv+2;jj++) lnval[i].data[jj] = alnLmax[jj];
	lnval[i].key = alnLmax[num_cv+1];
      }
    }

    // Sort values
    std::sort(lnval.begin(), lnval.end(), compare_val);

    build_next_generation(lnval);

  }

}
