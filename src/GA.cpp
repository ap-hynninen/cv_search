//
// Genetic Algorithm for collective variable (CV) finding
// (c) Antti-Pekka Hynninen, antti.pekka.hynninen@nrel.gov
//
#include <iostream>
#include <cassert>
#include <cmath>
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
GA::GA(Coord *coord, const int num_cv, const int ngenome, const double p_mutate, const int *hA, const int *hB) : 
  ngenome(ngenome), coord(coord), num_cv(num_cv), p_mutate(p_mutate) {
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

  int ncoord = coord->get_ncoord();

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

#ifdef USE_RANDOM
  std::uniform_int_distribution<> pick_coord(0, ncoord-1);
  std::uniform_real_distribution<> pick_type(0.0, 1.0);
#endif
  for (int i=istart;i < ngenome;i++) {
    for (int icv=0;icv < num_cv;icv++) {
      int a, b;
#ifdef USE_RANDOM
      a = pick_coord(rand_eng);
      b = pick_coord(rand_eng);
#else
      a = rand() % ncoord;
      b = rand() % ncoord;
#endif
      genome[i].get_gene(icv)->add(a);
      genome[i].get_gene(icv)->add(b);
      double r;
#ifdef USE_RANDOM
      r = pick_type(rand_eng);
#else
      r = (double)rand()/(double)RAND_MAX;
#endif
      if (r < 0.5) {
	// Build a triplet
	int c;
#ifdef USE_RANDOM
	c = pick_coord(rand_eng);
#else
	c = rand() % ncoord;
#endif
	genome[i].get_gene(icv)->add(c);
      }
    }	
  }

  /*
  int nprint = (ngenome >= 10) ? 10 : 1;
  for (int i=0;i < nprint;i++) {
    for (int icv=0;icv < num_cv;icv++) {
      genome[i].get_gene(icv)->print();
    }
  }
  */

}

//
// Advance GA
//
void GA::build_next_generation(std::vector<lnval_t> &lnval) {

  int ntop_genome = (ngenome >= 10) ? ngenome/10 : 1;
  int ncoord = coord->get_ncoord();

#ifdef USE_RANDOM
  std::uniform_int_distribution<> pick_genome(0,ntop_genome-1);
  std::uniform_int_distribution<> pick_coord(0,ncoord-1);
  std::uniform_real_distribution<> pick_action(0.0, 1.0);
#endif

  std::cout << "Best individual, fitness " << lnval[0].key << " :" << std::endl;
  for (int icv=0;icv < num_cv;icv++) {
    genome[lnval[0].ind].get_gene(icv)->print();
  }
  std::cout << "Rest of top 10:" << std::endl;
  int nprint = (ngenome >= 10) ? 10 : 1;
  for (int i=0;i < nprint;i++) {
    std::cout << lnval[i].key << std::endl;
  }

#define USE_RW
#ifdef USE_RW
  // Setup roulette wheel selection
  double* w = new double[ngenome];
  //for (int i=0;i < ngenome;i++) w[i] = 1.0/((-lnval[i].key)*(-lnval[i].key));
  // Set fac such that fac^ngenome = 0.01
  double fac = exp(log(0.01)/ngenome);
  w[0] = 1.0;
  for (int i=1;i < ngenome;i++) w[i] = fac*w[i-1];
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
      new_genome[i].get_gene(icv)->set(genome[a].get_gene(icv)->get_atom());
    }
  }

  for (int i=2;i < ngenome;i++) {
    double r;
#ifdef USE_RANDOM
    r = pick_action(rand_eng);
#else
    r = (double)rand()/(double)RAND_MAX;
#endif
    if (r < p_mutate) {
      // ---------------------------------------------
      // Mutate 
      // ---------------------------------------------
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
      rpair = rand() % ntop_genome;
#endif
#endif
      int a = lnval[rpair].ind;

      // Change one random atom
      int n = 0;
      for (int icv=0;icv < num_cv;icv++) {
	new_genome[i].get_gene(icv)->set(genome[a].get_gene(icv)->get_atom());
	n += genome[a].get_gene(icv)->get_natom();
      }
      int ipick;
#ifdef USE_RANDOM
      ipick = pick_coord(rand_eng) % n;
#else
      ipick = rand() % n;
#endif
      for (int icv=0;icv < num_cv;icv++) {
	if (ipick < genome[a].get_gene(icv)->get_natom()) {
	  int aa;
	  if (false) {
	    // Pick a neighbor of current atom
	    int ac = (*(new_genome[i].get_gene(icv)->get_atom()))[ipick];
	    int start = coord->get_nnlist_start(ac);
	    int end = coord->get_nnlist_end(ac);
	    int n = end - start + 1;
#ifdef USE_RANDOM
	    aa = pick_coord(rand_eng) % n;
#else
	    aa = rand() % n;
#endif
	    aa = coord->get_nnlist(start+aa);
	  } else {
	    // Pick a random atom
#ifdef USE_RANDOM
	    aa = pick_coord(rand_eng);
#else
	    aa = rand() % ncoord;
#endif
	  }
	  new_genome[i].get_gene(icv)->replace(ipick, aa);
	  break;
	}
	ipick -= genome[a].get_gene(icv)->get_natom();
      }

    } else {
      // ----------
      // Crossover
      // ----------
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
      rpair = rand() % ntop_genome;
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
      rpair = rand() % ntop_genome;
#endif
#endif
      int b = lnval[rpair].ind;

      for (int icv=0;icv < num_cv;icv++) {
#ifdef USE_RANDOM
	r = pick_action(rand_eng);
#else
	r = (double)rand()/(double)RAND_MAX;
#endif      
	if (r < 0.5) {
	  new_genome[i].get_gene(icv)->set(genome[a].get_gene(icv)->get_atom());
	} else {
	  new_genome[i].get_gene(icv)->set(genome[b].get_gene(icv)->get_atom());
	}
      }
    }
  }

  // Set genome = new_genome
  for (int i=0;i < ngenome;i++) {
    for (int icv=0;icv < num_cv;icv++) {
      genome[i].get_gene(icv)->set(new_genome[i].get_gene(icv)->get_atom());
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
	  genome[i].get_gene(icv)->eval(coord->get_coord(ishoot));
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
	bool skip_lm = false;
	for(int icv=0;icv < num_cv;icv++) {
	  skip_lm = (skip_lm || genome[i].get_gene(icv)->has_duplicate());
	}
	if (skip_lm) {
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

    std::cout << "----------------- Iteration " << iter << " -----------------" << std::endl;
    build_next_generation(lnval);

  }

}
