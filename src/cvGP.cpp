#include <iostream>
#include <fstream>
#include <strstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include "../include/struct.h"
#include "../gpc/include/gp.h"
#include "../gpc/include/gpconfig.h"
#include "../include/Coord.hpp"
#include "../include/cv_util.h"
#include "../include/LM.hpp"

// Functions
const int ADD  = -1;
const int SUB  = -2;
const int MUL  = -3;
const int SQRT = -4;
const int DIV  = -5;
const int SQR  = -6;
const int ADF0 = -7;
const int X1   = -8;
const int X2   = -9;
const int SQRSUB = -10;
const int ABS   = -11;
const int ANGLE = -12;
const int DIHE  = -13;

// Define class identifiers
const int cvGeneID=GPUserID;
const int cvGPID=GPUserID+1;
const int cvPopulationID=GPUserID+2;

Coord *cvCoord = NULL;
int ngenome;
int num_cv;
LM<double> **lm_th = NULL;
double **alnLmax_th = NULL;
int *cv = NULL;
float **genome_cv_th = NULL;
int nalist;
int *alist = NULL;
int nblist;
int *blist = NULL;
double **zA_th = NULL;
double **zB_th = NULL;
float *correct_cv = NULL;
int *num_illegal = NULL;

// Define configuration parameters and the neccessary array to
// read/write the configuration to a file. If you need more variables,
// just add them below and insert an entry in the configArray.
GPVariables cfg;
char *InfoFileName="data";
struct GPConfigVarInformation configArray[]=
{
  {"PopulationSize", DATAINT, &cfg.PopulationSize},
  {"NumberOfGenerations", DATAINT, &cfg.NumberOfGenerations},
  {"CreationType", DATAINT, &cfg.CreationType},
  {"CrossoverProbability", DATADOUBLE, &cfg.CrossoverProbability},
  {"CreationProbability", DATADOUBLE, &cfg.CreationProbability},
  {"MaximumDepthForCreation", DATAINT, &cfg.MaximumDepthForCreation},
  {"MaximumDepthForCrossover", DATAINT, &cfg.MaximumDepthForCrossover},
  {"SelectionType", DATAINT, &cfg.SelectionType},
  {"TournamentSize", DATAINT, &cfg.TournamentSize},
  {"DemeticGrouping", DATAINT, &cfg.DemeticGrouping},
  {"DemeSize", DATAINT, &cfg.DemeSize},
  {"DemeticMigProbability", DATADOUBLE, &cfg.DemeticMigProbability},
  {"SwapMutationProbability", DATADOUBLE, &cfg.SwapMutationProbability},
  {"ShrinkMutationProbability", DATADOUBLE, &cfg.ShrinkMutationProbability},
  {"SteadyState", DATAINT, &cfg.SteadyState},
  {"AddBestToNewPopulation", DATAINT, &cfg.AddBestToNewPopulation},
  {"InfoFileName", DATASTRING, &InfoFileName},
  {"", DATAINT, NULL}
};

class cvGP;

class cvGene : public GPGene
{
public:
  // Duplication (mandatory)
  cvGene (const cvGene& gpo) : GPGene (gpo) { }
  virtual GPObject& duplicate () { return *(new cvGene(*this)); }

  // Creation of own class objects (mandatory)
  cvGene (GPNode& gpo) : GPGene (gpo) {}
  virtual GPGene* createChild (GPNode& gpo) {
    return new cvGene (gpo); }

  // Tree evaluation (not mandatory, but somehow the trees must be
  // parsed to evaluate the fitness)
  float evaluate(float *coord, cvGP& gp, float arg0, float arg1);
  float3 evaluate_vec(float3 *coord);

  // Load and save (not mandatory)
  cvGene () {}
  virtual int isA () { return cvGeneID; }
  virtual GPObject* createObject() { return new cvGene; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Access children (not mandatory)
  cvGene* NthMyChild (int n) {
    return (cvGene*) GPContainer::Nth (n); }
};

class cvGP : public GP 
{
public:
  // Duplication (mandatory)
  cvGP (cvGP& gpo) : GP (gpo) { }
  virtual GPObject& duplicate () { return *(new cvGP(*this)); }

  // Creation of own class objects (mandatory)
  cvGP (int genes) : GP (genes) {}
  virtual GPGene* createGene (GPNode& gpo) {
    return new cvGene (gpo); }

  // Tree evaluation (mandatory)
  virtual void evaluate ();

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Load and save (not mandatory)
  cvGP () {}
  virtual int isA () { return cvGPID; }
  virtual GPObject* createObject() { return new cvGP; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);

  // Access trees (not mandatory)
  cvGene* NthMyGene (int n) {
    return (cvGene*) GPContainer::Nth (n); }
};


class cvPopulation : public GPPopulation
{
public:
  // Constructor (mandatory)
  cvPopulation (GPVariables& GPVar_, GPAdfNodeSet& adfNs_) : 
    GPPopulation (GPVar_, adfNs_) {}

  // Duplication (mandatory)
  cvPopulation (cvPopulation& gpo) : GPPopulation(gpo) {}
  virtual GPObject& duplicate () { return *(new cvPopulation(*this)); }

  // Creation of own class objects (mandatory)
  virtual GP* createGP (int numOfTrees) { return new cvGP (numOfTrees); }

  // Load and save (not mandatory)
  cvPopulation () {}
  virtual int isA () { return cvPopulationID; }
  virtual GPObject* createObject() { return new cvPopulation; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Access genetic programs (not mandatory)
  cvGP* NthMyGP (int n) {
    return (cvGP*) GPContainer::Nth (n); }
};

// This function is the divide with closure. Basically if you divide
// anything by zero you get an error so we have to stop this
// process. We check for a very small denominator and return a very
// high value.
inline double div (double x, double y)
{
  if (fabs (y)<1e-6)
    {
      if (x*y<0.0) 
	return -1e6;
      else
	return 1e6;
    }
  else
    return x/y;
}

inline float divf (float x, float y)
{
  if (fabsf (y)<1e-6)
    {
      if (x*y<0.0f) 
	return -1e6;
      else
	return 1e6;
    }
  else
    return x/y;
}

//
// Protected sqrt
//
inline float psqrtf(float x) {
  return sqrtf(fabsf(x));
}

//
// Evaluate
//
float cvGene::evaluate(float *coord, cvGP& gp, float arg0, float arg1) {
  float c, t, a0, a1;
  float3 v1, v2, v3;
  float v1abs, v3abs;
  float dx, dy, dz;

  if (isFunction ()) {
    switch (node->value ()) {
    case ABS:
      v1 = NthMyChild(0)->evaluate_vec((float3 *)coord);
      v2 = NthMyChild(1)->evaluate_vec((float3 *)coord);
      dx = v1.x - v2.x;
      dy = v1.y - v2.y;
      dz = v1.z - v2.z;
      c = sqrtf(dx*dx + dy*dy + dz*dz);
      break;
    case ANGLE:
      v1 = NthMyChild(0)->evaluate_vec((float3 *)coord);
      v2 = NthMyChild(1)->evaluate_vec((float3 *)coord);
      v3 = NthMyChild(2)->evaluate_vec((float3 *)coord);
      v1.x -= v2.x;
      v1.y -= v2.y;
      v1.z -= v2.z;
      v3.x -= v2.x;
      v3.y -= v2.y;
      v3.z -= v2.z;
      v1abs = sqrtf(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);
      v3abs = sqrtf(v3.x*v3.x + v3.y*v3.y + v3.z*v3.z);
      if (v1abs < 1.0e-14 || v3abs < 1.0e-14) {
	c = 0.0f;
      } else {
	c = acosf((v1.x*v3.x + v1.y*v3.y + v1.z*v3.z)/(v1abs*v3abs));
      }
      break;
    case ADD:
      c = NthMyChild(0)->evaluate(coord, gp, arg0, arg1) + NthMyChild(1)->evaluate(coord, gp, arg0, arg1);
      break;	
    case SUB:
      c = NthMyChild(0)->evaluate(coord, gp, arg0, arg1) - NthMyChild(1)->evaluate(coord, gp, arg0, arg1);
      break;	
    case MUL:
      c = NthMyChild(0)->evaluate(coord, gp, arg0, arg1) * NthMyChild(1)->evaluate(coord, gp, arg0, arg1);
      break;
    case SQR:
      t = NthMyChild(0)->evaluate(coord, gp, arg0, arg1);
      c = t*t;
      break;
    case ADF0:
      a0 = NthMyChild(0)->evaluate(coord, gp, arg0, arg1);
      a1 = NthMyChild(1)->evaluate(coord, gp, arg0, arg1);
      c = gp.NthMyGene(1)->evaluate(coord, gp, a0, a1);
      break;
    case SQRT:
      c = psqrtf(NthMyChild(0)->evaluate(coord, gp, arg0, arg1));
      break;
    case DIV:
      c = divf(NthMyChild(0)->evaluate(coord, gp, arg0, arg1),
	       NthMyChild(1)->evaluate(coord, gp, arg0, arg1));
      break;
    case SQRSUB:
      c = NthMyChild(0)->evaluate(coord, gp, arg0, arg1) - NthMyChild(1)->evaluate(coord, gp, arg0, arg1);
      c *= c;
      break;
    default: 
      GPExitSystem ("cvGene::evaluate", "Invalid node function");
    }
  }

  if (isTerminal ()) {
    switch (node->value()) {
    case X1:
      c = arg0;
      break;
    case X2:
      c = arg1;
      break;
    default:
      c = 0.0f; //coord[node->value()];
      break;
    }
  }

  return c;
}

float3 cvGene::evaluate_vec(float3 *coord) {
  if (isFunction()) {
    GPExitSystem("cvGene::evaluate_vec","No functions allowed!");
  }
  if (isTerminal()) {
    return coord[node->value()];
  }

  GPExitSystem("cvGene::evaluate_vec","Should not end up here");
  exit(1);
}

// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void cvGP::evaluate ()
{

  const int nshoot = cvCoord->get_nshoot();

  int tid = 0;
#ifdef USE_OPENMP
  // NOTE: omp_get_thread_num() returns 0 outside parallel regions
  tid = omp_get_thread_num();
#endif

  // Get pointers to memory for this thread
  float *genome_cv = genome_cv_th[tid];
  double *zA = zA_th[tid];
  double *zB = zB_th[tid];
  double *alnLmax = alnLmax_th[tid];
  LM<double> *lm = lm_th[tid];

  // Evaluate main tree for each shooting point
  for (int ishoot=0;ishoot < nshoot;ishoot++) {
    for (int icv=0;icv < num_cv;icv++) {
      genome_cv[ishoot + icv*nshoot] = 
	NthMyGene(icv)->evaluate((float *)cvCoord->get_coord(ishoot), *this, 0.0f, 0.0f);
      if (genome_cv[ishoot + icv*nshoot] == 0.0f) {
	stdFitness = 1.0e6;
	num_illegal[tid]++;
	return;
      }
    }
  }
  
  /*
  // ---------------------------
  // Try to reproduce |x12-x53|
  // ---------------------------
  stdFitness = 0.0;
  for (int ishoot=0;ishoot < nshoot;ishoot++) {
    stdFitness += fabs(genome_cv[ishoot] - correct_cv[ishoot]);
  }
  */

  // --------------------------------
  // Compute Likelihood maximization
  // --------------------------------

  // Prepare CVs
  // Get maximum and minimum among all shooting points
  for (int icv=0;icv < num_cv;icv++) {
    float qmin, qmax;
    getminmax(nshoot, &genome_cv[icv*nshoot], qmin, qmax);
    if (qmin == 0.0f && qmax == 0.0f) {
      stdFitness = 1.0e6;
      return;
    }
    // Reduce variables to [0, 1]
    for(int ishoot=0;ishoot < nshoot;ishoot++) {
      genome_cv[ishoot + icv*nshoot] =  (genome_cv[ishoot + icv*nshoot]-qmin)/(qmax-qmin);
    }
  }

  // alist[0...nalist-1] = list of shooting points for A
  // blist[0...nblist-1] = list of shooting points for B
  for (int icv=0;icv < num_cv;icv++) {
    for(int j=0;j < nalist;j++) {
      zA[j*num_cv + icv] = genome_cv[alist[j] + icv*nshoot];
    }
    for(int j=0;j < nblist;j++){
      zB[j*num_cv + icv] = genome_cv[blist[j] + icv*nshoot];
    }
  }
  
  double crit_move = 0.001;
  double crit_grad = 0.01;
  double crit_dlnL = 0.0001;
  double maxsize = 0.1;
  
  lm->calc_lm(false, num_cv, cv, nalist, nblist, num_cv, zA, zB,
	      crit_move, crit_grad, crit_dlnL, maxsize, alnLmax);
  
  //lnval[i].ind = i;
  //for (int jj=0;jj < num_cv+2;jj++) lnval[i].data[jj] = alnLmax[jj];
  //lnval[i].key = alnLmax[num_cv+1];
 
  stdFitness = 0.0 - alnLmax[num_cv+1];
}

// Create function and terminal set
void createNodeSet (GPAdfNodeSet& adfNs, const int ncoord)
{

  if (true) {
    // Reserve space for the node sets
    adfNs.reserveSpace (num_cv);

    // Value branches
    for (int icv=0;icv < num_cv;icv++) {
      GPNodeSet& ns0 = *new GPNodeSet (2 + ncoord);
      adfNs.put (icv, ns0);

      ns0.putNode (*new GPNode (ABS,  "ABS",  2));
      ns0.putNode (*new GPNode (ANGLE,"ANGLE",3));
      for (int i=0;i < ncoord;i++) {
	char istr[16];
	sprintf(istr,"%d",i);
	ns0.putNode (*new GPNode (i, istr));
      }
    }

  }

  if (false) {
    // Reserve space for the node sets
    adfNs.reserveSpace (1);
    
    // Now define the function and terminal set for each ADF and place
    // function/terminal sets into overall ADF container
    GPNodeSet& ns0=*new GPNodeSet (2 + 6);
    adfNs.put (0, ns0);

    ns0.putNode (*new GPNode (ADD,  "ADD",  2));
    /*
    ns0.putNode (*new GPNode (SUB,  "SUB",  2));
    ns0.putNode (*new GPNode (MUL,  "MUL",  2));
    ns0.putNode (*new GPNode (DIV,  "DIV",  2));
    ns0.putNode (*new GPNode (SQR,  "SQR",  1));
    ns0.putNode (*new GPNode (SQRT, "SQRT", 1));
    */
    ns0.putNode (*new GPNode (SQRSUB,  "SQRSUB",  2));
    // 12:
    ns0.putNode (*new GPNode (36, "12.x"));
    ns0.putNode (*new GPNode (37, "12.y"));
    ns0.putNode (*new GPNode (38, "12.z"));
    // 53:
    ns0.putNode (*new GPNode (159, "53.x"));
    ns0.putNode (*new GPNode (160, "53.y"));
    ns0.putNode (*new GPNode (161, "53.z"));

  }
  if (false) {
    // Reserve space for the node sets
    adfNs.reserveSpace (2);
    
    // Now define the function and terminal set for each ADF and place
    // function/terminal sets into overall ADF container
    GPNodeSet& ns0=*new GPNodeSet (7 + 4);
    GPNodeSet& ns1=*new GPNodeSet (6 + 2);
    adfNs.put (0, ns0);
    adfNs.put (1, ns1);
    
    // Value branch
    ns0.putNode (*new GPNode (ADD,  "ADD",  2));
    ns0.putNode (*new GPNode (SUB,  "SUB",  2));
    ns0.putNode (*new GPNode (MUL,  "MUL",  2));
    ns0.putNode (*new GPNode (SQR,  "SQR",  1));
    ns0.putNode (*new GPNode (ADF0, "ADF0", 2));
    ns0.putNode (*new GPNode (DIV,  "DIV",  2));
    ns0.putNode (*new GPNode (SQRT, "SQRT", 1));
    // 12:
    ns0.putNode (*new GPNode (36, "12.x"));
    ns0.putNode (*new GPNode (37, "12.y"));
    //ns0.putNode (*new GPNode (38, "12.z"));
    // 53:
    ns0.putNode (*new GPNode (159, "53.x"));
    ns0.putNode (*new GPNode (160, "53.y"));
    //ns0.putNode (*new GPNode (161, "53.z"));
    
    // ADF0
    ns1.putNode (*new GPNode (ADD,  "ADD",  2));
    ns1.putNode (*new GPNode (SUB,  "SUB",  2));
    ns1.putNode (*new GPNode (MUL,  "MUL",  2));
    ns1.putNode (*new GPNode (SQR,  "SQR",  1));
    ns1.putNode (*new GPNode (DIV,  "DIV",  2));
    ns1.putNode (*new GPNode (SQRT, "SQRT", 1));
    ns1.putNode (*new GPNode (X1,  "X1"));
    ns1.putNode (*new GPNode (X2,  "X2"));
  }

  /*
  for (int i=0;i < 3*ncoord;i++) {
    char istr[16];
    sprintf(istr,"%d",i);
    ns0.putNode (*new GPNode (i, istr));
  }
  */

}

// Randomly rotate and translate coordinates
void translate_rotate_coord(Coord *coord) {
  const int nshoot = coord->get_nshoot();
  const int ncoord = coord->get_ncoord();

  for (int ishoot=0;ishoot < nshoot;ishoot++) {
    float tx = ((double)rand()/(double)RAND_MAX - 0.5)*10.0;
    float ty = ((double)rand()/(double)RAND_MAX - 0.5)*10.0;
    float tz = ((double)rand()/(double)RAND_MAX - 0.5)*10.0;
    float3 *v = coord->get_coord(ishoot);
    for (int i=0;i < ncoord;i++) {
      v[i].x += tx;
      v[i].y += ty;
      v[i].z += tz;
    }
  }

}

void calc_correct_cv(Coord *coord) {
  const int nshoot = coord->get_nshoot();
  for (int ishoot=0;ishoot < nshoot;ishoot++) {
    float3 *v = coord->get_coord(ishoot);
    float x = v[12].x - v[53].x;
    float y = v[12].y - v[53].y;
    float z = v[12].z - v[53].z;
    correct_cv[ishoot] = x*x + y*y + z*z;
  }
}

void doGP(Coord *coord, const int num_cv_in, const int *hA, const int *hB) {

  //translate_rotate_coord(coord);
  correct_cv = new float[coord->get_nshoot()];
  calc_correct_cv(coord);

  for (int i=0;i < 10;i++) printf("%f\n",correct_cv[i]);

  num_cv = num_cv_in;
  cvCoord = coord;

  // Init GP system.
  GPInit (1, -1);
  
  // Declare the GP Variables, set defaults and read configuration
  // file.  The defaults will be overwritten by the configuration file
  // when read.  If it doesn't exist, the defaults will be written to
  // the file.
  GPConfiguration config (cout, "cv.ini", configArray);

  // Open the main output file for data and statistics file. First set
  // up names for data file.
  // Remember we should delete the string from the stream, well just a
  // few bytes
  ostrstream strOutFile, strStatFile, strTeXFile;
  strOutFile  << InfoFileName << ".dat" << ends;
  strStatFile << InfoFileName << ".stc" << ends;
  ofstream fout (strOutFile.str());
  ofstream bout (strStatFile.str());
  
  // Print the configuration to the files just opened
  fout << cfg << endl;
  cout << cfg << endl;

  // Create the adf function/terminal set and print it out.
  GPAdfNodeSet adfNs;
  createNodeSet (adfNs, coord->get_ncoord());
  cout << adfNs << endl;

  ngenome = cfg.PopulationSize;

  alist = new int[coord->get_nshoot()];
  blist = new int[coord->get_nshoot()];
  AsBs(coord->get_nshoot(), hA, hB, nalist, alist, nblist, blist);

  int nthread = 1;
#ifdef USE_OPENMP
#pragma omp parallel
  {
    if (omp_get_thread_num() == 0) nthread = omp_get_num_threads();
  }
#endif

  lm_th        = new LM<double>*[nthread];
  alnLmax_th   = new double*[nthread];
  genome_cv_th = new float*[nthread];
  zA_th        = new double*[nthread];
  zB_th        = new double*[nthread];

  for (int tid=0;tid < nthread;tid++) {
    lm_th[tid]        = new LM<double>(num_cv);
    alnLmax_th[tid]   = new double[num_cv+2];
    genome_cv_th[tid] = new float[coord->get_nshoot()*num_cv];
    zA_th[tid]        = new double[nalist*num_cv];
    zB_th[tid]        = new double[nblist*num_cv];
  }

  num_illegal = new int[nthread];
  for (int tid=0;tid < nthread;tid++) num_illegal[tid] = 0;
  cv = new int[num_cv];
  for(int icv=0;icv < num_cv;icv++) cv[icv] = icv;

  // Create a population with this configuration
  cout << "Creating initial population ..." << endl;
  cvPopulation* pop=new cvPopulation (cfg, adfNs);
  pop->create ();
  cout << "Ok." << endl;
  pop->createGenerationReport (1, 0, fout, bout);

  //pop->printOn(cout);

  // This next for statement is the actual genetic programming system
  // which is in essence just repeated reproduction and crossover loop
  // through all the generations .....
  cvPopulation* newPop=NULL;
  for (int gen=1; gen<=cfg.NumberOfGenerations; gen++) {
    for (int tid=0;tid < nthread;tid++) num_illegal[tid] = 0;
    // Create a new generation from the old one by applying the
    // genetic operators
    if (!cfg.SteadyState)
      newPop=new cvPopulation (cfg, adfNs);
    pop->generate (*newPop);
    
    // Delete the old generation and make the new the old one
    if (!cfg.SteadyState) {
      delete pop;
      pop=newPop;
    }
    
    // Create a report of this generation and how well it is doing
    pop->createGenerationReport (0, gen, fout, bout);

    int num_illegal_tot = 0;
    for (int tid=0;tid < nthread;tid++) num_illegal_tot += num_illegal[tid];
    std::cout << "Number of illegal trees = " << num_illegal_tot << std::endl;
  }
  
  cout << "\nResults are in " 
       << InfoFileName << ".dat," 
       << InfoFileName << ".stc." << endl;

  for (int tid=0;tid < nthread;tid++) {
    delete lm_th[tid];
    delete [] alnLmax_th[tid];
    delete [] genome_cv_th[tid];
    delete [] zA_th[tid];
    delete [] zB_th[tid];
  }

  delete [] lm_th;
  delete [] alnLmax_th;
  delete [] genome_cv_th;
  delete [] zA_th;
  delete [] zB_th;
  delete [] correct_cv;
  delete [] alist;
  delete [] blist;
  delete [] cv;
  delete [] num_illegal;
}
