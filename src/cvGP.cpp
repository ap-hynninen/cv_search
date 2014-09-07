#include "../gpc/include/gp.h"

// Define class identifiers
const int cvGeneID=GPUserID;
const int cvGPID=GPUserID+1;
const int cvPopulationID=GPUserID+2;

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
  double evaluate ();

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

// This function evaluates the fitness of a genetic tree.  We have the
// freedom to define this function in any way we like.  
double cvGene::evaluate ()
{
  double returnValue=0.0;

  /*
  switch (node->value ())
    {
    case FUNCTION1: 
      returnValue=NthMyChild(0)->evaluate ();
      break;

    case FUNCTION2: 
      returnValue=NthMyChild(0)->evaluate ();
      break;
      
    case TERMINAL1:
      returnValue=0.0;
      break;
      
    case TERMINAL2:
      returnValue=0.0;
      break;

    default: 
      GPExitSystem ("cvGene::evaluate", "Undefined node value");
    }
  */

  return returnValue;
}



// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void cvGP::evaluate ()
{
  // Evaluate main tree
  stdFitness=NthMyGene (0)->evaluate ();
}
