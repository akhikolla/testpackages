#include "result.h"

Result::Result()
{
  admissible = false;
  this->stoppingRulesNSC = new std::vector<StoppingRule>();
  this->curtailmentResults = new std::map<int, Curtailment>();
  this->name = new std::string();
  this->useCurtailment = false;
}

Result::Result(int n, int r, int n1, int r1, double alpha, double beta, double petP1, double petP0, double enP1, double enP0, int iD, double p0, double p1)
{
  this->n = n;
  this->r = r;
  this->n1 = n1;
  this->r1 = r1;
  this->alpha = alpha;
  this->beta = beta;
  this->petP0 = petP0;
  this->petP1 = petP1;
  this->enP0 = enP0;
  this->enP1 = enP1;
  this->admissible = false;
  this->stoppingRulesNSC = new std::vector<StoppingRule>();
  this->curtailmentResults = new std::map<int, Curtailment>();
  this->iD = iD;
  this->name = new std::string();
  this->p0 = p0;
  this->p1 = p1;
}

Result::~Result()
{
  delete stoppingRulesNSC;
  delete curtailmentResults;
  delete name;
}

void Result::setAdmissible(double start, double stop)
{
  this->admissibleStart = start;
  this->admissibleStop = stop;
  this->admissible = true;
}

void Result::setAdmissible(double start, double stop, std::string name)
{
  setAdmissible(start, stop);
  *(this->name) = name;
}

void Result::addCurtailmentResult(Curtailment curResult)
{    
    curtailmentResults->insert(std::pair<int, Curtailment>((int)(curResult.cut*100 +0.5), curResult));
    this->useCurtailment = true; 
}

void Result::setUseCurtailment(bool useCurtailment){this->useCurtailment = useCurtailment;}

void Result::setCut(int cut){ this->cut = cut;}

SEXP Result::get_R_Representation()
{
  using namespace Rcpp;
  List rRepres;
  List curResults = List::create();
  DataFrame stoppingRules;
  DataFrame solution = DataFrame::create(_["ID"] = NumericVector::create(iD),
                      _["r1"] = NumericVector::create(r1),
                      _["n1"] = NumericVector::create(n1),
                      _["r"] = NumericVector::create(r),
                      _["n"] = NumericVector::create(n),
                      _["p0"] = NumericVector::create(p0),
                      _["p1"] = NumericVector::create(p1),
                      _["enP0"] = NumericVector::create(enP0),
                      _["enP1"] = NumericVector::create(enP1),
                      _["petP0"] = NumericVector::create(petP0),
                      _["petP1"] = NumericVector::create(petP1),
                      _["Alpha"] = NumericVector::create(alpha),
                      _["Beta"] = NumericVector::create(beta),
                      _["Admissible"] = NumericVector::create(admissible),
                      _["Admiss_Start"] = NumericVector::create(admissible ? admissibleStart : NA_REAL),
                      _["Admiss_End"] = NumericVector::create(admissible ? admissibleStop : NA_REAL),
                      _["Type"] = CharacterVector::create(*name));

  rRepres.push_back(solution, "Solution");
  // If curtailment is used return all available information regarding (non-)stochastic curtailment.
  if(useCurtailment)
  {
    std::map<int, Result::Curtailment>::iterator it = curtailmentResults->begin();
    while( it != curtailmentResults->end())      
    {        
        Result::Curtailment cur = it->second;
        DataFrame tmp = DataFrame::create( _["Cut"] = NumericVector::create(cur.cut),
                                           _["En_SC"] = NumericVector::create(cur.en_sc),
                                           _["Pet_SC"] = NumericVector::create(cur.pet_sc),
                                           _["Type_1_Errorrate"] = NumericVector::create(cur.type1_errorRate),
                                           _["Type_2_Errorrate"] = NumericVector::create(cur.type2_errorRate),                                          
                                           _["En_lower"] = NumericVector::create(cur.en_lower), 
                                           _["En_upper"] = NumericVector::create(cur.en_upper),
                                           _["Pet_lower"] = NumericVector::create(cur.pet_lower),
                                           _["Pet_upper"] = NumericVector::create(cur.pet_upper),
                                           _["alpha_lower"] = NumericVector::create(cur.alpha_lower),
                                           _["alpha_upper"] = NumericVector::create(cur.alpha_upper),                                           
                                           _["beta_lower"] = NumericVector::create(cur.beta_lower),
                                           _["beta_upper"] = NumericVector::create(cur.beta_upper));

        std::vector<float> tmpRi;
        std::vector<float> tmpK;
        for(unsigned int i = 0; i < cur.stoppingRulesNSC->size(); i ++ ) 
        {
          tmpRi.push_back(cur.stoppingRulesNSC->at(i)[0]);
          tmpK.push_back(cur.stoppingRulesNSC->at(i)[1]);
        }
        DataFrame ri = DataFrame::create(_["Needed_responses"] = tmpRi);
        DataFrame k = DataFrame::create(_["Observed_patients"] = tmpK);
        stoppingRules = DataFrame::create(_["Needed_responses"] = ri, _["Observed_patients"] = k);
        List curResult = List::create(_["Details"] = tmp, _["Stoppingrules"] = stoppingRules);
        curResults.push_back(curResult);
        it++;
    }
  }

  rRepres.push_back(curResults, "Curtailment");

  return rRepres;

}

bool Result::getAdmissible(){return admissible;}

int Result::getN(){return n;}

int Result::getR(){return r;}

int Result::getN1(){return n1;}

int Result::getR1(){return r1;}

double Result::getEnP0(){return enP0;}

double Result::getEnP1(){return enP1;}

double Result::getPetP0(){return petP0;}

std::map<int, Result::Curtailment>* Result::getCurtailmentResults(){return curtailmentResults;}

