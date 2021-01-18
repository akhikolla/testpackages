#include "resultsub1.h"

ResultSub1::ResultSub1()
{
  admissible = false;
  useCurtailment = false;
  this->name = new std::string();
  curtailmentResults = new std::map<int, Curtailment_SubD1>();
}

ResultSub1::ResultSub1(int n, int r, int s, int n1, int r1, double alpha, double beta, double petP0, double enP0, int iD,
                double pc0, double pt0, double pc1, double pt1)
{
  this->n = n;
  this->r = r;
  this->s = s;
  this->n1 = n1;
  this->r1 = r1;
  this->alpha = alpha;
  this->beta = beta;
  this->petP0 = petP0;
  this->enP0 = enP0;
  this->admissible = false;
  this->iD = iD;
  this->pc0 = pc0;
  this->pt0 = pt0;
  this->pc1 = pc1;
  this->pt1 = pt1;
  this->name = new std::string();
  this->useCurtailment = false;
  curtailmentResults = new std::map<int, Curtailment_SubD1>();
}

ResultSub1::~ResultSub1()
{
  delete name;
  delete curtailmentResults;
}

bool ResultSub1::getAdmissible(){return admissible;}

void ResultSub1::setAdmissible(double start, double stop)
{
  this->admissibleStart = start;
  this->admissibleStop = stop;
  this->admissible = true;
}

void ResultSub1::setAdmissible(double start, double stop, std::string name)
{
  setAdmissible(start, stop);
  *(this->name) = name;
}

int ResultSub1::getN(){return this->n;}
void ResultSub1::setN(int n){this->n = n;}

int ResultSub1::getR(){return this->r;}
void ResultSub1::setR(int r){this->r = r;}

int ResultSub1::getS(){return this->s;}
void ResultSub1::setS(int s){this->s = s;}

int ResultSub1::getN1(){return this->n1;}
void ResultSub1::setN1(int n1){this->n1 = n1;}

int ResultSub1::getR1(){return this->r1;}
void ResultSub1::setR1(int r1){this->r1 = r1;}

double ResultSub1::getEnP0(){return this->enP0;}
void ResultSub1::setEnP0(double enP0){this->enP0 = enP0;}

double ResultSub1::getPetP0(){return petP0;}
void ResultSub1::setPetP0(double petP0){this->petP0 = petP0;}

void ResultSub1::addCurtailmentResult(Curtailment_SubD1 curResult)
{    
    curtailmentResults->insert(std::pair<int, Curtailment_SubD1>((int)(curResult.cut*100 +0.5), curResult));
    this->useCurtailment = true;
}

std::map<int, ResultSub1::Curtailment_SubD1>* ResultSub1::getCurtailmentResults(){return curtailmentResults;}

void ResultSub1::setUseCurtailment(bool useCurtailment){this->useCurtailment = useCurtailment;}

void ResultSub1::setCut(int cut){ this->cut = cut;}

int ResultSub1::getID(){return this->iD;}

SEXP ResultSub1::get_R_Representation()
{
  using namespace Rcpp;
  List rRepres;
  List curResults = List::create();  
  DataFrame stoppingRules;
  
  DataFrame solution = DataFrame::create(_["ID"] = NumericVector::create(iD),
                      _["r1"] = NumericVector::create(r1),
                      _["n1"] = NumericVector::create(n1),
                      _["r"] = NumericVector::create(r),
                      _["s"] = NumericVector::create(s),
                      _["n"] = NumericVector::create(n),
                      _["pc0"] = NumericVector::create(pc0),
                      _["pt0"] = NumericVector::create(pt0),
                      _["pc1"] = NumericVector::create(pc1),
                      _["pt1"] = NumericVector::create(pt1),
                      _["enP0"] = NumericVector::create(enP0),
                      _["petP0"] = NumericVector::create(petP0),
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
    std::map<int, Curtailment_SubD1>::iterator it;
    Curtailment_SubD1 cur;
    
    std::vector<StoppingRule_SubD1>::iterator stopRuleIt;
    StoppingRule_SubD1 sr;
    
    for(it = curtailmentResults->begin(); it != curtailmentResults->end(); it++)
    {
      cur = it->second;
      
      DataFrame tmp = DataFrame::create( _["Cut"] = NumericVector::create(cur.cut),
                                           _["En_SC"] = NumericVector::create(cur.en_sc),
                                           _["Pet_SC"] = NumericVector::create(cur.pet_sc),
                                           _["Type_1_Errorrate"] = NumericVector::create(cur.type1_errorRate),
                                           _["Type_2_Errorrate"] = NumericVector::create(cur.type2_errorRate));
      
      //get stopping rules for given curtailment
      std::vector<float> tmpT;
      std::vector<float> tmpU;
      std::vector<float> tmpEnrolled;
      
      for(stopRuleIt = cur.stoppingRulesNSC->begin(); stopRuleIt != cur.stoppingRulesNSC->end(); stopRuleIt++)
      {
        sr = *stopRuleIt;
        tmpT.push_back(sr.t_int);
        tmpU.push_back(sr.u_int);
        tmpEnrolled.push_back(sr.enrolled_int);
      }
      
      DataFrame t = DataFrame::create(_["Needed_responses_ep1"] = tmpT);
      DataFrame u = DataFrame::create(_["Needed_responses_ep2"] = tmpU);
      DataFrame enrolled = DataFrame::create(_["Enrolled_patients"] = tmpEnrolled);
      stoppingRules = DataFrame::create(_["Needed_responses_ep1"] = t, _["Needed_responses_ep2"] = u, _["Enrolled_patients"] = enrolled);
      
      List curResult = List::create(_["Details"] = tmp, _["Stoppingrules"] = stoppingRules);
      curResults.push_back(curResult);
    }
  }
  
  rRepres.push_back(curResults, "Curtailment");

  return rRepres;
}
