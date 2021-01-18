#ifndef Data_H
#define Data_H


#include <RcppArmadillo.h>
#define NDEBUG 1
#include <vector>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Data{
public:
  int m_nobs;
  vector< vector< Col<double> > > m_yi;
  vector< Col<int> > m_nbtimeobs;
  Col<int> nbseq;

  Data(){};
  ~Data(){};
  Data(const S4 & obj){
    this->m_nobs = obj.slot("nobs");
    this->nbseq = as< Col<int> >(obj.slot("nbseq"));
    this->m_nbtimeobs = as< vector< Col<int> >  >(obj.slot("nbtimeobs"));
    this->m_yi = as< vector< vector< Col<double> > > >(obj.slot("yi"));
  }
};
#endif
