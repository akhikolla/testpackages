#ifndef ParamSpecificLaw_H
#define ParamSpecificLaw_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class ParamSpecificLaw{
public:
  Col<double> m_eps, m_a, m_b;
  
  ParamSpecificLaw(){};
  ~ParamSpecificLaw(){};

  ParamSpecificLaw(const List mysep){
    this->m_eps = as< Col<double> >(mysep[0]);
    this->m_a = as< Col<double> >(mysep[1]);
    this->m_b = as< Col<double> >(mysep[2]);
  };


  Col<double> dlogspecific(const Col<double> & obs, const int & state){
    Col<double> out = obs * 0;
      uvec who = find(obs != 0);
      out.elem(who) =  log(1- m_eps(state))  + m_a(state) * log(m_b(state)) + (m_a(state) - 1) * log(obs.elem(who)) - m_b(state) * obs.elem(who) - lgamma(m_a(state));
      who = find(obs == 0);
      out.elem(who) = log( m_eps(state)) * ones(who.n_rows);
    return out;
  }

  Mat<double> dlogspecificAll(const Col<double> & obs, const int & statemax){
    Mat<double> out = zeros(obs.n_rows, statemax);
    for (int h=0; h<statemax; h++)
      out.col(h) = dlogspecific(obs, h);
    return out;
  }


  Mat<double> dspecificProbaCond(const Col<double> & obs, const int & state){
    Mat<double> out = zeros(obs.n_rows, 2);
    Col<double> tmp1 =ones(obs.n_rows);
    Col<double> tmp2 =zeros(obs.n_rows);
      tmp1 =zeros(obs.n_rows);
      uvec who = find(obs != 0);
      tmp1.elem(who) =  ones(who.n_rows);
      who = find(obs == 0);
      tmp2.elem(who) =  ones(who.n_rows);
      out.col(0) = tmp1;
      out.col(1) = tmp2;
    
    return out;
  }


};
#endif
