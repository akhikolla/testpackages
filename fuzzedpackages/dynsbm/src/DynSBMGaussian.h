/*
    This file is part of dynsbm.

    dysbm is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    dynsbm is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with dynsbm.  If not, see <http://www.gnu.org/licenses/>
 */
#ifndef DYNSBM_DYNSBMGAUSSIAN_H
#define DYNSBM_DYNSBMGAUSSIAN_H
#include<DynSBM.h>
#include <R.h>
#include <Rmath.h>
namespace dynsbm{
  class DynSBMGaussian
    : public DynSBM<double>{
  protected:
    double*** _muql;
    double* _sigma;
  public:
    DynSBMGaussian(int T, int N, int Q, const Rcpp::IntegerMatrix & present, bool isdirected = false, bool withselfloop = false)
      : DynSBM<double>(T,N,Q,present,isdirected,withselfloop) {
      allocate3D(_muql,_t,_q,_q);
      _sigma = new double[_t];
    }
    ~DynSBMGaussian(){
      deallocate3D(_muql,_t,_q,_q);
      delete[] _sigma;
    }
    double getMu(int t, int q, int l) const{
      return(_muql[t][q][l]);
    }
    double getSigma(int t) const{
      return(_sigma[t]);
    }
    virtual double logDensity(int t, int q, int l, double y) const{
      if(y>0.){ // testing y==0 for real values
          int give_log = 1;
          return(_1minusbetaql[t][q][l] // which is actually log(1-_betaql[t][q][l]))
		 + Rf_dnorm4(y, _muql[t][q][l], _sigma[t], give_log));
      } else{
          return(_betaql[t][q][l]); // which is actually log(_betaql[t][q][l]))
      }
    }
    virtual void updateTheta(double*** const Y);
    virtual void updateFrozenTheta(double*** const Y);
    friend class DynSBMGaussianAddEventFunctor;
  };

  
  class DynSBMGaussianAddEventFunctor{
    DynSBMGaussian& _dynsbm;
    double*** _sumql;
  public:
    DynSBMGaussianAddEventFunctor(DynSBMGaussian& dynsbm, double*** sumql)
      : _dynsbm(dynsbm), _sumql(sumql) {}
    void operator()(double proba, double y, int t, int q, int l){
      _dynsbm._muql[t][q][l] += proba*y;
      _sumql[t][q][l] += proba;
    }
  };
}
#endif
