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
#ifndef DYNSBM_DYNSBMBINARY_H
#define DYNSBM_DYNSBMBINARY_H
#include<DynSBM.h>
namespace dynsbm{
  class DynSBMBinary 
    : public DynSBM<int>{
  public:
    DynSBMBinary(int T, int N, int Q, const Rcpp::IntegerMatrix & present, bool isdirected = false, bool withselfloop = false)
      : DynSBM<int>(T,N,Q,present,isdirected,withselfloop) {}
    ~DynSBMBinary(){}
    virtual double logDensity(int t, int q, int l, int y) const{
      if(y==0){
	return(_betaql[t][q][l]); // trick: which is actually log(_betaql[t][q][l]))
      } else{
	return(_1minusbetaql[t][q][l]); // trick: which is actually log(1-_betaql[t][q][l]))
      }
    }
    virtual void updateTheta(int*** const Y);
    virtual void updateFrozenTheta(int*** const Y);
  };
  
  class DynSBMBinaryAddEventFunctor{
    DynSBMBinary& _dynsbm;
  public:
    DynSBMBinaryAddEventFunctor(DynSBMBinary& dynsbm)
      : _dynsbm(dynsbm) {}
    void operator()(double proba, int y, int t, int q, int l){}
  };
}
#endif
