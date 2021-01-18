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
#ifndef DYNSBM_EM_H
#define DYNSBM_EM_H
#include<DynSBM.h>
#include<iostream>
namespace dynsbm{
  template<class TDynSBM, typename Ytype>
  class EM{
  private:
    TDynSBM _model;
  public:
    EM(int T, int N, int Q, const Rcpp::IntegerMatrix& present, bool isdirected = false, bool withselfloop = false)
      : _model(T,N,Q,present,isdirected,withselfloop) {}
    ~EM(){};
    const TDynSBM& getModel() const{
      return _model;
    }
    void initialize(const std::vector<int>& clustering, Ytype*** const Y, bool frozen=false){ 
      _model.initTau(clustering);
      if(frozen)
	_model.updateFrozenTheta(Y);
      else
	_model.updateTheta(Y);	
      _model.initNotinformativeStationary();
      _model.initNotinformativeTrans();
    }
    int run(Ytype*** const Y, int nbit, int nbitFP, bool frozen){
      double prevlogl = _model.completedLoglikelihood(Y);
      //---- estimation
      int it = 0, nbiteff = 0;
      while(it<nbit){
	int itfp = 0;
	double prevloglfp = prevlogl;
	while (itfp<nbitFP){
	  _model.updateTau(Y);
	  if (itfp%3==0){ // saving time
	    double newloglfp = _model.completedLoglikelihood(Y);
	    if(fabs((prevloglfp-newloglfp)/prevloglfp)<1e-4){
	      itfp = nbitFP;
	    } else{
	      prevloglfp = newloglfp;
	      itfp = itfp+1;
	    }
	  } else
	    itfp = itfp+1;
#ifdef DEBUG
	  std::cerr<<"After EStep: "<<_model.completedLoglikelihood(Y)<<std::endl;
#endif
	}
	_model.updateTrans();
#ifdef DEBUG
	std::cerr<<"After MStep on trans: "<<_model.completedLoglikelihood(Y)<<std::endl;
#endif
	_model.updateStationary();
	if(frozen)
	  _model.updateFrozenTheta(Y);
	else
	  _model.updateTheta(Y);	  
#ifdef DEBUG
	std::cerr<<"After MStep on theta: "<<_model.completedLoglikelihood(Y)<<std::endl;
#endif
	double newlogl = _model.completedLoglikelihood(Y);
#ifdef DEBUG
	std::cerr<<"Testing the likelihood decrease: "<<prevlogl<<" -> "<<newlogl<<std::endl;
#endif
	nbiteff++;
	if(fabs((prevlogl-newlogl)/prevlogl)<1e-4){
#ifdef DEBUG
	  std::cerr<<"Stopping: criteria is reached"<<std::endl;
#endif
	  it=nbit;
	}
	if(prevlogl>newlogl){
#ifdef DEBUG
	  std::cerr<<"Stopping: increasing logl"<<std::endl;
#endif
	  it=nbit;
	}
	prevlogl = newlogl;
	it = it+1;
      }
      return(nbiteff);
    }
  };
}
#endif
