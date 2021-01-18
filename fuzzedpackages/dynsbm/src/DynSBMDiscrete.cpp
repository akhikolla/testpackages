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
#include<DynSBMDiscrete.h>
namespace dynsbm{
  void DynSBMDiscrete::updateTheta(int*** const Y){// M-step
    for(int t=0;t<_t;t++)
      for(int q=0;q<_q;q++)
        for(int l=0;l<_q;l++)
          for(int k=0;k<_k;k++) _multinomprobaql[t][q][l][k] = 0.;
    
    DynSBMDiscreteAddEventFunctor addEventFunctor(*this);
    updateThetaCore<DynSBMDiscreteAddEventFunctor>(Y, addEventFunctor);
    
    for(int t=0;t<_t;t++){// symmetrize+normalize
      for(int q=(_isdirected?0:1);q<_q;q++){
        for(int l=0;l<q;l++){
	  double summultinomprobaql = 0.;
	  for(int k=0;k<_k;k++) summultinomprobaql += _multinomprobaql[t][q][l][k];
	  if (summultinomprobaql>0)
	    for(int k=0;k<_k;k++){
	      _multinomprobaql[t][q][l][k] = _multinomprobaql[t][q][l][k]/summultinomprobaql;
	      if(!_isdirected) _multinomprobaql[t][l][q][k] = _multinomprobaql[t][q][l][k];
	    }
	}        	
	if(_isdirected)
	  for(int l=q+1;l<_q;l++){
	    double summultinomprobaql = 0.;
	    for(int k=0;k<_k;k++) summultinomprobaql += _multinomprobaql[t][q][l][k];
	    if (summultinomprobaql>0)
	      for(int k=0;k<_k;k++){
		_multinomprobaql[t][q][l][k] = _multinomprobaql[t][q][l][k]/summultinomprobaql;
	      }
	  }      
      }
    }
    for(int q=0;q<_q;q++){// symmetrize+normalize
      double summultinomprobaqq = 0.;
      for(int k=0;k<_k;k++) summultinomprobaqq += _multinomprobaql[0][q][q][k];
      if(summultinomprobaqq>0)
	for(int k=0;k<_k;k++) _multinomprobaql[0][q][q][k] = _multinomprobaql[0][q][q][k]/summultinomprobaqq;    
      for(int t=1;t<_t;t++)
        for(int k=0;k<_k;k++) _multinomprobaql[t][q][q][k] =  _multinomprobaql[0][q][q][k];      
    }
    correctMultinomproba();
  }
  
  void DynSBMDiscrete::updateFrozenTheta(int*** const Y){// M-step
    for(int t=0;t<_t;t++)
      for(int q=0;q<_q;q++)
        for(int l=0;l<_q;l++)
          for(int k=0;k<_k;k++) _multinomprobaql[t][q][l][k] = 0.;
    
    DynSBMDiscreteAddEventFunctor addEventFunctor(*this);
    updateFrozenThetaCore<DynSBMDiscreteAddEventFunctor>(Y, addEventFunctor);
    
    // symmetrize+normalize
    for(int q=(_isdirected?0:1);q<_q;q++){
      for(int l=0;l<q;l++){
	double summultinomprobaql = 0.;
	for(int k=0;k<_k;k++) summultinomprobaql += _multinomprobaql[0][q][l][k];
	if (summultinomprobaql>0)
	  for(int k=0;k<_k;k++){
	    _multinomprobaql[0][q][l][k] = _multinomprobaql[0][q][l][k]/summultinomprobaql;
	    if(!_isdirected) _multinomprobaql[0][l][q][k] = _multinomprobaql[0][q][l][k];
	  }
      }        	
      if(_isdirected)
	for(int l=q+1;l<_q;l++){
	  double summultinomprobaql = 0.;
	  for(int k=0;k<_k;k++) summultinomprobaql += _multinomprobaql[0][q][l][k];
	  if (summultinomprobaql>0)
	    for(int k=0;k<_k;k++){
	      _multinomprobaql[0][q][l][k] = _multinomprobaql[0][q][l][k]/summultinomprobaql;
	    }
	}      
    }
    for(int q=0;q<_q;q++){// symmetrize+normalize
      double summultinomprobaqq = 0.;
      for(int k=0;k<_k;k++) summultinomprobaqq += _multinomprobaql[0][q][q][k];
      if(summultinomprobaqq>0)
	for(int k=0;k<_k;k++) _multinomprobaql[0][q][q][k] = _multinomprobaql[0][q][q][k]/summultinomprobaqq;
    }
    // note the constraint multinomprobaql[t,q,l,.] = multinomprobaql[0,q,l,.]
    for(int t=1;t<_t;t++){
      for(int q=0;q<_q;q++)
	for(int l=0;l<_q;l++)
	  for(int k=0;k<_k;k++) _multinomprobaql[t][q][l][k] =  _multinomprobaql[0][q][l][k];
    }
    correctMultinomproba();
  }
  
  void DynSBMDiscrete::correctMultinomproba(){ // numerical issue : avoid too small value for multinomproba
    for(int t=0;t<_t;t++){
        for(int q=0;q<_q;q++){
            for(int l=0;l<_q;l++){
                for(int k=0;k<_k;k++){
                  if(_multinomprobaql[t][q][l][k]<precision)
                    _multinomprobaql[t][q][l][k] = precision;
                  if(_multinomprobaql[t][q][l][k]>(1-precision)){
                        _multinomprobaql[t][q][l][k] = 1-precision;
                  }
		  // trick: store log(_multinomprobaql[t][q][l][k])
		   _multinomprobaql[t][q][l][k] = log(_multinomprobaql[t][q][l][k]);
                }
            }
        }
    }
  }
}
