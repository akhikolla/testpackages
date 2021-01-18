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
#include<DynSBMGaussian.h>
#include<iostream>
namespace dynsbm{
  void DynSBMGaussian::updateTheta(double*** const Y){// M-step
    for(int t=0;t<_t;t++){
      _sigma[t] = 0.;
      for(int q=0;q<_q;q++)
	for(int l=0;l<_q;l++)
	  _muql[t][q][l] = 0;
    }
    // note the constraint _muql[t,q,q] = _muql[,q,q]
    double*** sumql;
    allocate3D<double>(sumql,_t,_q,_q);
    
    DynSBMGaussianAddEventFunctor addEventFunctor(*this,sumql);
    updateThetaCore<DynSBMGaussianAddEventFunctor>(Y, addEventFunctor);
    
    for(int t=0;t<_t;t++){// symmetrize+normalize 
      for(int q=(_isdirected?0:1);q<_q;q++){
	for(int l=0;l<q;l++){
	  if (sumql[t][q][l]>0){
	    _muql[t][q][l] = _muql[t][q][l] / sumql[t][q][l];
	    if(!_isdirected) _muql[t][l][q] = _muql[t][q][l];
	  }
	} 	
	if(_isdirected)
	  for(int l=q+1;l<_q;l++)
	    if (sumql[t][q][l]>0)
	      _muql[t][q][l] = _muql[t][q][l] / sumql[t][q][l];	
      }
    }
    for(int q=0;q<_q;q++){// symmetrize+normalize
      // note the constraint muql[t,q,q] = muql[,q,q]
      if (sumql[0][q][q]>0)
	_muql[0][q][q] = _muql[0][q][q] / sumql[0][q][q];
      for(int t=1;t<_t;t++)
	_muql[t][q][q] = _muql[0][q][q];
    }
    // _sigma - homoscedastic at each time
    for(int t=0;t<_t;t++){
      double sumtaumarginalt = 0.;
      for(int i=1;i<_n;i++){
	if (ispresent(t,i)){
	  // j<i
	  for(int j=0;j<i;j++){
	    if (ispresent(t,j))
	      for(int q=0;q<_q;q++){
		for(int l=0;l<=q;l++){
		  if (q!=l){
		    if(Y[t][i][j]>0){
		      _sigma[t] += tauMarginal(t,i,q)*tauMarginal(t,j,l)*(Y[t][i][j]-_muql[t][q][l])*(Y[t][i][j]-_muql[t][q][l]) +
			tauMarginal(t,i,l)*tauMarginal(t,j,q)*(Y[t][i][j]-_muql[t][q][l])*(Y[t][i][j]-_muql[t][q][l]);
		      sumtaumarginalt += tauMarginal(t,i,q)*tauMarginal(t,j,l) + tauMarginal(t,i,l)*tauMarginal(t,j,q);
		      if(_isdirected){
			_sigma[t] += tauMarginal(t,j,q)*tauMarginal(t,i,l)*(Y[t][j][i]-_muql[t][q][l])*(Y[t][j][i]-_muql[t][q][l]) +
			  tauMarginal(t,j,l)*tauMarginal(t,i,q)*(Y[t][j][i]-_muql[t][l][q])*(Y[t][j][i]-_muql[t][l][q]);
			sumtaumarginalt += tauMarginal(t,j,q)*tauMarginal(t,i,l) + tauMarginal(t,j,l)*tauMarginal(t,i,q);
		      }
		    }
		  } else{
		    if(Y[t][i][j]>0){
		      _sigma[t] += tauMarginal(t,i,q)*tauMarginal(t,j,q)*(Y[t][i][j]-_muql[t][q][q])*(Y[t][i][j]-_muql[t][q][q]);
		      sumtaumarginalt += tauMarginal(t,i,q)*tauMarginal(t,j,q);
		      if(_isdirected){
			_sigma[t] += tauMarginal(t,i,q)*tauMarginal(t,j,q)*(Y[t][j][i]-_muql[t][q][q])*(Y[t][j][i]-_muql[t][q][q]);
			sumtaumarginalt += tauMarginal(t,i,q)*tauMarginal(t,j,q);
		      }
		    }
		  }
		}
	      }
	  }
	  // j==i considered only if selfloop allowed
	  if (_withselfloop)
	    for(int q=0;q<_q;q++)
	      if(Y[t][i][i]>0){		  
		_sigma[t] += tauMarginal(t,i,q)*(Y[t][i][i]-_muql[t][q][q])*(Y[t][i][i]-_muql[t][q][q]);
		sumtaumarginalt += tauMarginal(t,i,q);
	      }
	}	  
      }
      _sigma[t] = sqrt(_sigma[t]/sumtaumarginalt);
    }
    deallocate3D<double>(sumql,_t,_q,_q);
  }
  
  void DynSBMGaussian::updateFrozenTheta(double*** const Y){// M-step
    for(int t=0;t<_t;t++){
      _sigma[t] = 0.;
      for(int q=0;q<_q;q++)
	for(int l=0;l<_q;l++)
	  _muql[t][q][l] = 0;
    }    
    double*** sumql;
    allocate3D<double>(sumql,_t,_q,_q);
    
    DynSBMGaussianAddEventFunctor addEventFunctor(*this,sumql);
    updateFrozenThetaCore<DynSBMGaussianAddEventFunctor>(Y, addEventFunctor);
    
    // symmetrize+normalize 
    for(int q=(_isdirected?0:1);q<_q;q++){
      for(int l=0;l<q;l++){
	if (sumql[0][q][l]>0){
	  _muql[0][q][l] = _muql[0][q][l] / sumql[0][q][l];
	  if(!_isdirected) _muql[0][l][q] = _muql[0][q][l];
	}
      } 	
      if(_isdirected)
	for(int l=q+1;l<_q;l++)
	  if (sumql[0][q][l]>0)
	    _muql[0][q][l] = _muql[0][q][l] / sumql[0][q][l];	
    }   
    for(int q=0;q<_q;q++){// symmetrize+normalize
      if (sumql[0][q][q]>0)
	_muql[0][q][q] = _muql[0][q][q] / sumql[0][q][q];
    }
    // note the constraint muql[t,q,l] = muql[0,q,l]
    for(int t=1;t<_t;t++)
      for(int q=0;q<_q;q++)
	for(int l=0;l<_q;l++)
	  _muql[t][q][l] = _muql[0][q][l];    
    // _sigma - homoscedastic at each time
    for(int t=0;t<_t;t++){
      double sumtaumarginalt = 0.;
      for(int i=1;i<_n;i++){
	if (ispresent(t,i)){
	  // j<i
	  for(int j=0;j<i;j++){
	    if (ispresent(t,j))
	      for(int q=0;q<_q;q++){
		for(int l=0;l<=q;l++){
		  if (q!=l){
		    if(Y[t][i][j]>0){
		      _sigma[0] += tauMarginal(t,i,q)*tauMarginal(t,j,l)*(Y[t][i][j]-_muql[t][q][l])*(Y[t][i][j]-_muql[t][q][l]) +
			tauMarginal(t,i,l)*tauMarginal(t,j,q)*(Y[t][i][j]-_muql[t][q][l])*(Y[t][i][j]-_muql[t][q][l]);
		      sumtaumarginalt += tauMarginal(t,i,q)*tauMarginal(t,j,l) + tauMarginal(t,i,l)*tauMarginal(t,j,q);
		      if(_isdirected){
			_sigma[0] += tauMarginal(t,j,q)*tauMarginal(t,i,l)*(Y[t][j][i]-_muql[t][q][l])*(Y[t][j][i]-_muql[t][q][l]) +
			  tauMarginal(t,j,l)*tauMarginal(t,i,q)*(Y[t][j][i]-_muql[t][l][q])*(Y[t][j][i]-_muql[t][l][q]);
			sumtaumarginalt += tauMarginal(t,j,q)*tauMarginal(t,i,l) + tauMarginal(t,j,l)*tauMarginal(t,i,q);
		      }
		    }
		  } else{
		    if(Y[t][i][j]>0){
		      _sigma[0] += tauMarginal(t,i,q)*tauMarginal(t,j,q)*(Y[t][i][j]-_muql[t][q][q])*(Y[t][i][j]-_muql[t][q][q]);
		      sumtaumarginalt += tauMarginal(t,i,q)*tauMarginal(t,j,q);
		      if(_isdirected){
			_sigma[0] += tauMarginal(t,i,q)*tauMarginal(t,j,q)*(Y[t][j][i]-_muql[t][q][q])*(Y[t][j][i]-_muql[t][q][q]);
			sumtaumarginalt += tauMarginal(t,i,q)*tauMarginal(t,j,q);
		      }
		    }
		  }
		}
	      }
	  }
	  // j==i considered only if selfloop allowed
	  if (_withselfloop)
	    for(int q=0;q<_q;q++)
	      if(Y[t][i][i]>0){		  
		_sigma[0] += tauMarginal(t,i,q)*(Y[t][i][i]-_muql[t][q][q])*(Y[t][i][i]-_muql[t][q][q]);
		sumtaumarginalt += tauMarginal(t,i,q);
	      }
	}	  
      }
      _sigma[0] = sqrt(_sigma[0]/sumtaumarginalt);
      // note the constraint sigma[t]=sigma[0]
      for(int t=1;t<_t;t++) _sigma[t] = _sigma[0];       
    }
    deallocate3D<double>(sumql,_t,_q,_q); 
  }
}
